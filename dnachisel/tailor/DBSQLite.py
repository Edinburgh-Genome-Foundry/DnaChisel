'''
Created on Nov 1, 2012

@author: jcg
'''

from sqlite3 import connect,Row
from time import strftime
from math import sqrt
from random import choice
from uuid import uuid4
from subprocess import check_output
from socket import gethostbyname,gethostname
from .tools import pick_random

class DBSQLite():
    '''
    Constructor
    '''
    def __init__(self, dbfile, designMethod, initialize=True, seedSequence=None):
        self.dbfile = dbfile
        self.designMethod = designMethod
        self.seedSequence = seedSequence
        self.seedId = str(uuid4().int)
        
        #SQL queries buffers
        self.des_solutions = {}
        self.des_solutions_sql = []
        self.gen_solutions_sql = []
        self.gen_solutions_id = {}        
        
        #Connect database
        self.con = connect(self.dbfile + ".sqlite")
        self.con.row_factory = Row
        self.con.isolation_level = None
        self.cur = self.con.cursor()
        
        #Initiate DB
        if initialize:
            self.DBInit()
        
        #Register worker
        self.worker_id = str(uuid4().int)        
        self.registerWorker()
        

    def DBInit(self):
        '''
        Initialize database
        input: dbfile - file to save DB        
        input: designMethod - a class Design with information about features name/type and design wanted  
        returns: Nothing
        '''
        pass        
        
        #Design Dynamic tables
        features_fields = ''.join([feature+"Level TEXT, " for feature in self.designMethod.featuresList])
        
        table_ds = "create table desired_solution(des_solution_id TEXT PRIMARY KEY,"+features_fields+\
                   "status TEXT,worker_id TEXT,start_time TEXT,FOREIGN KEY(worker_id) REFERENCES worker(worker_id));"
        
        features_values_fields = ''.join([feature+" "+self.designMethod.features[feature]['type']+", " for feature in self.designMethod.featuresList])
        features_level_position_fields = ''.join([feature+"Position "+self.designMethod.features[feature]['type']+", " for feature in self.designMethod.featuresList])
                           
        table_gs = "create table generated_solution(generated_solution_id TEXT PRIMARY KEY, des_solution_id TEXT, sequence TEXT,"+features_values_fields+features_fields+features_level_position_fields+\
                   "worker_id TEXT, FOREIGN KEY(worker_id) REFERENCES worker(worker_id));"

        
        #Create Tables
        self.cur.executescript("""
                        PRAGMA writable_schema = 1;
                        delete from sqlite_master where type = 'table';
                        PRAGMA writable_schema = 0;

                        VACUUM;
                        PRAGMA INTEGRITY_CHECK;
                        
                        create table worker(worker_id TEXT PRIMARY KEY, hostname TEXT, ip TEXT, time_start INTEGER, time_finish INTEGER);""" +\
                        table_ds + "\n" + \
                        table_gs + "\n")
                        
        #Populate Tables
                        
        #Desired Solutions DB
        n_features = self.designMethod.n_features
        all_comb = [tuple(d_sol.split('.')) + (d_sol,) for d_sol in self.designMethod.listDesigns]
        features_levels_fields = ''.join([feature+"Level, " for feature in self.designMethod.featuresList])
                 
        sql = "insert into desired_solution("+features_levels_fields+"des_solution_id,status, worker_id, start_time) \
                                     values("+"?,"*(n_features+1)+"'WAITING',NULL,NULL);"
        
        self.cur.executemany(sql,all_comb)
        
    
    def registerWorker(self):
        start_time = strftime("%Y-%m-%d %H:%M:%S %Z")
        hostname = check_output("hostname").rstrip()
        ip = gethostbyname(gethostname()).rstrip()
        
        self.cur.execute("insert into worker(worker_id, hostname, ip, time_start, time_finish) values (?,?,?,?,NULL);", 
                        (self.worker_id, hostname, ip, start_time) )
                
        self.cur.execute("select * from desired_solution")
        
        for row in self.cur.fetchall():
            key = str(row['des_solution_id'])                        
            self.des_solutions[key] = {'status': str(row['status']), 'des_solution_id': str(row['des_solution_id'])}
                        
        self.cur.execute("select generated_solution_id from generated_solution")
        
        for row in self.cur.fetchall():                        
            self.gen_solutions_id[str(row['generated_solution_id'])] = '1'            
        
        return 0
                
        
    def DBGetSolution(self, solution_id):
        '''
        Get solution given solution_id
        returns: a dictionary with a solution with all attributes
        '''
        pass
    
        self.cur.execute("select * from generated_solution where generated_solution_id=?",(solution_id,))
        return dict(self.cur.fetchone())
        
    def DBGetDesiredSolution(self):
        '''
        Get a desired solution that wasn't found yet
        returns: a dictionary with a desired solution or None
        '''
        pass
    
        #Insert buffered solutions into DB
        self.flushSQL()
        
        if self.designMethod.listDesigns == []:
            return None
              
        # get a desired solution from db
        self.cur.execute("select * from desired_solution where status='WAITING' order by random() LIMIT 1")
        des_solution = self.cur.fetchone()
        
        if des_solution != None:
            des_solution = dict(des_solution)
            #des_sol = des_solutions[0]
            # set worker as working on desired solution
            start_time = strftime("%Y-%m-%d %H:%M:%S %Z")
            self.cur.execute("update desired_solution set worker_id=?, status=?, start_time=? where des_solution_id = ?;", (self.worker_id, 'RUNNING', start_time, des_solution['des_solution_id'] ) )
        
            return des_solution
        else: #no more solutions waiting, so help doing the ones that are currently running
            self.cur.execute("select * from desired_solution where status='RUNNING' LIMIT 20")
            des_solutions = (self.cur.fetchall())                    
            
            if des_solutions != []:
                return dict(choice(des_solutions))
            else:
                return None
            
    def DBChangeStatusDesiredSolution(self, desired_solution_id, status='WAITING'):
        self.cur.execute("update desired_solution set worker_id=NULL, status=?, start_time=NULL where des_solution_id = ?;", (status, desired_solution_id ) )    
    
    def DBGetClosestSolution(self,desiredSolution):
        '''
        Get a solution that is closer to the desired solution
        returns: a dictionary with a solution with all attributes
        '''
        
        if self.designMethod.listDesigns == [] or desiredSolution == None:
            # if there is no desiredSolution return a random solution
            query = "SELECT generated_solution_id, sequence FROM generated_solution ORDER BY RANDOM() limit 1;"
            self.cur.execute(query)

            return self.cur.fetchone()
                
        #get closer solutions and distances from a sample of 1000
        query = "SELECT * FROM generated_solution AS r1  JOIN (SELECT (ABS(RANDOM() % (select count(*) from generated_solution))) as selid FROM generated_solution limit 5000) as r2 ON r1.rowid == r2.selid"        
        
        self.cur.execute(query,desiredSolution)                        
        all_solutions = (self.cur.fetchall())
        
        #print all_solutions
        distanceArray = [self.distanceBetweenSolutions(sol_i, desiredSolution) for sol_i in all_solutions]
        total_fit = sum([1/(dist+0.0001) for dist in distanceArray])
        p_array = [((1/dist)/total_fit) for dist in distanceArray]
        
        if all_solutions == []:
            return None
        else:
            rnd_sol_indx = pick_random(p_array)
            return dict(all_solutions[rnd_sol_indx])
        pass
    
    
    def DBCheckDesign(self, desired_solution_id):
        '''
        Get the status of a solution to design
        returns: a boolean with the result of status == 'DONE'        
        '''
        pass    
        self.cur.execute("select * from desired_solution where des_solution_id=?",(desired_solution_id,))
        return self.cur.fetchone()['status'] == "DONE"
    
    def DBInsertSolution(self,solution,desired_solution_id=""):
        '''
        Insert solution into database
        returns: Nothing
        '''
        pass
    
        if solution.solid in self.gen_solutions_id or solution.valid == False:
            return 0             
        else:
            self.gen_solutions_id[solution.solid] = '1'
            
        key = '.'.join([str(solution.levels[feature+'Level']) for feature in self.designMethod.featuresList]) 
        
        if not self.designMethod.listDesigns==[]: #RandomSampling mode does not have desired targets
            if desired_solution_id == "": #Worker found solution for something it WASN'T looking for                
                if key in self.des_solutions:            
                    desired_solution_id = str(self.des_solutions[key]['des_solution_id'])
                    if self.des_solutions[key]['status'] != 'DONE':
                        self.des_solutions[key]['status'] = 'DONE'
                        self.des_solutions_sql.append({'worker_id' : self.worker_id, 'status': 'DONE', 'des_solution_id' : desired_solution_id})                    
            else:
                self.des_solutions[key]['status'] = 'DONE'                                                  
                self.des_solutions_sql.append({'worker_id' : self.worker_id, 'status': 'DONE', 'des_solution_id' : desired_solution_id})
        else:
            desired_solution_id = key
                        
        #update generated solution table
        dict_with_values = {'generated_solution_id' : solution.solid, 
                            'des_solution_id': desired_solution_id,  
                            'sequence' :solution.sequence ,
                            'worker_id': self.worker_id}
        dict_with_values.update(solution.scores)
        dict_with_values.update(solution.levels)                
        dict_with_values.update({ (feature+'Position'): self.calculateRelativeLevel(feature,solution.levels[feature+'Level'],solution.scores[feature]) for feature in self.designMethod.features })
        
        self.gen_solutions_sql.append(dict_with_values)
        
    def DBCloseConnection(self):
        '''
        Closes connection to DB
        returns: Nothing
        '''
        #Insert buffered solutions into DB
        self.flushSQL()
        
        finish_time = strftime("%Y-%m-%d %H:%M:%S %Z")
        self.cur.execute("update worker set time_finish = ? where worker_id = ?;", (finish_time, self.worker_id) )
        self.cur.close()    
                                                       

    #############
    # Auxiliary functions
    
    def flushSQL(self):    
        
        # desired solutions
        sql = "update desired_solution set worker_id=:worker_id, status=:status where des_solution_id=:des_solution_id"
        self.cur.executemany(sql , self.des_solutions_sql)
        self.des_solutions_sql[:] = [] #empty list
        
        # generated solutions
        features_fields = ','.join([feature+", "+feature+"Level, "+feature+"Position" for feature in self.designMethod.features])
        features_values_fields = ','.join([":"+feature+", :"+feature+"Level, :"+feature+"Position" for feature in self.designMethod.features])
                           
        sql = "insert into generated_solution(generated_solution_id, des_solution_id, sequence, "+features_fields+",worker_id) \
                                    values(:generated_solution_id, :des_solution_id, :sequence, "+features_values_fields+", :worker_id);"                            
        self.cur.executemany(sql , self.gen_solutions_sql)
        self.gen_solutions_sql[:] = [] #empty list
        
    def distanceBetweenSolutions(self,sol1,levels_sol2):
                        
        euc_dist = 0

        for feature in self.designMethod.features:            
            if levels_sol2[feature+'Level']=='?' or sol1[feature+'Level']=='?':
                #d = int(max(self.designMethod.thresholds[feature].keys()))
                d=1
            elif int(levels_sol2[feature+'Level'])==int(sol1[feature+'Level']):
                d=0                                   
            else:  
                d=(int(levels_sol2[feature+'Level'])-int(sol1[feature+'Level']))
                rel_level = self.calculateRelativeLevel(feature,sol1[feature+'Level'],sol1[feature])
                if d > 0:
                    rel_level = 1-rel_level
                d = abs(d)+rel_level
            euc_dist += d**2
        
        euc_dist = sqrt(euc_dist)
        
        return euc_dist                    

    def calculateRelativeLevel(self,feature="",level=1,featureScore=0):
        
        if level == '?':
            return 0
        
        thresholds = self.designMethod.thresholds[feature][level]
        
        if isinstance(thresholds,tuple): 
            t_max = thresholds[1]
            t_min = thresholds[0]
            
            #TODO how to see how far a solution is when limits are infinity?
            if t_max==None:
                return 0
            elif t_min==None:
                return 0
            
            return float(abs(featureScore-t_min)/abs(t_max-t_min)) 
        
        return 0