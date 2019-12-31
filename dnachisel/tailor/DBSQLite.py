from math import sqrt
from random import choice
from socket import gethostbyname, gethostname
from sqlite3 import Row, connect
from subprocess import check_output
from time import strftime
from uuid import uuid4

from .tools import pick_random


class DBSQLite():
    '''
    Constructor
    '''
    def __init__(self,
                 dbfile,
                 designMethod,
                 initialize=True,
                 seedSequence=None):
        self.dbfile = dbfile
        self.designMethod = designMethod
        self.seedSequence = seedSequence
        self.seedId = str(uuid4().int)

        #SQL queries buffers
        self._des_solutions = {}
        self._des_solutions_sql = []
        self._gen_solutions_sql = []
        self._gen_solutions_id = {}

        #Connect database
        self._con = connect(self.dbfile + ".sqlite")
        self._con.row_factory = Row
        self._con.isolation_level = None
        self._cur = self._con.cursor()

        #Initiate DB
        if initialize:
            self._create_table()
            self._insert_ds_table()

        #Register worker
        self.worker_id = str(uuid4().int)
        self._insert_worker_table()
        self._registerWorker()



    def DBInsertSolution(self, solution, desired_solution_id=""):
        '''
        Insert solution into database
        returns: Nothing
        '''
        if solution.solid in self._gen_solutions_id or solution.valid == False:
            return 0
        else:
            self._gen_solutions_id[solution.solid] = '1'

        key = '.'.join([
            str(solution.levels[feature + 'Level'])
            for feature in self.designMethod.featuresList
        ])

        if not self.designMethod.listDesigns == []:  #RandomSampling mode does not have desired targets
            if desired_solution_id == "":  #Worker found solution for something it WASN'T looking for
                if key in self._des_solutions:
                    desired_solution_id = str(
                        self._des_solutions[key]['des_solution_id'])
                    if self._des_solutions[key]['status'] != 'DONE':
                        self._des_solutions[key]['status'] = 'DONE'
                        self._des_solutions_sql.append({
                            'worker_id':self.worker_id,
                            'status':'DONE',
                            'des_solution_id':desired_solution_id
                        })
            else:
                self._des_solutions[key]['status'] = 'DONE'
                self._des_solutions_sql.append({
                    'worker_id':self.worker_id,
                    'status':'DONE',
                    'des_solution_id':desired_solution_id
                })
        else:
            desired_solution_id = key

        #update generated solution table
        dict_with_values = {
            'generated_solution_id': solution.solid,
            'des_solution_id': desired_solution_id,
            'sequence': solution.sequence,
            'worker_id': self.worker_id
        }
        dict_with_values.update(solution.scores)
        dict_with_values.update(solution.levels)
        dict_with_values.update({
            (feature + 'Position'):
            self._calculateRelativeLevel(feature,
                                        solution.levels[feature + 'Level'],
                                        solution.scores[feature])
            for feature in self.designMethod.features
        })

        self._gen_solutions_sql.append(dict_with_values)


    def DBGetSolution(self, solution_id):
        '''
        Get solution given solution_id
        returns: a dictionary with a solution with all attributes
        '''

        self._cur.execute(
            "select * from generated_solution where generated_solution_id=?",
            (solution_id, ))
        return dict(self._cur.fetchone())

    def DBGetDesiredSolution(self):
        '''
        Get a desired solution that wasn't found yet
        returns: a dictionary with a desired solution or None
        '''

        #Insert buffered solutions into DB
        self._flushSQL()

        if self.designMethod.listDesigns == []:
            return None

        # get a desired solution from db
        self._cur.execute(
            """select * from desired_solution 
                where status='WAITING' 
                order by random() 
                LIMIT 1"""
        )
        des_solution = self._cur.fetchone()
        if des_solution != None:
            des_solution = dict(des_solution)
            #des_sol = des_solutions[0]
            # set worker as working on desired solution
            start_time = strftime("%Y-%m-%d %H:%M:%S %Z")
            self._cur.execute(
                "update desired_solution set worker_id=?, status=?, start_time=? where des_solution_id = ?;",
                (self.worker_id, 'RUNNING', start_time,
                 des_solution['des_solution_id']))

            return des_solution
        else:  #no more solutions waiting, so help doing the ones that are currently running
            self._cur.execute(
                "select * from desired_solution where status='RUNNING' LIMIT 20"
            )
            des_solutions = (self._cur.fetchall())

            if des_solutions != []:
                return dict(choice(des_solutions))
            else:
                return None

    def DBGetClosestSolution(self, desiredSolution):
        '''
        Get a solution that is closer to the desired solution
        returns: a dictionary with a solution with all attributes
        '''

        if self.designMethod.listDesigns == [] or desiredSolution == None:
            # if there is no desiredSolution return a random solution
            query = """SELECT generated_solution_id, sequence 
                        FROM generated_solution 
                        ORDER BY RANDOM() 
                        limit 1;"""
            self._cur.execute(query)

            return self._cur.fetchone()

        #get closer solutions and distances from a sample of 1000
        query = """SELECT * FROM generated_solution AS r1  
        JOIN 
        (SELECT (ABS(RANDOM() % (select count(*) from generated_solution))) as selid 
            FROM generated_solution limit 5000) as r2 
        ON r1.rowid == r2.selid"""

        self._cur.execute(query, desiredSolution)
        all_solutions = (self._cur.fetchall())

        #print all_solutions
        distanceArray = [
            self.distanceBetweenSolutions(sol_i, desiredSolution)
            for sol_i in all_solutions
        ]
        total_fit = sum([1 / (dist + 0.0001) for dist in distanceArray])
        p_array = [((1 / dist) / total_fit) for dist in distanceArray]

        if all_solutions == []:
            return None
        else:
            rnd_sol_indx = pick_random(p_array)
            return dict(all_solutions[rnd_sol_indx])


    def DBChangeStatusDesiredSolution(self,
                                      desired_solution_id,
                                      status='WAITING'):
        self._cur.execute(
            "update desired_solution set worker_id=NULL, status=?, start_time=NULL where des_solution_id = ?;",
            (status, desired_solution_id))

    def DBCheckDesign(self, desired_solution_id):
        '''
        Get the status of a solution to design
        returns: a boolean with the result of status == 'DONE'        
        '''
        self._cur.execute(
            "select * from desired_solution where des_solution_id=?",
            (desired_solution_id, ))
        return self._cur.fetchone()['status'] == "DONE"


    def DBCloseConnection(self):
        '''
        Closes connection to DB
        returns: Nothing
        '''
        #Insert buffered solutions into DB
        self._flushSQL()

        finish_time = strftime("%Y-%m-%d %H:%M:%S %Z")
        self._cur.execute(
            "update worker set time_finish = ? where worker_id = ?;",
            (finish_time, self.worker_id))
        self._cur.close()

    #############
    # init
    def _create_table(self):
        '''
        Initialize database
        input: dbfile
            file to save DB        
        input: designMethod
            a class Design with information about features name/type and design wanted
        returns: Nothing
        '''

        #Design Dynamic tables
        features_fields = ''.join([
            feature + "Level TEXT, "
            for feature in self.designMethod.featuresList
        ])

        table_ds = """create table desired_solution(
            des_solution_id TEXT PRIMARY KEY,
            {}
            status TEXT,
            worker_id TEXT,
            start_time TEXT,
            FOREIGN KEY(worker_id) REFERENCES worker(worker_id));""".format(features_fields)

        features_values_fields = ''.join([
            feature + " " + self.designMethod.features[feature]['type'] + ", "
            for feature in self.designMethod.featuresList
        ])
        features_level_position_fields = ''.join([
            feature + "Position " +
            self.designMethod.features[feature]['type'] + ", "
            for feature in self.designMethod.featuresList
        ])

        table_gs = "create table generated_solution("
        table_gs += "generated_solution_id TEXT PRIMARY KEY, des_solution_id TEXT, sequence TEXT,"
        table_gs += features_values_fields+features_fields+features_level_position_fields
        table_gs += "worker_id TEXT, FOREIGN KEY(worker_id) REFERENCES worker(worker_id));"

        table_worker = """create table worker(
                            worker_id TEXT PRIMARY KEY, 
                            hostname TEXT, 
                            ip TEXT, 
                            time_start INTEGER, 
                            time_finish INTEGER);
                            """

        #Create Tables
        self._cur.executescript("""
                        PRAGMA writable_schema = 1;
                        delete from sqlite_master where type = 'table';
                        PRAGMA writable_schema = 0;

                        VACUUM;
                        PRAGMA INTEGRITY_CHECK;
                        
                        """ + 
                        table_worker + 
                        table_ds + "\n" + 
                        table_gs + "\n"
                        )

        #Populate Tables
    
    def _insert_ds_table(self):
        #Desired Solutions DB
        n_features = self.designMethod.n_features
        all_comb = [
            tuple(d_sol.split('.')) + (d_sol, )
            for d_sol in self.designMethod.listDesigns
        ]
        features_levels_fields = ''.join([
            feature + "Level, " for feature in self.designMethod.featuresList
        ])

        sql = "insert into desired_solution(" + features_levels_fields 
        sql += "des_solution_id,status, worker_id, start_time) values(" 
        sql += "?," * (n_features + 1) + "'WAITING',NULL,NULL);"

        self._cur.executemany(sql, all_comb)

    def _insert_worker_table(self):
        start_time = strftime("%Y-%m-%d %H:%M:%S %Z")
        hostname = check_output("hostname").rstrip()
        ip = gethostbyname(gethostname()).rstrip()
        self._cur.execute(
            "insert into worker(worker_id, hostname, ip, time_start, time_finish)"
             " values (?,?,?,?,NULL);",
            (self.worker_id, hostname, ip, start_time))

    def _registerWorker(self):
        self._cur.execute("select * from desired_solution")
        for row in self._cur.fetchall():
            key = str(row['des_solution_id'])
            self._des_solutions[key] = {
                'status': str(row['status']),
                'des_solution_id': str(row['des_solution_id'])
            }

        self._cur.execute("select generated_solution_id from generated_solution")
        for row in self._cur.fetchall():
            self._gen_solutions_id[str(row['generated_solution_id'])] = '1'

        return 0

    #############
    # Auxiliary functions

    def _flushSQL(self):

        # desired solutions
        sql = "update desired_solution set worker_id=:worker_id, status=:status "
        sql += "where des_solution_id=:des_solution_id"
        self._cur.executemany(sql, self._des_solutions_sql)
        self._des_solutions_sql = []  #empty list

        # generated solutions
        features_fields = ','.join([
            '{0}, {0}Level, {0}Position'.format(feature)
            # feature + ", " + feature + "Level, " + feature + "Position"
            for feature in self.designMethod.features
        ])
        features_values_fields = ','.join([
            ':{0}, :{0}Level, :{0}Position'.format(feature)
            # ":" + feature + ", :" + feature + "Level, :" + feature + "Position"
            for feature in self.designMethod.features
        ])

        sql = """insert into generated_solution(
                generated_solution_id, 
                des_solution_id, 
                sequence, 
                {features_fields},
                worker_id
            ) 
            values(
                :generated_solution_id, 
                :des_solution_id, 
                :sequence, 
                {features_values_fields}, 
                :worker_id
            );""".format(
                features_fields = features_fields,
                features_values_fields =features_values_fields)
        self._cur.executemany(sql, self._gen_solutions_sql)
        self._gen_solutions_sql = []  #empty list

    def distanceBetweenSolutions(self, sol1, levels_sol2):

        euc_dist = 0

        for feature in self.designMethod.features:
            if levels_sol2[feature + 'Level'] == '?' or sol1[feature +
                                                             'Level'] == '?':
                #d = int(max(self.designMethod.thresholds[feature].keys()))
                d = 1
            elif int(levels_sol2[feature + 'Level']) == int(sol1[feature +
                                                                 'Level']):
                d = 0
            else:
                d = (int(levels_sol2[feature + 'Level']) -
                     int(sol1[feature + 'Level']))
                rel_level = self._calculateRelativeLevel(
                    feature, sol1[feature + 'Level'], sol1[feature])
                if d > 0:
                    rel_level = 1 - rel_level
                d = abs(d) + rel_level
            euc_dist += d**2

        euc_dist = sqrt(euc_dist)

        return euc_dist

    def _calculateRelativeLevel(self, feature="", level=1, featureScore=0):

        if level == '?':
            return 0

        thresholds = self.designMethod.thresholds[feature][level]

        if isinstance(thresholds, tuple):
            t_max = thresholds[1]
            t_min = thresholds[0]

            #TODO how to see how far a solution is when limits are infinity?
            if t_max == None:
                return 0
            elif t_min == None:
                return 0

            return abs(featureScore - t_min) / abs(t_max - t_min)

        return 0
