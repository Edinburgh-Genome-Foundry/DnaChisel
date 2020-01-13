from math import sqrt
from random import choice
from socket import gethostbyname, gethostname
from sqlite3 import Row, connect
from subprocess import check_output
from time import strftime
from uuid import uuid4

from .mathtools import pick_random


class DBSQLite():
    '''
    Constructor
    '''
    def __init__(self,
                 dbfile,
                 design_space,
                 initialize=True,
                 seedSequence=None):
        self.dbfile = dbfile
        self.design_space = design_space
        self.seedSequence = seedSequence
        self.seedId = str(uuid4().int)

        #SQL queries buffers
        self._desired_designs = {}
        self._desired_designs_sql = []
        self._generated_solutions_sql = []
        self._generated_solutions_id = {}

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

    #############
    # init
    #############
    def _create_table(self):
        '''
        Initialize database
        input: dbfile
            file to save DB        
        input: design_space
            a class Design with information about features name/type and design wanted
        returns: Nothing
        '''

        #Design Dynamic tables
        features_fields = ''.join([
            feature + "_Level TEXT, "
            for feature in self.design_space.feature_label
        ])

        table_ds = """create table desired_solution(
            des_solution_id TEXT PRIMARY KEY,
            {}
            status TEXT,
            worker_id TEXT,
            start_time TEXT,
            FOREIGN KEY(worker_id) REFERENCES worker(worker_id));""".format(features_fields)

        features_values_fields = ''.join([
            feature + " " + value_type + ", "
            for feature,value_type in zip(
                self.design_space.feature_label,
                self.design_space.value_type)
        ])
        features_level_position_fields = ''.join([
            feature + "_Position REAL" + ", "
            for feature in self.design_space.feature_label
        ])

        table_gs = """create table generated_solution(
            generated_solution_id TEXT PRIMARY KEY, 
            des_solution_id TEXT, 
            sequence TEXT,
            {}
            worker_id TEXT, FOREIGN KEY(worker_id) REFERENCES worker(worker_id));
        """.format(features_values_fields+features_fields+features_level_position_fields)

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
        n_features = len(self.design_space.feature_label)
        all_comb = [
            tuple(d_sol.split('.')) + (d_sol, )
            for d_sol in self.design_space.designs_list
        ]
        features_levels_fields = ''.join([
            feature + "_Level, " for feature in self.design_space.feature_label
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
            self._desired_designs[key] = {
                'status': str(row['status']),
                'des_solution_id': str(row['des_solution_id'])
            }

        self._cur.execute("select generated_solution_id from generated_solution")
        for row in self._cur.fetchall():
            self._generated_solutions_id[str(row['generated_solution_id'])] = '1'

        return 0


    ###############
    # main operation
    ###############
    def DBInsertSolution(self, solution, desired_solution_id=""):
        '''
        Insert solution into database
        returns: Nothing
        '''
        if solution.solution_id in self._generated_solutions_id :
            return 0
        else:
            self._generated_solutions_id[solution.solution_id] = '1'

        key = '.'.join([
            str(solution.levels[feature + '_Level'])
            for feature in self.design_space.feature_label
        ])

        if not self.design_space.designs_list == []:  #RandomSampling mode does not have desired targets
            if desired_solution_id == "":  #Worker found solution for something it WASN'T looking for
                if key in self._desired_designs:
                    desired_solution_id = str(
                        self._desired_designs[key]['des_solution_id'])
                    if self._desired_designs[key]['status'] != 'DONE':
                        self._desired_designs[key]['status'] = 'DONE'
                        self._desired_designs_sql.append({
                            'worker_id':self.worker_id,
                            'status':'DONE',
                            'des_solution_id':desired_solution_id
                        })
            else:
                self._desired_designs[key]['status'] = 'DONE'
                self._desired_designs_sql.append({
                    'worker_id':self.worker_id,
                    'status':'DONE',
                    'des_solution_id':desired_solution_id
                })
        else:
            desired_solution_id = key

        #update generated solution table
        dict_with_values = {
            'generated_solution_id': solution.solution_id,
            'des_solution_id': desired_solution_id,
            'sequence': solution.sequence,
            'worker_id': self.worker_id
        }
        dict_with_values.update(solution.scores)
        dict_with_values.update(solution.levels)
        dict_with_values.update({
            (self.design_space.feature_label[i] + '_Position'):
            self.design_space.range_set_list[i].calculate_relative_level(
                solution.scores[self.design_space.feature_label[i]]
            )
            for i in range(len(self.design_space.feature_label))
        })

        self._generated_solutions_sql.append(dict_with_values)


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

        if self.design_space.designs_list == []:
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
        ! wrong, need fix
        '''

        if self.design_space.designs_list == [] or desiredSolution == None:
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
        # if len(all_solutions) > 0:print(dict(all_solutions[0]))
        distanceArray = [ 
            self.design_space._distance_between_designs(
                    dict(sol_i),
                    desiredSolution['des_solution_id'])
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
    # Auxiliary functions

    def _flushSQL(self):

        # desired solutions
        sql = "update desired_solution set worker_id=:worker_id, status=:status "
        sql += "where des_solution_id=:des_solution_id"
        self._cur.executemany(sql, self._desired_designs_sql)
        self._desired_designs_sql = []  #empty list

        # generated solutions
        features_fields = ','.join([
            '{0}, {0}_Level, {0}_Position'.format(feature)
            for feature in self.design_space.feature_label
        ])
        features_values_fields = ','.join([
            ':{0}, :{0}_Level, :{0}_Position'.format(feature)
            for feature in self.design_space.feature_label
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
        self._cur.executemany(sql, self._generated_solutions_sql)
        self._generated_solutions_sql = []  #empty list

