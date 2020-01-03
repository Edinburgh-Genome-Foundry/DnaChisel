from random import choice
from time import time
from uuid import uuid4

from .DBSQLite import DBSQLite
from .tools import hammingDistance

class SequenceDesigner(object):
    '''
    Initializes class that design sequences based on a design method
    '''
    def __init__(self, name, seed, design_space, dbfile, createDB=True):

        self.name = name
        self.design_space = design_space
        self.dbfile = dbfile
        self.max_iterations = 100  #maximum number of tries allowed to find desired solution
        self.max_sol_counter = 100000

        self.solutionsHash = {}
        
        self.dbconnection = DBSQLite(
            dbfile=dbfile,
            design_space=design_space,
            initialize=createDB,
            seedSequence=seed
        )

    def define_problem(self,sequence,solution_id):
        '''
        define a design problem
        '''
        raise NotImplementedError(
            'define_problem must be implemented before using it')
        
    def run(self):
        start_time = time()
        sol_counter = 1  # seed solution
        last_counter = 1
        last_timepoint = time()
        accepted = 1
        initial_dist = 0
        
        seed_sequence = self.dbconnection.seedSequence
        solution_id=self.dbconnection.seedId
        
        self._solution = self.define_problem(seed_sequence,solution_id)
        mutation_space = self._solution.mutation_space
        
        self.dbconnection.DBInsertSolution(self._solution)
        self.solutionsHash[self._solution.solution_id] = self._solution
        
        all_combinations_found = False
        while not all_combinations_found and sol_counter <= self.max_sol_counter:
            iteration = 0

            if time() - last_timepoint >= 1:  #Print statistics every 1 second
                print(
                    "time elapsed: %.2f (s) \t solutions generated: %d \t rate (last period.): %0.2f sol/s  \t rate (overall): %0.2f sol/s"
                    % ((time() - start_time), sol_counter,
                       (sol_counter - last_counter) /
                       (time() - last_timepoint), sol_counter /
                       (time() - start_time)))
                last_counter = sol_counter
                last_timepoint = time()
            
            # Retrieve some desired solution (i.e. a particular combination of features that was not yet found)
            desired_solution = self.dbconnection.DBGetDesiredSolution()

            if desired_solution == None:  #There are no more desired solutions

                if self.design_space.designs_list != []:  #All desired combinations were found
                    all_combinations_found = True
                    break
            else:
                initial_dist = self.design_space._distance_between_designs(
                    self._solution.scores, 
                    desired_solution['des_solution_id'])
                print("looking for combination: ",
                      desired_solution['des_solution_id'])
                desired_solution_id = desired_solution['des_solution_id']
                
            if choice([True, True, True, True, True, True, True, False]):
                closestSolution = self.dbconnection.DBGetClosestSolution(
                    desired_solution)
            else:
                closestSolution = self.dbconnection.DBGetClosestSolution(None)
                
            if closestSolution != None:
                print ("SolutionIterator: Found close sequence, starting from here...")
                if closestSolution[
                        'generated_solution_id'] in self.solutionsHash:
                    parent = self.solutionsHash[
                        closestSolution['generated_solution_id']]
                else:
                    parent = self.define_problem( 
                        closestSolution['sequence'],
                        closestSolution['generated_solution_id'],
                        )
                solution = parent
            else:
                print ("SolutionIterator: Starting from init sequence")
                parent = self._solution
                solution = parent
            
            found = False
            old_solution = solution

            # Sequence evolution cycle
            while (not solution.is_match_design(desired_solution)
                   and iteration != self.max_iterations
                   and not found 
                   and not all_combinations_found):

                if solution != parent:
                    self.dbconnection.DBInsertSolution(solution)
                    self.solutionsHash[
                        solution.solution_id] = solution  ### Cache for rapid access
                    sol_counter += 1


                if self.design_space.designs_list != []:
                    dist_old = self.design_space._distance_between_designs(
                        old_solution.scores, 
                        desired_solution['des_solution_id'])

                    dist_cur = self.design_space._distance_between_designs(
                        solution.scores, 
                        desired_solution['des_solution_id'])

                    if dist_old < dist_cur:
                        solution = old_solution


                old_solution = solution
                
                
                new_sequence = old_solution.mutation_space.apply_random_mutations(
                    n_mutations = choice([1,2]),
                    sequence=old_solution.sequence)
                
                solution =  self.define_problem(new_sequence,str(uuid4().int))

                # No solution found
                if solution == None or solution.sequence == None:
                    solution = None
                    break
                #go to next iteration
                iteration += 1
                

                #check if my desired solution was already found
                if self.design_space.designs_list != [] and iteration % (
                        self.max_iterations / 2) == 0:
                    found = self.dbconnection.DBCheckDesign(
                        desired_solution_id)


                if self.design_space.designs_list == []:
                    #Stops when number generated solutions is equal to the desired sample size
                    if sol_counter >= self.design_space.ndesigns:
                        all_combinations_found = True
                        print("RandomSampling: %s solutions generated." %
                              (sol_counter))                
                
            #insert solution in the DB
            if (solution != None and solution.is_match_design(desired_solution_id)
                    and solution != parent ):
                print("Solution found... inserting into DB...")
                self.dbconnection.DBInsertSolution(solution,
                                                   desired_solution_id)
                self.solutionsHash[solution.solution_id] = solution
                sol_counter += 1
            elif found == True:
                print("Solution already found by other worker")
            else:
                if self.design_space.designs_list != [] and not all_combinations_found:
                    print("No solution could be found...")
                    self.dbconnection.DBChangeStatusDesiredSolution(
                        desired_solution_id, 'WAITING')
                

        #set worker as finished
        self.dbconnection.DBCloseConnection()

        if len(self.design_space.designs_list) == 1:
            print("\n###########################")
            print("# Optimized solution:")
            print("# ID: ", solution.solution_id)
            print("# Sequence: ", solution.sequence)
            print("# Scores: ", [
                feat + ": " + str(solution.scores[feat])
                for feat in self.design_space.feature_label
            ])
            print("# Levels: ", [
                feat + "_Level: " + str(solution.levels[feat + "_Level"])
                for feat in self.design_space.feature_label
            ])
            print("# Number of generated solutions: ", sol_counter)
            print("# Distance to seed: ",
                  hammingDistance(self._solution.sequence, solution.sequence))
            print("###########################\n")

        print("Program finished, all combinations were found...")
        
        return (sol_counter, hammingDistance(self._solution.sequence,
                                             solution.sequence), initial_dist)
