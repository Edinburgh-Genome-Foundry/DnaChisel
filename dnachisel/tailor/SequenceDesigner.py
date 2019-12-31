from tailor.DBSQLite import DBSQLite
from tailor.Solution import Solution
from time import time
from random import choice
from tailor.tools import hammingDistance, pick_random_tuple
from math import exp, sqrt
import sys


class SequenceDesigner(object):
    '''
    Initializes class that design sequences based on a design method
    '''
    def __init__(self, name, seed, design, dbfile, createDB=True):

        self.name = name
        self.designMethod = design
        self.dbfile = dbfile
        self.solutionsHash = {}
        self.max_iterations = 100  #maximum number of tries allowed to find desired solution
        self.max_sol_counter = 100000
        
        # run 
        self._statistics_period = 1
        # temperature
        self.temp_ON = True
        self.temp_init = 100
        self.temp_schedule = 0.79

        self.dbconnection = DBSQLite(dbfile=dbfile,
                                     designMethod=design,
                                     initialize=createDB,
                                     seedSequence=seed)

    def configureSolution(self, solution):
        '''
        Solution configuration
        '''
        raise NotImplementedError(
                'configureSolution must be implemented before using it')

    def validateSolution(self, solution):
        '''
        Solution validation tests
        '''
        solution.valid = True

    def additionalConfigurationPreMutation(self, solution):
        '''
        This method is executed before the mutation iteration happens 
        and can be used to set additional mutational properties
        '''
        pass

    def additionalConfigurationPostMutation(self, solution):
        '''
        This method is executed after the mutation iteration happens 
        and can be used to set additional mutational properties
        '''
        pass

    def calculateRelativeLevel(self, feature="", level=1, featureScore=0):
        thresholds = self.designMethod.thresholds[feature][level]

        if isinstance(thresholds, tuple):
            t_max = thresholds[1]
            t_min = thresholds[0]

            #TODO how to see how far a solution is when limits are infinity?
            if t_max == None:
                return 0
            elif t_min == None:
                return 0

            return float(abs(featureScore - t_min) / abs(t_max - t_min))

        return 0

    def distanceBetweenSolutions(self, sol1, levels_sol2):

        levels_sol1 = sol1.levels

        euc_dist = 0

        for feature in self.designMethod.features:
            if levels_sol2[feature +
                           'Level'] == '?' or levels_sol1[feature +
                                                          'Level'] == '?':
                d = 1
            elif int(levels_sol2[feature + 'Level']) == int(
                    levels_sol1[feature + 'Level']):
                d = 0
            else:
                d = (int(levels_sol2[feature + 'Level']) -
                     int(levels_sol1[feature + 'Level']))
                rel_level = self.calculateRelativeLevel(
                    feature, levels_sol1[feature + 'Level'],
                    sol1.scores[feature])
                if d > 0:
                    rel_level = 1 - rel_level
                d = abs(d) + rel_level
            euc_dist += d**2

        euc_dist = sqrt(euc_dist)

        return euc_dist

    def run(self, selection="directional"):

        start_time = time()
        sol_counter = 1 # seed solution
        last_counter = 1
        last_timepoint = time()
        accepted = 1
        initial_dist = 0

        # seed solution
        master = Solution(
            sol_id=self.dbconnection.seedId,
            sequence=self.dbconnection.seedSequence,
            design=self.designMethod
        )
        self.configureSolution(master)
        self.validateSolution(master)
        solution = master

        if not master.valid:
            raise Exception("Seed inserted is not a valid sequence!")

        self.dbconnection.DBInsertSolution(master)
        self.solutionsHash[master.solid] = master

        all_combinations_found = False

        while not all_combinations_found and sol_counter <= self.max_sol_counter:
            iteration = 0

            if time() - last_timepoint >= self._statistics_period:  #Print statistics every 1 period
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

                if self.designMethod.listDesigns != []:  #All desired combinations were found
                    all_combinations_found = True
                    break
            else:
                initial_dist = self.distanceBetweenSolutions(master, desired_solution)
                print("looking for combination: ",
                      desired_solution['des_solution_id'])
                desired_solution_id = desired_solution['des_solution_id']
            """
            parent = master
            solution = parent
            """

            #Choose stochastically a close solution to the desired one
            if choice([True, True, True, True, True, True, True, False]):
                closestSolution = self.dbconnection.DBGetClosestSolution(
                    desired_solution)
            else:
                closestSolution = self.dbconnection.DBGetClosestSolution(None)
#             print(closestSolution)

            if closestSolution != None:
                #print "SolutionIterator: Found close sequence, starting from here..."
                if closestSolution[
                        'generated_solution_id'] in self.solutionsHash:
                    parent = self.solutionsHash[
                        closestSolution['generated_solution_id']]
                else:
                    parent = Solution(
                        sol_id=closestSolution['generated_solution_id'],
                        sequence=closestSolution['sequence'],
                        design=self.designMethod)
                    self.configureSolution(parent)
                    self.validateSolution(parent)

                solution = parent
            else:
                #print "SolutionIterator: Starting from master sequence"
                parent = master
                solution = parent

            found = False
            old_solution = solution

            # Sequence evolution cycle
            while (
                not solution.checkSolution(desired_solution) 
                and solution.valid 
                and iteration != self.max_iterations 
                and not found 
                and not all_combinations_found
            ):

                if solution != parent:
                    self.dbconnection.DBInsertSolution(solution)
                    self.solutionsHash[
                        solution.solid] = solution  ### Cache for rapid access
                    sol_counter += 1

                #generate next solution
                if selection == "neutral":
                    solution = choice([parent, solution])
                    
                elif selection == "directional":
                    if self.designMethod.listDesigns != []:
                        dist_old = self.distanceBetweenSolutions(
                            old_solution, desired_solution)
                        dist_cur = self.distanceBetweenSolutions(
                            solution, desired_solution)

                        if dist_old < dist_cur:
                            solution = old_solution

                
                elif selection == "temperature":
                    if self.designMethod.listDesigns != []:
                        dist_old = self.distanceBetweenSolutions(
                            old_solution, desired_solution)
                        dist_cur = self.distanceBetweenSolutions(
                            solution, desired_solution)

                        delta = dist_cur - dist_old

                        try:
                            prob = min([
                                exp(-delta /
                                    (self.temp_init *
                                     (self.temp_schedule**iteration))), 1
                            ])
                        except OverflowError:
                            prob = 0 if delta > 0 else 1

                        solution = pick_random_tuple([(old_solution, 1 - prob),
                                                      (solution, prob)])

                        if solution != old_solution and delta > 0:
                            accepted = accepted + 1

                    pass
                else:
                    #use current solution as parent for next round of mutations
                    sys.stderr.write(
                        "Selection option selected is not available, using 'directional' instead...\n"
                    )
                    selection = 'directional'

                self.additionalConfigurationPreMutation(solution) 

                old_solution = solution
                solution = solution.mutate(desired_solution)

                # No solution found
                if solution == None or solution.sequence == None:
                    solution = None
                    break

                self.configureSolution(solution)
                self.validateSolution(solution)

                self.additionalConfigurationPostMutation(solution)

                #go to next iteration
                iteration += 1

                #check if my desired solution was already found
                if self.designMethod.listDesigns != [] and iteration % (
                        self.max_iterations / 2) == 0:
                    found = self.dbconnection.DBCheckDesign(
                        desired_solution_id)

                if self.designMethod.listDesigns == []:
                    #Stops when number generated solutions is equal to the desired sample size
                    if sol_counter >= self.designMethod.nDesigns:
                        all_combinations_found = True
                        print("RandomSampling: %s solutions generated." %
                              (sol_counter))

            #insert solution in the DB
            if (
                solution != None 
                and solution.checkSolution(desired_solution) 
                and solution != parent 
                and solution.valid
            ):
                print("Solution found... inserting into DB...")
                self.dbconnection.DBInsertSolution(solution,desired_solution_id)
                self.solutionsHash[solution.solid] = solution
                sol_counter += 1
            elif found == True:
                print("Solution already found by other worker")
            else:
                if self.designMethod.listDesigns != [] and not all_combinations_found:
                    print("No solution could be found...")
                    self.dbconnection.DBChangeStatusDesiredSolution(
                        desired_solution_id, 'WAITING')

        #set worker as finished
        self.dbconnection.DBCloseConnection()

        if len(self.designMethod.listDesigns) == 1:
            print("\n###########################")
            print("# Optimized solution:")
            print("# ID: ", solution.solid)
            print("# Sequence: ", solution.sequence)
            print("# Scores: ", [
                feat + ": " + str(solution.scores[feat])
                for feat in self.designMethod.featuresList
            ])
            print("# Levels: ", [
                feat + "Level: " + str(solution.levels[feat + "Level"])
                for feat in self.designMethod.featuresList
            ])
            print("# Number of generated solutions: ", sol_counter)
            print("# Distance to seed: ",
                  hammingDistance(master.sequence, solution.sequence))
            print("###########################\n")

        print("Program finished, all combinations were found...")
        return (
            sol_counter, 
            hammingDistance(master.sequence,solution.sequence), 
            initial_dist)
