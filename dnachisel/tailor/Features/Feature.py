from uuid import uuid4
from ..mutation import randomMutationOperator
from ..Solution import Solution

class Feature(object):
    '''
    A master class for Solution Features to be manipulated
    '''
    def __init__(self, featureObject=None, solution=None, label=""):
        
        if featureObject == None: #create a new instance of object
            self.scores             = {}        
            self.solution           = solution
            self.label              = label
            self.targetInstructions = {}
            self.subfeatures        = {}
            self.level              = None          
        else: #copy instance
            self.solution           = featureObject.solution
            self.label              = featureObject.label
            self.targetInstructions = featureObject.targetInstructions
            self.subfeatures        = {}
            self.scores             = featureObject.scores
            self.level              = featureObject.level

    def set_scores(self):
        '''
        Call to the function calculating the score with proper parameters
        To be implemented in subclasses
        '''
        pass
    
    def add_subfeature(self,feature):
        self.subfeatures[feature.label+feature.__class__.__name__] = feature
    
    
    def getTargets(self,desiredSolution):
        '''
        given a desired solution, returns all the features that need to be modified
        '''
        targets = []
        #evaluate if goal has achieved for base class
        if self.defineTarget(desiredSolution):
            targets.append(self)
                    
        #evaluate if goal has achieved for children classes
        for subfeature in self.subfeatures.values():
            targets.extend(subfeature.getTargets(desiredSolution))
                
        return targets

    def defineTarget(self,desiredSolution):
        '''
        Function that determines if a target wasn't hit and, if not, updates target instructions 
        '''
        if desiredSolution == None:
            return True
        
        #check if there is a target
        if (self.label+self.__class__.__name__+"Level") not in  desiredSolution:
            return False
        else:
            target_level = desiredSolution[self.label+self.__class__.__name__+"Level"]
            
            if target_level == 0:
                return False        
                        
            if  target_level != self.level:
                
                level_info = self.solution.designMethod.thresholds[self.label+self.__class__.__name__][target_level]
                
                if isinstance(level_info, tuple): #Then it's a numeric range                
                    if  level_info[0]-self.scores[self.label+self.__class__.__name__] > 0:
                        self.targetInstructions['direction'] = '+' #increase
                    elif level_info[0]-self.scores[self.label+self.__class__.__name__] < 0:
                        self.targetInstructions['direction'] = '-' #decrease
                elif isinstance(level_info, int): #numeric variable
                    if  level_info-self.scores[self.label+self.__class__.__name__] > 0:
                        self.targetInstructions['direction'] = '+' #increase
                    elif level_info-self.scores[self.label+self.__class__.__name__] < 0:
                        self.targetInstructions['direction'] = '-' #decrease
                elif isinstance(level_info, str): #nominal variable
                    self.targetInstructions['direction'] = level_info 
                else:
                    self.targetInstructions['direction'] = 'NA' #not applicable                                    
                    
                return True

            return False

        
    def set_level(self):
        '''
        define levels and update solution levels dictionary (only works for Numeric scores)
        '''
        
        if self.solution.designMethod != None: #Design mode
            if (self.label+self.__class__.__name__) in self.solution.designMethod.thresholds:
                
                for level_name in self.solution.designMethod.thresholds[self.label+self.__class__.__name__].keys():
                    
                    level_info = self.solution.designMethod.thresholds[self.label+self.__class__.__name__][level_name]                    
                    
                    if self.scores[self.label+self.__class__.__name__]==None or self.scores[self.label+self.__class__.__name__]=="NA":
                        self.level = "NA"
                    elif isinstance(level_info, tuple): #Then it's a numeric range
                        if level_info[0] == None and self.scores[self.label+self.__class__.__name__] <= level_info[1]: #to minus infinity
                            self.level = level_name
                        elif level_info[1] == None and self.scores[self.label+self.__class__.__name__] >= level_info[0]: #to plus infinity
                            self.level = level_name
                        elif level_info[0]!=None and level_info[1]!=None and round(self.scores[self.label+self.__class__.__name__],4) >= level_info[0] and round(self.scores[self.label+self.__class__.__name__],4) <= level_info[1]:                        
                            self.level = level_name
                    elif isinstance(level_info, (list,set)): #Then your level is a set of possible states
                        if self.scores[self.label+self.__class__.__name__] in level_info:
                            self.level = level_name                
                    else: #Then your level is a literal (either string or number)
                        if self.scores[self.label+self.__class__.__name__] == level_info:
                            self.level = level_name
                            
                if self.level == None :
                    #if isinstance(level_info, tuple):
                    #    sys.stderr.write("Feature: Level not defined... " + self.label+self.__class__.__name__ + " -> " + str(self.scores[self.label+self.__class__.__name__]) + "\n")
                    self.level = '?'
        else: #Analysis mode
            self.level = 'NA'
        
        return 

    def randomMutation(self,pos=None,n_mut=[1,2],mutable_region=None):
        if mutable_region == None:
            if self.mutable_region == None:
                mutable_region = self.solution.mutable_region
            else:
                mutable_region = self.mutable_region
                                      
        new_seq = randomMutationOperator(self.solution.sequence, self.solution.keep_aa, mutable_region, self.solution.cds_region, pos,n_mut=n_mut)
        
        return Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.solution.cds_region, keep_aa = self.solution.keep_aa, mutable_region = self.solution.mutable_region, parent=self.solution, design=self.solution.designMethod)

        
    def mutate(self,mutable_region=None):
        '''
        Specify how to call operator to mutate the sequence
        '''
        pass
        return self.randomMutation(mutable_region=mutable_region)
    
