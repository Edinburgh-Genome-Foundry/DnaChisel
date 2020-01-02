import sys
from random import choice
from uuid import uuid4

from .mutation import randomMutationOperator


class Solution:
    '''
    A Solution encapsulates a sequence and their inherent attributes:
        solution_id - ID for Solution
        seqence - sequence for Solution
        cds_region - a tuple indicating the location of (Begin,End) of CDS sequence (this will be necessary in the design mode if one want to contrain mutations).
        mutable_region - a list with all positions that can be mutated
        parent - Solution from which the current Solution was derived
        
    '''
    def __init__(self,
                 solution_id=0,
                 sequence="",
                 cds_region=None,
                 keep_aa=False,
                 mutable_region=None,
                 parent=None,
                 design=None):

        if sequence == None:
            sys.stderr.write("Tried to create a solution with sequence NULL\n")
            self.sequence = None
            return None

        #check if solution is in DB
        self.mutable_region = mutable_region
        self.cds_region = cds_region
        self.keep_aa = keep_aa
        self.solid = solution_id
        self.parent = parent
        self.sequence = sequence.lower()
        self.scores = {}
        self.levels = {}
        self.features = {}
        self.designMethod = design
        self.valid = True

    def add_feature(self, feature):
        featureLabel = feature.label + feature.__class__.__name__
        if featureLabel not in self.features:
            self.features[featureLabel] = feature
            #update scores
            self.scores.update(feature.scores)
            #update levels
            if feature.level != None:
                self.levels[featureLabel + "Level"] = feature.level
            for subfeature in feature.subfeatures.values():
                self.add_feature(subfeature)
        else:
            sys.stderr.write("Feature label already exists!")

        return

    def checkSolution(self, desiredSolution):

        if desiredSolution == None:
            return False

        same = True
        for feature in self.designMethod.features.keys():
            key = feature + "Level"
            same = same & (desiredSolution[key] == 0
                           or desiredSolution[key] == self.levels[key])

        return same

    def mutate(self, desiredSolution=None, random=False):

        if desiredSolution==None or random or self.designMethod.listDesigns == [] or self.features == {}:
            return self.randomMutation()
        else:
            # get features with targets
            mutable = []
            for feature in self.features.values():
                if feature.defineTarget(desiredSolution):
                    mutable.append(feature)

            if mutable == []:
                return None

            rm = choice(mutable)
            #tomutatefeatures = [k.label+k.__class__.__name__ for k in mutable]
            #print tomutatefeatures
            #print [self.scores[k] for k in tomutatefeatures]
            #print [self.levels[k+"Level"] for k in tomutatefeatures]
            #print [desiredSolution[k+"Level"] for k in tomutatefeatures]
            #print "mutating... " + rm.label+rm.__class__.__name__

            return rm.mutate()
            #return choice(mutable).randomMutation()
            #return self.randomMutation()

    def randomMutation(self, pos=None, n_mut=[1, 2]):
        new_seq = randomMutationOperator(self.sequence, self.keep_aa,
                                         self.mutable_region, self.cds_region,
                                         pos, n_mut)
        return Solution(solution_id=str(uuid4().int),
                        sequence=new_seq,
                        cds_region=self.cds_region,
                        keep_aa=self.keep_aa,
                        mutable_region=self.mutable_region,
                        parent=self,
                        design=self.designMethod)
