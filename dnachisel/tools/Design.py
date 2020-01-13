from itertools import product
from math import sqrt

from .range import RangeSet


class DesignSpace:
    '''
    a =   DesignSpace(
        ['cai','gc'],
        [RangeSet([(1, Range(0,3)),(3, Range(3,6)),(2, Range(6,9))]),   
         RangeSet([('a', Range(0,0.3)),('b', Range(0.3,0.6)),('c', Range(0.6,1))]) ],
        ['REAL','REAL']
    )
    '''
    __slot__ = ['feature_label','range_set_list','value_type','designs_list','ndesigns']
    def __init__(self,feature_label,feature_specification,range_set_list,value_type,**args):
        self.feature_label = feature_label
        self.feature_specification = feature_specification
        self.range_set_list = range_set_list
        self.value_type = value_type
        self.designs_list = []
        self.ndesigns  = 0
    
    def __repr__(self):
        text = '''DesignSpace(
            feature_label=%s,
            RangeSet = %s,
            value_type=%s
            designs_list = %s)'''%(
            self.feature_label,self.range_set_list,self.value_type,self.designs_list )
        return text
        
    
    def __getitem__(self,label):
        index = self.feature_label.index(label)
        return self.range_set_list[index]
    
    def _get_value_type(self,label):
        index = self.feature_label.index(label)
        return self.value_type[index]
        
    def _distance_between_designs(self,scores,desired_design_label):
        
        desired_label = desired_design_label.split('.')
        desired_label_dict = {
            self.feature_label[i]: desired_label[i]
            for i in range(len(desired_label))
        }
        
        dists = [
            self[feature]._calculate_distance(scores[feature],desired_label_dict[feature])
            for feature in self.feature_label
        ]
        
        return sqrt(sum([ i*i for i in dists])) 

class CustomDesign(DesignSpace):
    '''
    Class encoding a custom design (as many targets as you want)
    '''
    __slot__ = ['feature_label','range_set_list','value_type','designs_list','ndesigns']
    def __init__(self,feature_label,feature_specification,range_set_list,value_type,**args):

        DesignSpace.__init__(self, feature_label,feature_specification,range_set_list,value_type)
        self.designs_list = args['targets']
        self.nDesigns = self.designs_list.__len__()


class Optimization(DesignSpace):
    '''
    Class encoding a single-target design (optimization) 
    '''
    __slot__ = ['feature_label','range_set_list','value_type','designs_list','ndesigns']
    def __init__(self,feature_label,feature_specification,range_set_list,value_type,**args):

        DesignSpace.__init__(self, feature_label,feature_specification,range_set_list,value_type)
        self.designs_list = [args['target']]
        self.nDesigns = self.designs_list.__len__()


class RandomSampling(DesignSpace):
    '''
    Class encoding a random sampling design
    '''
    __slot__ = ['feature_label','range_set_list','value_type','designs_list','ndesigns']
    def __init__(self,feature_label,feature_specification,range_set_list,value_type,**args):

        DesignSpace.__init__(self, feature_label,feature_specification,range_set_list,value_type)
        self.designs_list = []
        self.nDesigns = args['sample_size']


class FullFactorial(DesignSpace):
    '''
    Class encoding a multi factorial design
    '''
    __slot__ = ['feature_label','range_set_list','value_type','designs_list','ndesigns']
    def __init__(self,feature_label,feature_specification,range_set_list,value_type,**args):

        DesignSpace.__init__(self, feature_label,feature_specification,range_set_list,value_type)
        self.designs_list = self._computeCombinations()
        self.nDesigns = self.designs_list.__len__()

    def _computeCombinations(self):
        range_labels =[list(x) for x in self.range_set_list]
        return ['.'.join(x) for x in  list(product(*range_labels))]