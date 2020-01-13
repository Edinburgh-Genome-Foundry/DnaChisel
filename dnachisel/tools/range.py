class Range:
    """
    range operation for design

    Range(0,6,) means   0 <= x < 6
    """
    def __init__(self,
        start = None,
        end = None,
        bound = 'left' 
        ):
        self.start = start
        self.end = end
        self.bound = bound
    
    def __repr__(self):
        return 'Range(start = %s, end = %s, bound = %s)'%(
        self.start, self.end, self.bound)
        
    def __contains__(self,number):
        return self.start <= number < self.end
    
    def calculate_relative_level(self,number):
        if number not in self:
            raise ValueError('%s not in %s'%(number,self))
        return (number - self.start) / (self.end - self.start )
    
    def _is_overlap_with(self,other):
        first ,second = self, other
        if first.start > second.start:
            first, second = second, first
        
        if first.end <= second.first:
            return False
        return True





class RangeSet:
    """
    RangeSet([
        ('1', Range(0,3)),
        ('2', Range(3,6)),
        ('3', Range(6,9)),
        ])
    """
    def __init__(self,range_list):
        
        self.level_label = [str(i[0]) for i in range_list]
        self.level_range = [i[1] for i in range_list]
        self.level_num = len(self.level_label)
    
    def __repr__(self):
        text = 'RangeSet(['
        for i in range(self.level_num):
            text += "\n\t('%s',%s)"%(self.level_label[i],self.level_range[i])
        text += '\n\t])'
        return text
    
    def __iter__(self):
        for label in self.level_label:
            yield label


    def _index(self,number):
        for i in range(self.level_num):
            if number in self.level_range[i]:
                return i

    def get_level(self,number):
        return self.level_label[self._index(number)]    
    
    def get_range(self,number):
        return self.level_range[self._index(number)]    
    
    def calculate_relative_level(self,number):
        return self.get_range(number).calculate_relative_level(number)


    def _calculate_distance(self, number,desired_label):
        index = self._index(number)
        _range = self.level_range[ index ]
        relative_level = _range.calculate_relative_level(number)
        
        try:
            dist = self.level_label.index(str(desired_label)) - index
        except:
            raise ValueError('%s is not a legal level'%desired_label)
        
        if dist == 0:
            return 0
        elif dist > 0:
            return dist + 1 - relative_level
        else:
            return relative_level - dist