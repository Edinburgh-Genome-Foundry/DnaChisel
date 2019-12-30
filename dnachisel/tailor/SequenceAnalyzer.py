'''
Created on Nov 1, 2012

@author: jcg
'''

import sys
from .Solution import Solution
from csv import DictReader

class SequenceAnalyzer(object):
    
    '''
    Initializes class that analyzes sequence features 
    '''
    
    def __init__(self, input_file, input_type, sep = ","):
        
        if input_type == "CSV":
            self.list_of_input_sequences = self.readCSV(input_file,sep)
        elif input_type == "FASTA": 
            self.list_of_input_sequences = self.readFASTA(input_file)
        # elif input_file == "GENBANK":
        #     self.list_of_input_sequences = self.readGENBANK(input_file)
        else:
            sys.stderr.write("The input type entered is not supported, please use one of the following: [CSV,FASTA,GENBANK]") 
    
    def readCSV(self,input_file,sep=","):
        pass
        list_seq = []
        
        reader = DictReader(open(input_file), delimiter=sep, quotechar='"')
        
        for l in reader:            
            list_seq.append(l)
            
        return list_seq

    def readFASTA(self,input_file):
        pass
        list_seq = []
        
        reader = open(input_file)
        
        name = ""
        seq  = ""
        
        for l in reader:
            l=l.rstrip()
            
            if l[0] == ">":
                name = l.split(' ')[0][1:]
            else:
                seq = l                            
                list_seq.append({ 'name' : name , 'sequence' : seq})
            
        return list_seq
    
    def readGENBANK(self):
        pass        
            
    def configureSolution(self, solution):
        pass
    
    def outputStart(self):
        pass
    
    def output(self,solution):
        pass
        
    def run(self):
        
        self.outputStart()
        
        for sequence in self.list_of_input_sequences:
            sol_id = sequence['name']
            seq = sequence['sequence']
            
            solution = Solution(sol_id = sol_id, sequence = seq)
            self.configureSolution(solution)
            
            self.output(solution)
            