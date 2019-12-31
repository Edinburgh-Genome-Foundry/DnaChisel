from random import randint
from uuid import uuid4

from ..Solution import Solution
from ..tools import analyze_ntcontent, appendLabelToDict, codon2aa_table
from .Specification import Specification


def SimpleNtContentOperator(seq,
                            direction=0,
                            nucleotide=[],
                            mutable_region=None,
                            cds_region=(0, 9),
                            keep_aa=True):
    if direction == 0:
        return seq

    mutated = False
    seq = seq.lower()

    #check direction to decide possible mutations
    if direction == '-':
        #select mutable position based on presence of nucleotide(s)
        mutable_positions = [
            pos for pos in mutable_region if seq[pos] in set(nucleotide)
        ]
    elif direction == '+':
        #select mutable position based on absence of nucleotide(s)
        mutable_positions = [
            pos for pos in mutable_region if not (seq[pos] in set(nucleotide))
        ]

    while not mutated and mutable_positions.__len__() != 0:
        #print mutable_positions
        rnd_pos = mutable_positions.pop(
            randint(0,
                    mutable_positions.__len__() - 1))

        if direction == '-':
            possible_mutations = list(set('atcg') - set(nucleotide))
        elif direction == '+':
            possible_mutations = list(nucleotide)

        while not mutated and possible_mutations.__len__() != 0:

            new_seq = seq[:rnd_pos] + possible_mutations.pop(
                randint(0,
                        possible_mutations.__len__() - 1)) + seq[rnd_pos + 1:]

            if keep_aa == True and rnd_pos >= cds_region[
                    0] and rnd_pos <= cds_region[1]:
                #check if AA remains the same
                codon_p = (rnd_pos - cds_region[0]) / 3
                initial_codon = seq[(cds_region[0] +
                                     codon_p * 3):(cds_region[0] +
                                                   codon_p * 3 + 3)]
                final_codon = new_seq[(cds_region[0] +
                                       codon_p * 3):(cds_region[0] +
                                                     codon_p * 3 + 3)]

                #print "initial codon: " + str(initial_codon) + " AA: " + codon2aa_table[initial_codon]
                #print "final codon: " + str(final_codon) + " AA: " + codon2aa_table[final_codon]

                if codon2aa_table[initial_codon] == codon2aa_table[
                        final_codon]:
                    mutated = True
            else:
                mutated = True

    if mutated == False:
        return None

    return new_seq



class NucleotideContent(Specification):
    """
    Nucleotide Content Feature
        solution - solution where nucleotide content should be computed
        label - some label to append to the name
        hi_range - start and end position to calculate nucleotide content - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept        
    """
    def __init__(self,
     nucleotideContentObject = None, 
     solution=None, label="", 
     args = { 
         'ntcontent_range' : (0,9),  
         'mutable_region' : None,
         'cds_region' : None,
         'keep_aa' : True 
         }
    ):
        if nucleotideContentObject == None: #create new instance
            #General properties of feature
            Specification.__init__(self, solution=solution, label=label)
            #Specifics of this Feature
            self.ntcontent_range    = args['ntcontent_range']
            self.sequence           = solution.sequence[self.ntcontent_range[0]:self.ntcontent_range[1]+1]
            self.mutable_region = args.get('mutable_region',solution.mutable_region)
            self.cds_region = args.get('cds_region' , solution.cds_region)
            self.keep_aa = args.get('keep_aa',solution.keep_aa)
            self.set_scores()
            self.set_level()
        else:
            Specification.__init__(self, nucleotideContentObject)
            self.ntcontent_range    = nucleotideContentObject.ntcontent_range
            self.sequence           = nucleotideContentObject.sequence
            self.mutable_region     = nucleotideContentObject.mutable_region
            self.cds_region         = nucleotideContentObject.cds_region
            self.keep_aa            = nucleotideContentObject.keep_aa
            self.scores             = nucleotideContentObject.scores
    
    def set_scores(self, scoring_function = analyze_ntcontent):
        self.scores = appendLabelToDict(scoring_function(self.sequence), self.label)
            
    def mutate(self, operator=SimpleNtContentOperator):
        if not self.targetInstructions:
            return None
        new_seq = operator(self.solution.sequence, self.targetInstructions['direction'], self.nucleotides, self.mutable_region, self.cds_region, keep_aa=self.keep_aa)
        if not new_seq:
            return None                
        return Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_region = self.cds_region, mutable_region = list(self.mutable_region), parent=self.solution, design=self.solution.designMethod)

class NucleotideContentAT(NucleotideContent):
    """
    Check AT content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['a','t']
        self.set_level()
        
class NucleotideContentGC(NucleotideContent):
    """
    Check GC content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['g','c']
        self.set_level()
        
class NucleotideContentA(NucleotideContent):
    """
    Check A content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)                
        self.nucleotides = ['a']
        self.set_level()
        
class NucleotideContentT(NucleotideContent):
    """
    Check T content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['t']
        self.set_level()
        
class NucleotideContentG(NucleotideContent):
    """
    Check G content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['g']
        self.set_level()
        
class NucleotideContentC(NucleotideContent):
    """
    Check C content
    """
    def __init__(self,  nucleotideContentObject):
        NucleotideContent.__init__(self,nucleotideContentObject)
        self.nucleotides = ['c']
        self.set_level()
