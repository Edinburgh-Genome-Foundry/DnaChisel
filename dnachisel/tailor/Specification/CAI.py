import sys
from math import exp, log
from random import choice, randint
from uuid import uuid4

from ..tools import (aa2codon_table, analyzeCodons, codon2aa_table,
                     default_cai_table,analyze_cai)
from .Specification import Specification
from ..Solution import Solution


def SimpleCAIOperator(sequence,
                      cai_range,
                      keep_aa,
                      mutable_region,
                      cds_regions,
                      direction='+',
                      cai_table=default_cai_table):
    '''
        Operator that given a sequence, mutates the sequence to change CAI
            sequnce: sequence 
            cai_range - start and end position to calculate cai - a tuple in the form (start, end)  
            mutable_region - a list with all bases that can be mutated
            cds_regions - a list of pairs with begin and end of CDSs - example: [(0,100), (150,300)]            
            direction: either increase ('+') or decrease ('-') CAI 
    '''
    mutated = False
    mutableCodonsPosition = [
        c for c in range(cai_range[0], cai_range[1], 3)
        if set([c, c + 1, c + 2]).issubset(mutable_region)
    ]

    if len(mutableCodonsPosition) == 0:
        sys.stderr.write(
            "SimpleCAIOperator: No codons available for mutation\n")
        return None

    result = analyzeCodons(sequence, mutableCodonsPosition)

    codons = (result[0])
    codons_cai = (result[1])
    codons_ind = list(range(0, codons.__len__()))

    while not mutated and codons_ind.__len__() != 0:

        rnd_ind = codons_ind.pop(randint(0, codons_ind.__len__() - 1))
        rnd_codon = codons[rnd_ind]
        rnd_codon_cai = codons_cai[rnd_ind]

        #select alternative codons
        if keep_aa == True and direction == '+':
            alt_codons = [
                c for c in aa2codon_table[codon2aa_table[rnd_codon]]
                if c != rnd_codon and cai_table[c] > rnd_codon_cai
                and codon2aa_table[c] != 'stop'
            ]
        elif keep_aa == True and direction == '-':
            alt_codons = [
                c for c in aa2codon_table[codon2aa_table[rnd_codon]]
                if c != rnd_codon and cai_table[c] < rnd_codon_cai
                and codon2aa_table[c] != 'stop'
            ]
        elif keep_aa == False and direction == '+':
            alt_codons = list(
                k for k, v in cai_table.iteritems()
                if v > rnd_codon_cai and codon2aa_table[k] != 'stop')
        elif keep_aa == False and direction == '-':
            alt_codons = list(
                k for k, v in cai_table.iteritems()
                if v < rnd_codon_cai and codon2aa_table[k] != 'stop')

        if alt_codons.__len__() != 0:
            mutated = True
            new_codon = choice(alt_codons)
            #print "new: " + str(new_codon)

    if mutated == False:
        sys.stderr.write("SimpleCAIOperator: Not able to mutate sequence\n")
        return None
    else:
        #print "CAI operator: old_codon -> " + str(rnd_codon)
        #print "CAI operator: new_codon -> " + str(new_codon)
        real_codon_pos = mutableCodonsPosition[rnd_ind]
        codon_position = (real_codon_pos - cai_range[0]) // 3
        all_codons = analyzeCodons(sequence,
                                   list(range(cai_range[0], cai_range[1] + 1, 3)))[0]
        all_codons[codon_position] = new_codon

        new_seq = sequence[:cai_range[0]] + ''.join(
            c for c in all_codons) + sequence[cai_range[1] + 1:]
        return new_seq


class CAI(Specification):
    """
    CAI Feature
        solution - solution where CAI should be computed
        label - some label to append to the name
        cai_range - start and end position to calculate CAI - a tuple in the form (start, end)  
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept
    """
    def __init__(
            self,

            solution=None,
            label="",
            args={
                'cai_range': (0, 59),
                'mutable_region': None,
                'cds_region': None,
                'keep_aa': True
            }):

        #General properties of feature
        Specification.__init__(self, solution=solution, label=label)
        #Specifics of this Feature
        self.cai_range = args['cai_range']
        self.sequence = solution.sequence[self.cai_range[0]:(
            self.cai_range[1] + 1)]
        self.mutable_region = args.get('mutable_region',solution.mutable_region)
        self.cds_region = args.get('cds_region', solution.cds_region)
        self.keep_aa = args.get('keep_aa',solution.keep_aa) 
        self.set_scores()
        self.set_level()


    def set_scores(self, scoring_function=analyze_cai):
        self.scores[self.label + "CAI"] = scoring_function(self.sequence)

    def mutate(self, operator=SimpleCAIOperator):
        if not self.targetInstructions:
            return None
        new_seq = operator(self.solution.sequence, self.cai_range,
                           self.keep_aa, self.mutable_region, self.cds_region,
                           self.targetInstructions['direction'])
        if not new_seq:
            return None
        return Solution(sol_id=str(uuid4().int),
                        sequence=new_seq,
                        cds_region=self.cds_region,
                        mutable_region=self.mutable_region,
                        parent=self.solution,
                        design=self.solution.designMethod)
