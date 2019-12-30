import sys
from random import choice, randint

from .tools import aa2codon_table, analyzeCodons, codon2aa_table, randomMutation


def mutateCDS(sequence,
              keep_aa,
              mutableCodonsPosition,
              cds_region,
              pos=None,
              n_mut=[1, 2]):
    if keep_aa == True:
        result = analyzeCodons(sequence, mutableCodonsPosition)

        n_mutations = choice(n_mut)

        codons = (result[0])
        codons_ind = list(range(0, codons.__len__()))

        mutated = False
        while codons_ind.__len__() != 0 and n_mutations > 0:
            rnd_ind = codons_ind.pop(randint(0, codons_ind.__len__() - 1))
            rnd_codon = codons[rnd_ind]
            alt_codons = [
                c for c in aa2codon_table[codon2aa_table[rnd_codon]]
                if c != rnd_codon and codon2aa_table[c] != 'stop'
            ]
            if alt_codons.__len__() != 0:
                mutated = True
                n_mutations -= 1
                new_codon = choice(alt_codons)
                real_codon_pos = mutableCodonsPosition[rnd_ind]
                codon_position = (real_codon_pos - cds_region[0]) / 3
                all_codons = analyzeCodons(
                    sequence, list(range(cds_region[0], cds_region[1] + 1, 3)))[0]
                all_codons[codon_position] = new_codon
                new_seq = sequence[:cds_region[0]] + ''.join(
                    c for c in all_codons) + sequence[cds_region[1] + 1:]
                sequence = new_seq

        if mutated == False:
            sys.stderr.write(
                "RandomMutator: Not able to mutate sequence keeping AA\n")
            return None
        else:
            return new_seq


def mutateAll(sequence,
              keep_aa,
              mutable_region,
              cds_region,
              pos=None,
              n_mut=[1, 2]):

    #####
    # Not necessary to keep AA
    #
    n_mutations = choice(n_mut)

    if mutable_region == []:
        return None

    while n_mutations != 0:
        n_mutations -= 1
        if pos != None:
            intersect_mut = list(set(mutable_region) & set(pos))
            if intersect_mut != []:
                position_to_mutate = choice(intersect_mut)
            else:
                position_to_mutate = choice(mutable_region)
        else:
            position_to_mutate = choice(mutable_region)
        mutation = randomMutation(sequence[position_to_mutate])

        new_seq = sequence[:position_to_mutate] + mutation + sequence[
            position_to_mutate + 1:]
        sequence = new_seq

    return new_seq


def randomMutationOperator(sequence,
                           keep_aa,
                           mutable_region,
                           cds_region,
                           pos=None,
                           n_mut=[1, 2]):
    '''
        Operator that given a sequence, mutates the sequence randomly
            sequnce: sequence   
            mutable_region - a list with all bases that can be mutated
            cds_regions - a list of pairs with begin and end of CDSs - example: [(0,100), (150,300)]             
    '''
    # print(mutable_region)
    # print(cds_region)
    mutableCodonsPosition = [
        c for c in range(cds_region[0], cds_region[1], 3)
        if set([c, c + 1, c + 2]).issubset(mutable_region)
    ]
    mutableUTRPosition = list(
        set(mutable_region) - set(range(cds_region[0], cds_region[1])))

    if mutableCodonsPosition == [] and mutableUTRPosition == []:
        sys.stderr.write(
            "randomMutationOperator: No codons available for mutation\n")
        return None
    else:
        if keep_aa == True:
            if (mutableUTRPosition == []) or (mutableCodonsPosition != []
                                              and choice([True, False])):
                return mutateCDS(sequence, keep_aa, mutableCodonsPosition,
                                 cds_region, pos, n_mut)
            else:
                return mutateAll(sequence, keep_aa, mutableUTRPosition,
                                 cds_region, pos, n_mut)
        else:
            return mutateAll(sequence, keep_aa, mutable_region, cds_region,
                             pos, n_mut)
