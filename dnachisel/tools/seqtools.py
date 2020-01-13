def analyze_ntcontent(seq):

    seq = seq.replace('u', 't')
    nuc_freq = {
        'NucleotideContentAT': 0,
        'NucleotideContentGC': 0,
        'NucleotideContentA': 0,
        'NucleotideContentT': 0,
        'NucleotideContentG': 0,
        'NucleotideContentC': 0
    }

    for i in range(len(seq)):
        nuc_freq['NucleotideContent' + seq[i].upper()] += 1

    nuc_freq['NucleotideContentA'] /= float(seq.__len__())
    nuc_freq['NucleotideContentT'] /= float(seq.__len__())
    nuc_freq['NucleotideContentG'] /= float(seq.__len__())
    nuc_freq['NucleotideContentC'] /= float(seq.__len__())
    nuc_freq['NucleotideContentAT'] = (nuc_freq['NucleotideContentA'] +
                                       nuc_freq['NucleotideContentT'])
    nuc_freq['NucleotideContentGC'] = (nuc_freq['NucleotideContentG'] +
                                       nuc_freq['NucleotideContentC'])

    return nuc_freq