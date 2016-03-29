import re
from Bio.Restriction.Restriction_Dictionary import rest_dict

class DNAPattern:

    def __init__(self, pattern, size=None, name=None):
        if size is None:
            size = len(pattern)
        self.pattern = pattern
        self.size = size
        self.name = name

    def find_matches(self, sequence, window=None):
        if window is not None:
            wstart, wend = window
            subsequence = sequence[wstart:wend]
            return [
                (wstart + start, wstart + end)
                for (start, end) in self.find_matches(subsequence)
            ]
        return [
            (match.start(), match.end())
            for match in re.finditer(self.pattern, sequence)
        ]

    def __str__(self):
        return self.pattern + ("" if self.name is None else
                               " (%s)" % self.name)


def homopolymer_pattern(letter, number):
    return DNAPattern(number * letter)


def enzyme_pattern(enzyme_name):
    enzyme_site = rest_dict[enzyme_name]["site"]
    return DNAPattern(enzyme_site, name=enzyme_name)
