from .biotools import reverse_complement
from Bio.SeqFeature import FeatureLocation

class Location:

    def __init__(self, start, end, strand=None):
        self.start = start
        self.end = end
        self.strand = strand

    def overlap_region(self, other_location):
        """Return the overlap span between two locations (None if None).
        """

        if other_location.start < self.start:
            self, other_location = other_location, self

        if self.start <= other_location.start <= self.end:
            start = other_location.start
            end = min(self.end, other_location.start)
            strand = self.strand
            return Location(start, end, strand)
        else:
            return None

    def extended(self, extension_length, upper_limit=None):
        lower = max(0, self.start - extension_length)
        upper = self.end + extension_length
        if upper_limit is not None:
            upper = min(upper_limit, upper)
        return Location(lower, upper, self.strand)

    def extract_sequence(self, sequence):
        result = sequence[self.start:self.end]
        if self.strand == -1:
            return reverse_complement(result)
        else:
            return result

    def to_tuple(self):
        return (self.start, self.end, self.strand)

    def __geq__(self, other):
        return self.to_tuple() >= other.to_tuple()

    def __add__(self, number):
        return Location(self.start + number, self.end + number, self.strand)

    def __repr__(self):
        result = "%d-%d" % (self.start, self.end)
        if self.strand is not None:
            result += "(%s)" % ({1: "+", -1: "-"}[self.strand])
        return result

    def __len__(self):
        return self.end - self.start

    @staticmethod
    def from_biopython_location(location):
        return Location(location.start, location.end, location.strand)

    def to_biopython_location(self):
        return FeatureLocation(self.start, self.end, self.strand)


class MultiLocation:

    def __init__(self, locations):
        self.locations = locations
