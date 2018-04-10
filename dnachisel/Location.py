"""Implements the Location class.

The class has useful methods for finding overlaps between locations, extract
a subsequence from a sequence, etc.
"""

from .biotools import reverse_complement
from Bio.SeqFeature import SeqFeature, FeatureLocation

class Location:
    """Represent a segment of a sequence, with a start, end, and strand.

    Warning: we use Python's splicing notation, so Location(5, 10) represents
    sequence[5, 6, 7, 8, 9] which corresponds to nucleotides number
    6, 7, 8, 9, 10.

    The data structure is similar to a Biopython's FeatureLocation, but with
    different methods for different purposes.

    Parameters
    ----------
    start
      Lowest position index of the segment

    end
      Highest position index of the segment

    strand
      Either 1 or -1 for sense or anti-sense orientation.

    """
    __slots__ = ['strand', 'start', 'end']

    def __init__(self, start, end, strand=None):
        """Initialize."""
        self.start = start
        self.end = end
        self.strand = strand

    def overlap_region(self, other_location):
        """Return the overlap span between two locations (None if None)."""
        strand = self.strand
        if other_location.start < self.start:
            self, other_location = other_location, self

        if other_location.start < self.end:
            start = other_location.start
            end = min(self.end, other_location.end)
            strand = self.strand
            return Location(start, end, strand)
        else:
            return None

    def extended(self, extension_length, lower_limit=0, upper_limit=None,
                 left=True, right=True):
        """Extend the location of a few basepairs on each side."""

        if left:
            lower = max(lower_limit, self.start - extension_length)
        else:
            lower = self.start

        if right:
            upper = self.end + extension_length
            if upper_limit is not None:
                upper = min(upper_limit, upper)
        else:
            upper = self.end

        return Location(lower, upper, self.strand)

    def extract_sequence(self, sequence):
        """Return the subsequence read at the given location."""
        result = sequence[self.start:self.end]
        if self.strand == -1:
            return reverse_complement(result)
        else:
            return result

    def to_tuple(self):
        """Return (start, end, strand)."""
        return (self.start, self.end, self.strand)

    @property
    def indices(self):
        result = list(range(self.start, self.end))
        return result if (self.strand != -1) else result[::-1]

    @staticmethod
    def from_tuple(some_tuple, default_strand=0):
        if len(some_tuple) == 2:
            start, end = some_tuple
            strand = default_strand
        else:
            start, end, strand = some_tuple
        return Location(start, end, strand)

    def __geq__(self, other):
        """Greater than."""
        return self.to_tuple() >= other.to_tuple()

    def __lt__(self, other):
        """Lower than."""
        return self.to_tuple() < other.to_tuple()

    def __add__(self, number):
        """Return the location shifted by the number"""
        return Location(self.start + number, self.end + number, self.strand)

    def __repr__(self):
        """Represent"""
        result = "%d-%d" % (self.start, self.end)
        if self.strand is not None:
            result += {1: "(+)", -1: "(-)", 0: ""}[self.strand]
        return result

    def __len__(self):
        """Size of the location"""
        return self.end - self.start

    @staticmethod
    def merge_overlapping_locations(locations):
        """Return a list of locations obtained by mergin all overlapping."""
        if len(locations) == 0:
            return locations
        locations = sorted(locations)
        new_locations = [locations[0]]
        for loc in locations[1:]:
            if new_locations[-1].overlap_region(loc) is not None:
                new_locations[-1].end = max(new_locations[-1].end, loc.end)
            else:
                new_locations.append(loc)
        return new_locations

    @staticmethod
    def from_biopython_location(location):
        """Return a DnaChisel Location from a Biopython location."""
        start, end, strand = [
            None if e is None else int(e)
            for e in [location.start, location.end, location.strand]
        ]
        return Location(start, end, strand)

    def to_biopython_location(self):
        """Return a Biopython FeatureLocation equivalent to the location."""
        start, end, strand = [
            None if e is None else int(e)
            for e in [self.start, self.end, self.strand]
        ]
        return FeatureLocation(start, end, strand)

    def to_biopython_feature(self, feature_type="misc_feature", **qualifiers):
        """Return a Biopython SeqFeature with same location and custom
        qualifiers."""
        return SeqFeature(self.to_biopython_location(),
                          type=feature_type,
                          qualifiers=qualifiers)


class MultiLocation:
    "TODO: Implement when necessary"
    def __init__(self, locations):
        self.locations = locations
