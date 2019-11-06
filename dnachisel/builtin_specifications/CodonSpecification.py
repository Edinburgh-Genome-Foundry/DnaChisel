"""Implements core specification VoidSpecification."""
from ..Specification import Specification

# from .VoidSpecification import VoidSpecification
from ..Location import Location


class CodonSpecification(Specification):
    """Special class for dealing with codon.

    In particular, this class implements a specific localized method

    """

    def codon_index_to_location(self, index):
        if self.location.strand >= 0:
            return Location(
                start=self.location.start + 3 * index,
                end=self.location.start + 3 * (index + 1),
                strand=1
            )
        else:
            return Location(
                start=self.location.end - 3 * (index + 1),
                end=self.location.end - 3 * index,
                strand=-1,
            )

    def localized(self, location, problem=None, with_righthand=True):
        """Generic localization method for codon specifications.

        Calls the class'  ``.localized_on_window`` method at the end.

        """
        if self.location is not None:
            overlap = self.location.overlap_region(location)
            if overlap is None:
                return None
            else:
                # return self
                o_start, o_end = overlap.start, overlap.end
                w_start, w_end = self.location.start, self.location.end

                if self.location.strand != -1:
                    start_codon = int((o_start - w_start) / 3)
                    end_codon = int((o_end - w_start - 1) / 3) + 1
                    new_location = Location(
                        start=w_start + 3 * start_codon,
                        end=min(w_end, w_start + 3 * (end_codon)),
                        strand=self.location.strand,
                    )
                else:
                    start_codon = int((w_end - o_end) / 3)
                    end_codon = int((w_end - o_start - 1) / 3) + 1
                    new_location = Location(
                        start=max(w_start, w_end - 3 * (end_codon)),
                        end=w_end - 3 * start_codon,
                        strand=self.location.strand,
                    )
                return self.localized_on_window(
                    new_location, start_codon, end_codon
                )
        else:
            return self
