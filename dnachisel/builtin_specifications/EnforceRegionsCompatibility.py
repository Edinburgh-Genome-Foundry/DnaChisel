from collections import Counter
import itertools

from ..Specification import Specification
from .VoidSpecification import VoidSpecification
from ..SpecEvaluation import SpecEvaluation

class EnforceRegionsCompatibility(Specification):
    max_possible_score = 0

    def __init__(self, regions_locations_locations, compatibility_condition,
                 label='', boost=1.0):
        self.regions_locations = regions_locations
        self.compatibility_condition = compatibility_condition
        self.boost = boost

    def evaluate(self, problem):
        incompatible_regions_locations_pairs = []
        for (r1, r2) in itertools.combinations(self.regions_locations, 2):
            if not self.compatibility_condition(r1, r2, problem):
                incompatible_regions_locations_pairs.append((r1, r2))

        all_regions_locations_with_incompatibility = [
            region
            for incompatibles_pair in incompatible_regions_locations_pairs
            for region in incompatibles_pair
        ]
        counter = Counter(all_regions_locations_with_incompatibility)
        all_regions_locations_with_incompatibility = sorted(
            list(set(all_regions_locations_with_incompatibility)),
            key=counter.get
        )

        score = -len(incompatible_regions_locations_pairs)
        if score == 0:
            message = "All compatible !"
        else:
            message = "Found the following imcompatibilities: %s" % (
                incompatible_regions_locations_pairs
            )
        return SpecEvaluation(
            self, problem,
            score=score,
            locations=all_regions_locations_with_incompatibility,
            message=message
        )

    def localized(self, location):
        if any(location.overlap_region(rl) for rl in self.regions_locations):
            return self
        else:
            return VoidSpecification(parent=self)

    def __repr__(self):
        return "CompatRegions(%s%s...)" % str(self.label,
                                              self.regions_locations[0])
    def __str__(self):
        return "CompatRegions(%s%s...)" % str(self.label,
                                              self.regions_locations[0])
