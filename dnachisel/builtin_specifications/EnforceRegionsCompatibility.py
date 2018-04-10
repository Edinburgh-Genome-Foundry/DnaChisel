from collections import Counter
import itertools

from ..Specification import Specification
from ..Location import Location
from .VoidSpecification import VoidSpecification
from ..SpecEvaluation import SpecEvaluation

class EnforceRegionsCompatibility(Specification):
    max_possible_score = 0

    def __init__(self, locations, compatibility_condition,
                 condition_label='', boost=1.0):
        self.locations = [
            Location.from_tuple(location)
            for location in locations
        ]
        self.compatibility_condition = compatibility_condition
        self.condition_label = condition_label
        self.boost = boost

    def evaluate(self, problem):
        incompatible_locations_pairs = []
        for (r1, r2) in itertools.combinations(self.locations, 2):
            if not self.compatibility_condition(r1, r2, problem):
                incompatible_locations_pairs.append((r1, r2))

        all_locations_with_incompatibility = [
            region
            for incompatibles_pair in incompatible_locations_pairs
            for region in incompatibles_pair
        ]
        counter = Counter(all_locations_with_incompatibility)
        all_locations_with_incompatibility = sorted(
            list(set(all_locations_with_incompatibility)),
            key=counter.get
        )

        score = -len(incompatible_locations_pairs)
        if score == 0:
            message = "All compatible !"
        else:
            message = "Found the following imcompatibilities: %s" % (
                incompatible_locations_pairs
            )
        return SpecEvaluation(
            self, problem,
            score=score,
            locations=all_locations_with_incompatibility,
            message=message
        )

    def localized(self, location, problem=None):
        if any(location.overlap_region(rl) for rl in self.locations):
            return self
        else:
            return VoidSpecification(parent=self)

    def __repr__(self):
        return "CompatRegions(%s%s...)" % (self.condition_label,
                                           self.locations[0])
    def __str__(self):
        return "CompatRegions(%s%s...)" % (self.condition_label,
                                           self.locations[0])

    def label_parameters(self):
        params = [('locations', ", ".join([str(l) for l in self.locations]))]
        if self.condition_label not in [None, '']:
            params.append(('condition', self.condition_label))
        return params
