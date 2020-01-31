"""Ensure that the location is primer-frienly.

AllowPrimer is a good example SpecificationSet that will create specifications
of type RepeatedKmerPattern, UniquifyAllKmers, AvoidPattern,
EnforceMeltingTemperature, AvoidHeterodimerization.
"""

from ..Specification.SpecificationSet import SpecificationSet
from ..Location import Location
from ..SequencePattern import RepeatedKmerPattern
from .UniquifyAllKmers import UniquifyAllKmers
from .AvoidPattern import AvoidPattern
from .EnforceMeltingTemperature import EnforceMeltingTemperature
from .AvoidHeterodimerization import AvoidHeterodimerization


class AllowPrimer(SpecificationSet):
    """Enforce various specifications for enabling primers at the location.

    This is useful for making sure that you will be able to conduct a PCR or
    sanger sequencing with a primer annealing at a particular location of
    the sequence.

    Shorthand for annotations: "primer".

    Parameters
    ----------

    location
      The exact location where the primer will anneal, i.e. the subsequence
      under this location will be the sequence

    tmin, tmax
      Minimum and maximum acceptable melting temperatures in Celcius, for
      instance 55 and 70. When the used as an optimization objective the
      "target" will be (tmin + tmax)/2.
    max_homology_length
      Maximal length of any homology between the subsequence at that location
      and anywhere else in the whole sequence.

    avoid_heterodim_with
      List of ATGC strings representing external (primer) sequences with which
      the optimised location should have no annealing.

    max_heterodim_tm
      Max melting temperature of the heterodimerization between this subsegment
      and the other sequences in "avoid_heterodim_with"

    avoided_repeats
      Properties of the repeated patterns avoided. List of pairs (K, N)
      meaning "avoid K-mers repeated N times in a row".
    """

    shorthand_name = "primer"

    def __init__(
        self,
        location=None,
        tmin=50,
        tmax=70,
        max_homology_length=6,
        avoid_heterodim_with=None,
        max_heterodim_tm=5,
        avoided_repeats=((2, 5), (3, 4), (4, 3)),
    ):
        location = Location.from_data(location)
        specs = {
            "unique_sequence": UniquifyAllKmers(
                k=max_homology_length, location=location
            ),
            "melting_temperature": EnforceMeltingTemperature(
                mini=tmin, maxi=tmax, location=location
            ),
            **{
                "repeats_%d_%d"
                % (k, n): AvoidPattern(
                    RepeatedKmerPattern(k, n), location=location
                )
                for (k, n) in avoided_repeats
            },
        }
        if avoid_heterodim_with is not None:
            specs["avoid_heterodimerization"] = AvoidHeterodimerization(
                other_primers_sequences=avoid_heterodim_with,
                tmax=max_heterodim_tm,
                location=location,
            )
        self.register_specifications(specs)
