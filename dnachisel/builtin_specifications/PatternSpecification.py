from ..Specification import Specification
from .VoidSpecification import VoidSpecification
from ..SequencePattern import enzyme_pattern, DnaNotationPattern

class PatternSpecification(Specification):
    """Class for Specifications such as presence or absence of a pattern.

    The particularity of the PatternSpecifications is that they will either
    infer or ask for the length of the associated pattern and use this to
    localize the specification efficiently when performing local optimization or
    solving.

    Parameters
    ----------

    pattern
      A SequencePattern or DnaNotationPattern

    enzyme
      Enzyme name, can be provided instead of pattern or dna_pattern

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)``


    """
    shrink_when_localized = True
    priority = 1 # higher than normal

    def __init__(self, pattern=None, location=None, boost=1.0, enzyme=None):
        """Initialize."""
        # if dna_pattern is not None:
        #     pattern = dna_pattern
        if enzyme is not None:
            pattern = enzyme_pattern(enzyme)
        if isinstance(pattern, str):
            pattern = DnaNotationPattern(pattern)
        self.pattern = pattern
        self.location = location
        self.enzyme = enzyme
        self.boost = boost

    def localized(self, location):
        """Localize the pattern to the given location. Taking into account the
        specification's own location, and the size of the pattern."""
        pattern_size = self.pattern.size
        if self.location is None:
            new_location = location.extended(pattern_size - 1)
        else:
            if self.location.overlap_region(location) is None:
                return VoidSpecification(parent_specification=self)
            else:
                if not self.shrink_when_localized:
                    return self
                extended_location = location.extended(pattern_size - 1)
                new_location = self.location.overlap_region(extended_location)

        return self.copy_with_changes(location=new_location)

    def label_parameters(self):
        return [('enzyme'. self.enzyme) if (self.enzyme is not None)
                  else (self.pattern.sequence
                        if hasattr(self.pattern, 'sequence')
                        else str(self.pattern))]
