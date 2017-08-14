from ..Specification import Specification
from ..SpecEvaluation import SpecEvaluation


class VoidSpecification(Specification):
    """Void Specifications are a special Specifications that always pass.

    Void Specifications are generally obtained when a spec is "voided"
    by a localization. For instance, if we are optimizing the segment (10,50)
    of a DNA segment, the Specification EnforceTranslation([300,500]) does not
    apply as it concerns a Gene that is in a completely different segment.
    Therefore the localized version of EnforceTranslation will be void.

    Note: the initializer accepts starred arguments/keyword arguments to make
    it easy to void any other specification by replacing the class to Void.
    Particularly useful when importing an optimization problem from genbank.
    """
    best_possible_score = 0
    __name__ = "Void"

    def __init__(self, parent_specification=None, message='default', boost=1.0,
                 *args, **kwargs):
        """Initialize."""
        self.parent_specification = parent_specification
        if message == 'default':
            message = "Pass (not relevant in this context)"
        self.message = message
        self.boost = boost

    def evaluate(self, problem):
        """Pass the test, always."""
        return SpecEvaluation(self, problem, score=0,
                              message=self.message,
                              locations=None)

    def label_parameters(self):
        return [str(self.parent_specification)]
