from .SequencePattern import SequencePattern
from .DnaNotationPattern import DnaNotationPattern
from .EnzymeSitePattern import EnzymeSitePattern
from .HomopolymerPattern import HomopolymerPattern
from .RepeatedKmerPattern import RepeatedKmerPattern
from .PSSMPattern import PSSMPattern

SequencePattern.registered_string_pattern_classes = [
    HomopolymerPattern,
    RepeatedKmerPattern,
    EnzymeSitePattern,
    DnaNotationPattern,
]

__all__ = [
    'SequencePattern',
    'DnaNotationPattern',
    'EnzymeSitePattern',
    'HomopolymerPattern',
    'RepeatedKmerPattern',
    'PSSMPattern'
]