"""Text and number formatting operations"""
from copy import deepcopy
import json

import numpy as np


def round_all_numbers_in_dict(d, rounding_digits=2, outplace=True):
    """ Return a new version of dict d with all floats rounded to N digits."""
    if outplace:
        d = deepcopy(d)
    for k, v in d.items():
        if isinstance(v, float):
            d[k] = np.round(v, rounding_digits)
        if isinstance(v, dict):
            round_all_numbers_in_dict(v, rounding_digits, outplace=False)
    return d


def dict_to_pretty_string(d, rounding_digits=2, indent=2):
    """Return a nicely JSON-like formatted string to print a dict."""
    d = round_all_numbers_in_dict(d, rounding_digits)
    formatted_text = json.dumps(d, indent=indent)
    for char in '{}",':
        formatted_text = formatted_text.replace(char, "")
    return formatted_text

def score_to_formatted_string(score, characters=9):
    """Transform a number (score) into a best-format string.

    The format will be either int (2234), float (10.234) or engineering
    (1.20E5), whichever is shorter. The score is then padded with left
    whitespaces to obtained the desired number of ``characters``."""
    raw = str(int(score) if (int(score) == score) else score)
    as_float = "%.02f" % score
    as_eng = "%.02E." % score
    return min([raw, as_float, as_eng], key=len).rjust(characters)