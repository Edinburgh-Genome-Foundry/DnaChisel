"""Generic methods for grouping locations and sets of indices"""


def windows_overlap(window1, window2):
    """Return the overlap span between two windows.

    Parameters
    ----------

    window1, window2
      Each window is a couple of the form (start, end) indicating the range of
      a segment of integers.

    Returns
    -------

    None
      In case the two windows do not overlap.
    [start, end]
      The coordinates of the overlap segment if there is one.
    """
    start1, end1 = window1
    start2, end2 = window2

    if start2 < start1:
        return windows_overlap(window2, window1)

    if start1 <= start2 <= end1:
        return [start2, min(end1, end2)]
    else:
        return None


def subdivide_window(window, max_span):
    """Subdivide a window (start, end) into windows of size < max_span
    (start, i_1), (i_1, i_2), ... (i_n, end)"""
    start, end = window
    inds = list(range(start, end, max_span)) + [end]
    return list(zip(inds, inds[1:]))


def group_nearby_indices(indices, max_gap=None, max_group_spread=None):
    """Return a list of groups of the different indices.

    Indices are considered from smaller to larger and placed into groups

    Parameters
    ----------
    max_gap
      Maximal allowed difference between two consecutive numbers of a group

    max_group_spread
      Maximal allowed difference between the smallest and largest elements
      of a group.
    """
    if len(indices) == 0:
        return []
    indices = sorted(indices)
    current_group = [indices[0]]
    groups = [current_group]
    for ind in indices[1:]:
        gap_small_enough = (max_gap is None) or (
            ind - current_group[-1] < max_gap
        )
        spread_small_enough = (max_group_spread is None) or (
            ind - current_group[0] < max_group_spread
        )
        if gap_small_enough and spread_small_enough:
            current_group.append(ind)
        else:
            current_group = [ind]
            groups.append(current_group)
    return groups


def group_nearby_segments(segments, max_start_gap=None, max_start_spread=None):
    """Return a list of groups of the different indices.

    Indices are considered from smaller to larger and placed into groups

    Parameters
    ----------
    max_gap
      Maximal allowed difference between two consecutive numbers of a group

    max_group_spread
      Maximal allowed difference between the smallest and largest elements
      of a group.
    """
    if len(segments) == 0:
        return []
    segments = sorted(segments)
    current_group = [segments[0]]
    groups = [current_group]
    for seg in segments[1:]:
        gap_small_enough = (max_start_gap is None) or (
            seg[0] - current_group[-1][0] < max_start_gap
        )
        spread_small_enough = (max_start_spread is None) or (
            seg[0] - current_group[0][0] < max_start_spread
        )
        if gap_small_enough and spread_small_enough:
            current_group.append(seg)
        else:
            current_group = [seg]
            groups.append(current_group)
    return groups
