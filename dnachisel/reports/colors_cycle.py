import itertools

try:
    import matplotlib.cm as cm

    MATPLOTLIB_AVAILABLE = True
except:
    MATPLOTLIB_AVAILABLE = False


def colors_cycle(lightness_factor=1.0, color_shift=0):
    """Returns an iterator over a range of colors"""
    if MATPLOTLIB_AVAILABLE:
        cycle = itertools.cycle(
            [cm.Paired(color_shift + 0.21 + 0.21 * i % 1.0) for i in range(30)]
        )
        return (
            "#%02x%02x%02x"
            % tuple([int(255 * c * lightness_factor) for c in rgb_tuple[:3]])
            for rgb_tuple in cycle
        )
    else:
        return itertools.cycle(
            [
                "#f1cccc",
                "#f1e5cc",
                "#e3f1cc",
                "#ccf1e3",
                "#ccd7f1",
                "#e0ccf1",
                "f1cce7",
            ]
        )
