"""Generate the PNGs of all annotations

Run with:

>>> python3 generate_annotations.py
"""
import dna_features_viewer as dfv
import pandas
import matplotlib.pyplot as plt
import tqdm


def parse_locations(s):
    if "," not in s:
        start, end = map(int, s.strip().strip("+").split("-"))
        return [(start, end, 1 if s.endswith("+") else 0)]
    return [
        location
        for sub in s.split(",")
        for location in parse_locations(sub)
    ]


def graphic_features_from_row(row):
    if row.spec.startswith("@"):
        color = "#4b8bb2"
        fontdict = dict(size=14, color="white", weight="bold", family="Ubuntu")
    else:
        color = "#ededa7"
        fontdict = dict(size=14, color="black", family="Ubuntu")
    

    return [
        dfv.GraphicFeature(
            start=start,
            end=end,
            strand=strand,
            label=row.spec,
            thickness=23,
            box_color=None,
            color=color,
            fontdict=fontdict,
        )
        for start, end, strand in parse_locations(row.location)
    ]


def plot_row_feature(row):

    grecord = dfv.GraphicRecord(
        sequence=row.sequence,
        features=graphic_features_from_row(row),
        labels_spacing=0,
    )

    ax, bla = grecord.plot(annotate_inline=True, with_ruler=False)
    grecord.plot_sequence(ax, background=None)

    if row.highlights != "none":
        for (start, end, _) in parse_locations(row.highlights):
            ax.fill_between(
                [start - 0.5, end - 0.5], -0.9, -0.5, facecolor="r", alpha=0.15
            )
    ax.set_ylim(top=0.7)
    ax.figure.set_figheight(1)
    return ax


data = pandas.read_csv("./examples.csv")
for i, row in tqdm.tqdm(list(data.iterrows())):
    ax = plot_row_feature(row)
    ax.figure.savefig(
        "./%s.png" % row.filename, bbox_inches="tight", dpi=120
    )
    plt.close(ax.figure)
