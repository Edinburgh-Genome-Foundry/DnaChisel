import dna_features_viewer as dfv
import pandas
import matplotlib.pyplot as plt


def parse_locations(s):
    return [tuple(map(int, sub.strip().split("-"))) for sub in s.split(",")]


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
            strand=1,
            label=row.spec,
            thickness=23,
            box_color=None,
            color=color,
            fontdict=fontdict,
        )
        for start, end in parse_locations(row.location)
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
        for (start, end) in parse_locations(row.highlights):
            ax.fill_between(
                [start - 1.5, end - 0.5], -0.9, -0.5, facecolor="r", alpha=0.15
            )
    ax.set_ylim(top=0.7)
    ax.figure.set_figheight(1)
    return ax


data = pandas.read_csv("./examples.csv")
for i, row in data.iterrows():
    ax = plot_row_feature(row)
    ax.figure.savefig(
        "./%s.png" % row.filename, bbox_inches="tight", dpi=200
    )
    plt.close(ax.figure)
