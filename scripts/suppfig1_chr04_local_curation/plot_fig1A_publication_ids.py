#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch, Patch
from matplotlib.path import Path as MplPath


CHR04_LENGTH = 15_289_010
EU_CHR2_LENGTH = 16_479_789

COMPONENTS = [
    {
        "query_id": "Chr14_A1",
        "display_id": "ptg000025l",
        "full_source_id": "ptg000025l_length_6435254",
        "source_start": 1,
        "source_end": 6_435_254,
        "chr04_start": 1,
        "chr04_end": 6_435_254,
        "label_y": 1.015,
        "color": "#0072B2",
    },
    {
        "query_id": "Chr14_A2",
        "display_id": "ptg000027l",
        "full_source_id": "ptg000027l_length_1065501",
        "source_start": 1,
        "source_end": 1_065_501,
        "chr04_start": 6_435_255,
        "chr04_end": 7_500_755,
        "label_y": 1.055,
        "color": "#56B4E9",
    },
    {
        "query_id": "Chr14_E",
        "display_id": "ptg000015l\n57,001–1,569,383",
        "full_source_id": "ptg000015l_length_1569383",
        "source_start": 57_001,
        "source_end": 1_569_383,
        "chr04_start": 7_500_756,
        "chr04_end": 9_013_138,
        "label_y": 1.015,
        "color": "#009E73",
    },
    {
        "query_id": "Chr14_G",
        "display_id": "ptg000020l",
        "full_source_id": "ptg000020l_length_369822",
        "source_start": 1,
        "source_end": 369_822,
        "chr04_start": 9_013_139,
        "chr04_end": 9_382_960,
        "label_y": 1.105,
        "color": "#E69F00",
    },
    {
        "query_id": "Chr14_H1",
        "display_id": "ptg000003l",
        "full_source_id": "ptg000003l_length_2916050",
        "source_start": 1,
        "source_end": 2_916_050,
        "chr04_start": 9_382_961,
        "chr04_end": 12_299_010,
        "label_y": 1.025,
        "color": "#CC79A7",
    },
    {
        "query_id": "Chr14_H3",
        "display_id": "ptg000004l\n1–2,990,000",
        "full_source_id": "ptg000004l_length_3055668",
        "source_start": 1,
        "source_end": 2_990_000,
        "chr04_start": 12_299_011,
        "chr04_end": 15_289_010,
        "label_y": 1.015,
        "color": "#D55E00",
    },
]

COMPONENT_BY_QUERY = {
    component["query_id"]: component
    for component in COMPONENTS
}


def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def i(row, key):
    return int(float(row[key]))


def f(row, key):
    return float(row[key])


def add_ribbon(
    ax,
    query_start,
    query_end,
    target_start,
    target_end,
    strand,
    color,
    query_y,
    target_y,
):
    middle_y = (query_y + target_y) / 2

    if strand == "+":
        target_left = target_start
        target_right = target_end
        edgecolor = "none"
        linestyle = "-"
        linewidth = 0
        alpha = 0.22
    else:
        target_left = target_end
        target_right = target_start
        edgecolor = "#4D4D4D"
        linestyle = "--"
        linewidth = 0.7
        alpha = 0.28

    vertices = [
        (query_start, query_y),
        (query_start, middle_y),
        (target_left, middle_y),
        (target_left, target_y),

        (target_right, target_y),

        (target_right, middle_y),
        (query_end, middle_y),
        (query_end, query_y),

        (query_start, query_y),
    ]

    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,

        MplPath.LINETO,

        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,

        MplPath.CLOSEPOLY,
    ]

    ax.add_patch(
        PathPatch(
            MplPath(vertices, codes),
            facecolor=color,
            edgecolor=edgecolor,
            linestyle=linestyle,
            linewidth=linewidth,
            alpha=alpha,
            zorder=1,
        )
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--blocks", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--min-length", type=int, default=20_000)
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--title", default="")

    args = parser.parse_args()

    input_path = Path(args.blocks)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    rows = read_rows(input_path)

    selected = []

    for row in rows:
        if row["query_id"] not in COMPONENT_BY_QUERY:
            continue

        if i(row, "minimum_aligned_length") < args.min_length:
            continue

        if f(row, "identity_pct") < args.min_identity:
            continue

        selected.append(row)

    if not selected:
        raise SystemExit("No alignment blocks passed the filters.")

    # Draw larger blocks first, with smaller blocks remaining visible on top.
    selected.sort(
        key=lambda row: i(row, "minimum_aligned_length"),
        reverse=True,
    )

    selected_tsv = Path(f"{out_prefix}.selected_blocks.tsv")
    write_tsv(
        selected_tsv,
        selected,
        list(selected[0].keys()),
    )

    mapping_rows = []

    for component in COMPONENTS:
        mapping_rows.append({
            "figure_display_id": component["display_id"].replace("\n", " "),
            "full_source_record": component["full_source_id"],
            "source_start_1based": component["source_start"],
            "source_end_1based": component["source_end"],
            "final_record": "Chr04",
            "chr04_start_1based": component["chr04_start"],
            "chr04_end_1based": component["chr04_end"],
        })

    mapping_tsv = Path(f"{out_prefix}.source_segment_mapping.tsv")
    write_tsv(
        mapping_tsv,
        mapping_rows,
        list(mapping_rows[0].keys()),
    )

    block_counts = Counter(
        COMPONENT_BY_QUERY[row["query_id"]]["full_source_id"]
        for row in selected
    )

    forward_counts = Counter(
        COMPONENT_BY_QUERY[row["query_id"]]["full_source_id"]
        for row in selected
        if row["final_chr04_strand_vs_target"] == "+"
    )

    reverse_counts = Counter(
        COMPONENT_BY_QUERY[row["query_id"]]["full_source_id"]
        for row in selected
        if row["final_chr04_strand_vs_target"] == "-"
    )

    summary_rows = []

    for component in COMPONENTS:
        source_id = component["full_source_id"]

        summary_rows.append({
            "full_source_record": source_id,
            "total_blocks": block_counts[source_id],
            "same_orientation_blocks": forward_counts[source_id],
            "reverse_orientation_blocks": reverse_counts[source_id],
            "minimum_block_length": args.min_length,
            "minimum_identity_pct": args.min_identity,
        })

    summary_tsv = Path(f"{out_prefix}.summary.tsv")
    write_tsv(
        summary_tsv,
        summary_rows,
        list(summary_rows[0].keys()),
    )

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig, ax = plt.subplots(figsize=(11.8, 4.5))

    query_y = 0.78
    target_y = 0.17
    bar_height = 0.055

    # Alignment ribbons
    for row in selected:
        component = COMPONENT_BY_QUERY[row["query_id"]]

        add_ribbon(
            ax=ax,
            query_start=i(row, "chr04_start_1based") / 1_000_000,
            query_end=i(row, "chr04_end_1based") / 1_000_000,
            target_start=i(row, "target_start_1based") / 1_000_000,
            target_end=i(row, "target_end_1based") / 1_000_000,
            strand=row["final_chr04_strand_vs_target"],
            color=component["color"],
            query_y=query_y,
            target_y=target_y + bar_height,
        )

    # Curated Chr04 source-segment bar
    for component in COMPONENTS:
        start_mb = (component["chr04_start"] - 1) / 1_000_000
        end_mb = component["chr04_end"] / 1_000_000
        width_mb = end_mb - start_mb
        center = (start_mb + end_mb) / 2

        ax.add_patch(
            Rectangle(
                (start_mb, query_y),
                width_mb,
                bar_height,
                facecolor=component["color"],
                edgecolor="white",
                linewidth=0.7,
                zorder=5,
            )
        )

        # Publication label placed above the segment.
        ax.annotate(
            component["display_id"],
            xy=(center, query_y + bar_height),
            xytext=(center, component["label_y"]),
            ha="center",
            va="bottom",
            fontsize=6.6,
            color="#303030",
            linespacing=0.95,
            arrowprops={
                "arrowstyle": "-",
                "color": "#777777",
                "linewidth": 0.45,
                "shrinkA": 1,
                "shrinkB": 1,
            },
            zorder=7,
            annotation_clip=False,
        )

    ax.text(
        0,
        1.145,
        "Curated Chr04 (15.289 Mb; source-contig-derived segments shown above)",
        ha="left",
        va="bottom",
        fontsize=9.3,
        fontweight="bold",
    )

    # wcHirVerb1 chromosome 2 bar
    eu_length_mb = EU_CHR2_LENGTH / 1_000_000

    ax.add_patch(
        Rectangle(
            (0, target_y),
            eu_length_mb,
            bar_height,
            facecolor="#BDBDBD",
            edgecolor="#4D4D4D",
            linewidth=0.65,
            zorder=5,
        )
    )

    ax.text(
        eu_length_mb / 2,
        target_y + bar_height / 2,
        "wcHirVerb1 chromosome 2 (OZ237780.1; 16.480 Mb)",
        ha="center",
        va="center",
        fontsize=8,
        fontweight="bold",
        color="black",
        zorder=6,
    )

    maximum_length_mb = max(
        CHR04_LENGTH,
        EU_CHR2_LENGTH,
    ) / 1_000_000

    ax.set_xlim(-0.15, maximum_length_mb + 0.15)
    ax.set_ylim(0.02, 1.20)

    ticks = list(range(0, 17, 2))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(x) for x in ticks], fontsize=8)
    ax.set_xlabel("Sequence coordinate (Mb)", fontsize=9)

    ax.set_yticks([])

    for spine in ["top", "left", "right"]:
        ax.spines[spine].set_visible(False)

    ax.spines["bottom"].set_color("#666666")
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(
        axis="x",
        width=0.6,
        length=3,
        color="#666666",
    )

    if args.title:
        ax.set_title(
            args.title,
            loc="left",
            fontsize=11,
            fontweight="bold",
            pad=10,
        )

    legend_handles = [
        Patch(
            facecolor="#777777",
            edgecolor="none",
            alpha=0.22,
            label="same orientation",
        ),
        Patch(
            facecolor="#777777",
            edgecolor="#4D4D4D",
            linestyle="--",
            linewidth=0.7,
            alpha=0.28,
            label="reverse orientation",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        loc="upper right",
        frameon=False,
        fontsize=7,
        ncol=2,
        bbox_to_anchor=(1.0, 1.045),
    )

    note = (
        "NUCmer one-to-one blocks (delta-filter -1); "
        f"displayed blocks ≥{args.min_length / 1000:.0f} kb and "
        f"≥{args.min_identity:.0f}% identity; "
        f"n={len(selected)} blocks. "
        "Thresholds were applied for visualization only."
    )

    ax.text(
        0.5,
        -0.18,
        note,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7,
        color="#555555",
    )

    fig.subplots_adjust(
        left=0.055,
        right=0.985,
        top=0.90,
        bottom=0.22,
    )

    fig.savefig(
        f"{out_prefix}.pdf",
        bbox_inches="tight",
    )

    fig.savefig(
        f"{out_prefix}.svg",
        bbox_inches="tight",
    )

    fig.savefig(
        f"{out_prefix}.png",
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Selected blocks: {len(selected)}")
    print(f"Wrote: {selected_tsv}")
    print(f"Wrote: {mapping_tsv}")
    print(f"Wrote: {summary_tsv}")
    print(f"Wrote: {out_prefix}.pdf")
    print(f"Wrote: {out_prefix}.svg")
    print(f"Wrote: {out_prefix}.png")


if __name__ == "__main__":
    main()
