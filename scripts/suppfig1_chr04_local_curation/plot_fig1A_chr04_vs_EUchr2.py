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
from matplotlib.lines import Line2D


CHR04_LENGTH = 15_289_010
EU_CHR2_LENGTH = 16_479_789
TARGET_ID = "OZ237780.1"

COMPONENTS = [
    {
        "label": "A1",
        "query_id": "Chr14_A1",
        "start": 1,
        "end": 6_435_254,
        "color": "#0072B2",
    },
    {
        "label": "A2",
        "query_id": "Chr14_A2",
        "start": 6_435_255,
        "end": 7_500_755,
        "color": "#56B4E9",
    },
    {
        "label": "E",
        "query_id": "Chr14_E",
        "start": 7_500_756,
        "end": 9_013_138,
        "color": "#009E73",
    },
    {
        "label": "G",
        "query_id": "Chr14_G",
        "start": 9_013_139,
        "end": 9_382_960,
        "color": "#E69F00",
    },
    {
        "label": "H1",
        "query_id": "Chr14_H1",
        "start": 9_382_961,
        "end": 12_299_010,
        "color": "#CC79A7",
    },
    {
        "label": "H3",
        "query_id": "Chr14_H3",
        "start": 12_299_011,
        "end": 15_289_010,
        "color": "#D55E00",
    },
]

COMPONENT_BY_ID = {x["query_id"]: x for x in COMPONENTS}


def read_blocks(path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_tsv(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def to_float(row, key):
    return float(row[key])


def to_int(row, key):
    return int(float(row[key]))


def add_ribbon(
    ax,
    query_start,
    query_end,
    target_start,
    target_end,
    strand,
    color,
    y_query=0.79,
    y_target=0.25,
):
    """
    Draw a curved alignment ribbon.

    query_start/query_end and target_start/target_end are in Mb.
    Reverse-orientation blocks are drawn as twisted ribbons with
    a dashed outline.
    """
    middle_y = (y_query + y_target) / 2.0

    if strand == "+":
        lower_left = target_start
        lower_right = target_end
        edgecolor = "none"
        linestyle = "-"
        linewidth = 0.0
        alpha = 0.22
    else:
        lower_left = target_end
        lower_right = target_start
        edgecolor = "#4D4D4D"
        linestyle = "--"
        linewidth = 0.35
        alpha = 0.16

    vertices = [
        (query_start, y_query),
        (query_start, middle_y),
        (lower_left, middle_y),
        (lower_left, y_target),

        (lower_right, y_target),

        (lower_right, middle_y),
        (query_end, middle_y),
        (query_end, y_query),

        (query_start, y_query),
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

    patch = PathPatch(
        MplPath(vertices, codes),
        facecolor=color,
        edgecolor=edgecolor,
        linestyle=linestyle,
        linewidth=linewidth,
        alpha=alpha,
        zorder=1,
    )
    ax.add_patch(patch)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--blocks",
        required=True,
        help="fig1A_backbone6_nucmer_blocks_all.tsv",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=20_000,
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=95.0,
    )
    parser.add_argument(
        "--title",
        default="Final Chr04 versus wcHirVerb1 chromosome 2",
    )
    args = parser.parse_args()

    input_path = Path(args.blocks)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    rows = read_blocks(input_path)

    selected = []

    for row in rows:
        query_id = row["query_id"]

        if query_id not in COMPONENT_BY_ID:
            continue

        if to_int(row, "minimum_aligned_length") < args.min_length:
            continue

        if to_float(row, "identity_pct") < args.min_identity:
            continue

        selected.append(row)

    if not selected:
        raise SystemExit("No alignment blocks passed the selected filters.")

    selected.sort(
        key=lambda row: to_int(row, "minimum_aligned_length"),
        reverse=True,
    )

    output_fields = list(selected[0].keys())
    selected_tsv = Path(f"{out_prefix}.selected_blocks.tsv")
    write_tsv(selected_tsv, selected, output_fields)

    block_counts = Counter(row["label"] for row in selected)
    forward_counts = Counter(
        row["label"]
        for row in selected
        if row["final_chr04_strand_vs_target"] == "+"
    )
    reverse_counts = Counter(
        row["label"]
        for row in selected
        if row["final_chr04_strand_vs_target"] == "-"
    )

    summary_tsv = Path(f"{out_prefix}.summary.tsv")

    with summary_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "label",
            "total_blocks",
            "forward_blocks",
            "reverse_blocks",
            "minimum_block_length",
            "minimum_identity_pct",
        ])

        for component in COMPONENTS:
            label = component["label"]
            writer.writerow([
                label,
                block_counts[label],
                forward_counts[label],
                reverse_counts[label],
                args.min_length,
                args.min_identity,
            ])

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    fig, ax = plt.subplots(figsize=(11.6, 4.3))

    query_y = 0.82
    target_y = 0.18
    bar_height = 0.055

    # Draw alignment ribbons first
    for row in selected:
        component = COMPONENT_BY_ID[row["query_id"]]

        add_ribbon(
            ax=ax,
            query_start=to_int(row, "chr04_start_1based") / 1_000_000,
            query_end=to_int(row, "chr04_end_1based") / 1_000_000,
            target_start=to_int(row, "target_start_1based") / 1_000_000,
            target_end=to_int(row, "target_end_1based") / 1_000_000,
            strand=row["final_chr04_strand_vs_target"],
            color=component["color"],
            y_query=query_y,
            y_target=target_y + bar_height,
        )

    # Final Chr04 component bar
    for component in COMPONENTS:
        start_mb = (component["start"] - 1) / 1_000_000
        end_mb = component["end"] / 1_000_000
        width_mb = end_mb - start_mb

        rect = Rectangle(
            (start_mb, query_y),
            width_mb,
            bar_height,
            facecolor=component["color"],
            edgecolor="white",
            linewidth=0.7,
            zorder=5,
        )
        ax.add_patch(rect)

        center = (start_mb + end_mb) / 2.0

        if component["label"] == "G":
            ax.text(
                center,
                query_y + bar_height / 2,
                "G",
                ha="center",
                va="center",
                fontsize=7,
                fontweight="bold",
                color="black",
                zorder=6,
            )
        else:
            ax.text(
                center,
                query_y + bar_height / 2,
                component["label"],
                ha="center",
                va="center",
                fontsize=7.5,
                fontweight="bold",
                color="white",
                zorder=6,
            )

        ax.text(
            center,
            query_y + bar_height + 0.027,
            f"{component['label']} ({width_mb:.2f} Mb)",
            ha="center",
            va="bottom",
            fontsize=6.5,
            color="#333333",
            zorder=6,
        )

    # EU chromosome 2 bar
    eu_length_mb = EU_CHR2_LENGTH / 1_000_000

    eu_bar = Rectangle(
        (0, target_y),
        eu_length_mb,
        bar_height,
        facecolor="#BDBDBD",
        edgecolor="#4D4D4D",
        linewidth=0.65,
        zorder=5,
    )
    ax.add_patch(eu_bar)

    ax.text(
        eu_length_mb / 2,
        target_y + bar_height / 2,
        "wcHirVerb1 chromosome 2",
        ha="center",
        va="center",
        fontsize=8,
        fontweight="bold",
        color="black",
        zorder=6,
    )

    # Sequence labels
    ax.text(
        0,
        query_y + 0.145,
        "Final Chr04: A1–A2–E–G–H1–H3 (15.289 Mb)",
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )

    ax.text(
        0,
        target_y - 0.065,
        "OZ237780.1 (16.480 Mb)",
        ha="left",
        va="top",
        fontsize=8,
        color="#333333",
    )

    # Axes
    maximum_length_mb = max(
        CHR04_LENGTH,
        EU_CHR2_LENGTH,
    ) / 1_000_000

    ax.set_xlim(-0.15, maximum_length_mb + 0.15)
    ax.set_ylim(0.03, 1.04)

    ticks = list(range(0, 17, 2))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(x) for x in ticks], fontsize=8)
    ax.set_xlabel("Sequence coordinate (Mb)", fontsize=9)

    ax.set_yticks([])

    for spine in ["top", "left", "right"]:
        ax.spines[spine].set_visible(False)

    ax.spines["bottom"].set_color("#666666")
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(axis="x", width=0.6, length=3, color="#666666")

    # Title
    if args.title:
        ax.set_title(
            args.title,
            loc="left",
            fontsize=11,
            fontweight="bold",
            pad=10,
        )

    # Legend
    orientation_handles = [
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
            alpha=0.16,
            label="reverse orientation",
        ),
    ]

    ax.legend(
        handles=orientation_handles,
        loc="upper right",
        frameon=False,
        fontsize=7,
        ncol=2,
        bbox_to_anchor=(1.0, 1.025),
    )

    # Reproducibility note
    note = (
        "NUCmer one-to-one blocks (delta-filter -1); "
        f"displayed blocks ≥{args.min_length / 1000:.0f} kb and "
        f"≥{args.min_identity:.0f}% identity; "
        f"n={len(selected)} blocks. "
        "Filtering was applied for visualization only."
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
        top=0.87,
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
    print(f"Wrote: {summary_tsv}")
    print(f"Wrote: {out_prefix}.pdf")
    print(f"Wrote: {out_prefix}.svg")
    print(f"Wrote: {out_prefix}.png")


if __name__ == "__main__":
    main()
