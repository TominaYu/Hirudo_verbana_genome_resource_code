#!/usr/bin/env python3

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Patch


BASE = Path(
    "project_compare_EU_Hverb/"
    "09_suppfig_chr_correspondence/v2_chr_only"
)

COORDS = (
    BASE
    / "02_tables/wcHirVerb1_vs_our_primary.1to1.coords.tsv"
)

EU_MAP = (
    BASE
    / "02_tables/wcHirVerb1_chromosome_accessions.tsv"
)

OUTDIR = BASE / "03_figures"
TABLEDIR = BASE / "02_tables"

OUTDIR.mkdir(parents=True, exist_ok=True)

VISUAL_MIN_BLOCK = 20_000
VISUAL_MIN_IDENTITY = 95.0

COLLINEAR_COLOR = "#4C78A8"
OPPOSITE_COLOR = "#F58518"
BAR_COLOR = "#C8C8C8"
TEXT_COLOR = "#222222"


def union_length(intervals):
    if not intervals:
        return 0

    intervals = sorted(
        (min(start, end), max(start, end))
        for start, end in intervals
    )

    merged = [list(intervals[0])]

    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])

    return sum(end - start + 1 for start, end in merged)


def chr_number(value):
    match = re.search(r"(\d+)$", value)

    if not match:
        raise ValueError(f"Could not parse chromosome number: {value}")

    return int(match.group(1))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


# ------------------------------------------------------------
# EU accession mapping
# ------------------------------------------------------------

eu_metadata = {}

with EU_MAP.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        eu_metadata[row["eu_chr"]] = {
            "number": int(row["eu_chr_number"]),
            "accession": row["accession"],
            "length": int(row["length_bp"]),
        }


# ------------------------------------------------------------
# NUCmer coordinates
# show-coords -rclTH columns:
# S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ COVR COVQ REF QRY
# ------------------------------------------------------------

blocks = []

with COORDS.open() as handle:
    for line_number, raw in enumerate(handle, start=1):
        if not raw.strip():
            continue

        fields = raw.rstrip("\n").split()

        if len(fields) < 13:
            raise ValueError(
                f"Expected at least 13 fields at line {line_number}, "
                f"found {len(fields)}"
            )

        s1 = int(fields[0])
        e1 = int(fields[1])
        s2 = int(fields[2])
        e2 = int(fields[3])
        len1 = int(fields[4])
        len2 = int(fields[5])
        identity = float(fields[6])
        ref_length = int(fields[7])
        query_length = int(fields[8])
        ref_id = fields[11]
        query_id = fields[12]

        if ref_id not in eu_metadata:
            raise ValueError(f"Unexpected reference ID: {ref_id}")

        if not re.fullmatch(r"Chr(?:0[1-9]|1[0-4])", query_id):
            raise ValueError(f"Unexpected query ID: {query_id}")

        same_orientation = (
            (e1 - s1) * (e2 - s2) >= 0
        )

        blocks.append({
            "eu_chr": ref_id,
            "our_chr": query_id,
            "eu_start": min(s1, e1),
            "eu_end": max(s1, e1),
            "our_start": min(s2, e2),
            "our_end": max(s2, e2),
            "eu_length": ref_length,
            "our_length": query_length,
            "aligned_bp": min(len1, len2),
            "identity_pct": identity,
            "same_orientation": same_orientation,
        })


# ------------------------------------------------------------
# Pairwise statistics using all one-to-one blocks
# ------------------------------------------------------------

pair_blocks = defaultdict(list)

for block in blocks:
    pair_blocks[(block["eu_chr"], block["our_chr"])].append(block)

pair_rows = []

for (eu_chr, our_chr), rows in pair_blocks.items():
    aligned_bp = sum(row["aligned_bp"] for row in rows)

    same_bp = sum(
        row["aligned_bp"]
        for row in rows
        if row["same_orientation"]
    )

    reverse_bp = aligned_bp - same_bp

    eu_covered = union_length([
        (row["eu_start"], row["eu_end"])
        for row in rows
    ])

    our_covered = union_length([
        (row["our_start"], row["our_end"])
        for row in rows
    ])

    eu_length = rows[0]["eu_length"]
    our_length = rows[0]["our_length"]

    pair_rows.append({
        "eu_chr": eu_chr,
        "eu_chr_number": eu_metadata[eu_chr]["number"],
        "eu_accession": eu_metadata[eu_chr]["accession"],
        "our_chr": our_chr,
        "aligned_bp": aligned_bp,
        "block_count": len(rows),
        "same_orientation_bp": same_bp,
        "reverse_orientation_bp": reverse_bp,
        "dominant_orientation": (
            "same" if same_bp >= reverse_bp else "reverse"
        ),
        "eu_length": eu_length,
        "our_length": our_length,
        "eu_union_covered_bp": eu_covered,
        "our_union_covered_bp": our_covered,
        "eu_union_coverage_pct": 100 * eu_covered / eu_length,
        "our_union_coverage_pct": 100 * our_covered / our_length,
    })

pair_rows.sort(
    key=lambda row: (
        row["eu_chr_number"],
        chr_number(row["our_chr"]),
    )
)

write_tsv(
    TABLEDIR / "all_chromosome_pair_statistics.tsv",
    pair_rows,
    list(pair_rows[0].keys()),
)


# ------------------------------------------------------------
# Reciprocal-best pair assignment
# ------------------------------------------------------------

best_our_for_eu = {}

for eu_chr in eu_metadata:
    candidates = [
        row for row in pair_rows
        if row["eu_chr"] == eu_chr
    ]

    if not candidates:
        raise SystemExit(f"No alignment found for {eu_chr}")

    best_our_for_eu[eu_chr] = max(
        candidates,
        key=lambda row: row["aligned_bp"],
    )


our_ids = [f"Chr{x:02d}" for x in range(1, 15)]
best_eu_for_our = {}

for our_chr in our_ids:
    candidates = [
        row for row in pair_rows
        if row["our_chr"] == our_chr
    ]

    if not candidates:
        raise SystemExit(f"No alignment found for {our_chr}")

    best_eu_for_our[our_chr] = max(
        candidates,
        key=lambda row: row["aligned_bp"],
    )


reciprocal_pairs = []

for eu_chr, row in best_our_for_eu.items():
    our_chr = row["our_chr"]
    reciprocal = (
        best_eu_for_our[our_chr]["eu_chr"] == eu_chr
    )

    if reciprocal:
        reciprocal_pairs.append(row)


if len(reciprocal_pairs) != 14:
    raise SystemExit(
        "Expected 14 reciprocal-best chromosome pairs, "
        f"found {len(reciprocal_pairs)}. "
        "Inspect all_chromosome_pair_statistics.tsv."
    )


# Order by our Chr01–Chr14 for manuscript consistency
reciprocal_pairs.sort(
    key=lambda row: chr_number(row["our_chr"])
)

write_tsv(
    TABLEDIR / "reciprocal_best_chromosome_pairs.tsv",
    reciprocal_pairs,
    list(reciprocal_pairs[0].keys()),
)


# ------------------------------------------------------------
# Plotting subset
# ------------------------------------------------------------

visual_blocks = [
    block for block in blocks
    if block["aligned_bp"] >= VISUAL_MIN_BLOCK
    and block["identity_pct"] >= VISUAL_MIN_IDENTITY
]

write_tsv(
    TABLEDIR / "visualized_one_to_one_blocks.tsv",
    visual_blocks,
    list(visual_blocks[0].keys()),
)


# ------------------------------------------------------------
# Panel A: correspondence heatmap
# ------------------------------------------------------------

eu_order = [f"EU_chr{x:02d}" for x in range(1, 15)]
our_order = [f"Chr{x:02d}" for x in range(1, 15)]

matrix = [
    [
        sum(
            block["aligned_bp"]
            for block in blocks
            if block["eu_chr"] == eu_chr
            and block["our_chr"] == our_chr
        ) / 1_000_000
        for our_chr in our_order
    ]
    for eu_chr in eu_order
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

fig_h, ax_h = plt.subplots(figsize=(7.2, 6.4))

image = ax_h.imshow(
    matrix,
    cmap="Blues",
    aspect="equal",
)

ax_h.set_xticks(range(14))
ax_h.set_xticklabels(
    our_order,
    rotation=45,
    ha="right",
    fontsize=8,
)

ax_h.set_yticks(range(14))
ax_h.set_yticklabels(
    [f"wcHirVerb1 chr{x}" for x in range(1, 15)],
    fontsize=8,
)

ax_h.set_xlabel("Our primary assembly")
ax_h.set_ylabel("wcHirVerb1 assembly")

colorbar = fig_h.colorbar(image, ax=ax_h, fraction=0.045, pad=0.04)
colorbar.set_label("One-to-one aligned sequence (Mb)")

ax_h.set_xticks(
    [x - 0.5 for x in range(1, 14)],
    minor=True,
)
ax_h.set_yticks(
    [y - 0.5 for y in range(1, 14)],
    minor=True,
)

ax_h.grid(
    which="minor",
    color="white",
    linewidth=0.7,
)

ax_h.tick_params(which="minor", bottom=False, left=False)

fig_h.tight_layout()

for extension, kwargs in [
    ("pdf", {}),
    ("svg", {}),
    ("png", {"dpi": 600}),
]:
    fig_h.savefig(
        OUTDIR / f"SuppFig_chr_correspondence_A_heatmap.{extension}",
        bbox_inches="tight",
        **kwargs,
    )

plt.close(fig_h)


# ------------------------------------------------------------
# Panel B: 14 pairwise chromosome plots
# ------------------------------------------------------------

fig, axes = plt.subplots(
    7,
    2,
    figsize=(14, 17),
)

axes = axes.flatten()

for ax, pair in zip(axes, reciprocal_pairs):
    eu_chr = pair["eu_chr"]
    our_chr = pair["our_chr"]
    accession = pair["eu_accession"]
    eu_length = int(pair["eu_length"])
    our_length = int(pair["our_length"])

    flip_query = (
        pair["dominant_orientation"] == "reverse"
    )

    panel_blocks = [
        block for block in visual_blocks
        if block["eu_chr"] == eu_chr
        and block["our_chr"] == our_chr
    ]

    maximum_length = max(eu_length, our_length)

    ax.set_xlim(0, maximum_length)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ref_y = 0.73
    query_y = 0.16
    bar_height = 0.075

    ax.add_patch(
        Rectangle(
            (0, ref_y),
            eu_length,
            bar_height,
            facecolor=BAR_COLOR,
            edgecolor="none",
            zorder=3,
        )
    )

    ax.add_patch(
        Rectangle(
            (0, query_y),
            our_length,
            bar_height,
            facecolor=BAR_COLOR,
            edgecolor="none",
            zorder=3,
        )
    )

    # Draw large blocks first
    panel_blocks.sort(
        key=lambda row: row["aligned_bp"],
        reverse=True,
    )

    for block in panel_blocks:
        ref_start = block["eu_start"]
        ref_end = block["eu_end"]

        query_start = block["our_start"]
        query_end = block["our_end"]

        if flip_query:
            displayed_query_start = (
                our_length - query_end + 1
            )
            displayed_query_end = (
                our_length - query_start + 1
            )
        else:
            displayed_query_start = query_start
            displayed_query_end = query_end

        collinear_after_display = (
            block["same_orientation"] != flip_query
        )

        color = (
            COLLINEAR_COLOR
            if collinear_after_display
            else OPPOSITE_COLOR
        )

        if collinear_after_display:
            polygon_points = [
                (ref_start, ref_y),
                (ref_end, ref_y),
                (displayed_query_end, query_y + bar_height),
                (displayed_query_start, query_y + bar_height),
            ]
        else:
            polygon_points = [
                (ref_start, ref_y),
                (ref_end, ref_y),
                (displayed_query_start, query_y + bar_height),
                (displayed_query_end, query_y + bar_height),
            ]

        ax.add_patch(
            Polygon(
                polygon_points,
                closed=True,
                facecolor=color,
                edgecolor="none",
                alpha=0.48,
                zorder=1,
            )
        )

    eu_number = pair["eu_chr_number"]

    ax.text(
        0,
        0.97,
        (
            f"wcHirVerb1 chr{eu_number} "
            f"({accession}; {eu_length / 1e6:.2f} Mb)"
        ),
        ha="left",
        va="top",
        fontsize=8.5,
        fontweight="bold",
    )

    query_label = (
        f"Our {our_chr} ({our_length / 1e6:.2f} Mb)"
    )

    if flip_query:
        query_label += " — reversed for display"

    ax.text(
        0,
        0.02,
        query_label,
        ha="left",
        va="bottom",
        fontsize=8.5,
        fontweight="bold",
    )

    ax.text(
        maximum_length,
        0.97,
        (
            f"{pair['aligned_bp'] / 1e6:.2f} Mb aligned"
        ),
        ha="right",
        va="top",
        fontsize=7,
        color="#555555",
    )


legend_handles = [
    Patch(
        facecolor=COLLINEAR_COLOR,
        edgecolor="none",
        alpha=0.55,
        label="dominant collinear orientation",
    ),
    Patch(
        facecolor=OPPOSITE_COLOR,
        edgecolor="none",
        alpha=0.55,
        label="opposite local orientation",
    ),
]

fig.legend(
    handles=legend_handles,
    loc="upper center",
    ncol=2,
    frameon=False,
    fontsize=8,
    bbox_to_anchor=(0.5, 0.995),
)

fig.suptitle(
    "Chromosome-scale correspondence between wcHirVerb1 and the curated primary assembly",
    fontsize=13,
    fontweight="bold",
    y=0.999,
)

fig.text(
    0.5,
    0.006,
    (
        "NUCmer one-to-one blocks (delta-filter -1); "
        f"displayed blocks ≥{VISUAL_MIN_BLOCK / 1000:.0f} kb "
        f"and ≥{VISUAL_MIN_IDENTITY:.0f}% identity. "
        "Chromosome orientation was standardized per pair for visualization."
    ),
    ha="center",
    va="bottom",
    fontsize=7,
    color="#555555",
)

fig.subplots_adjust(
    left=0.045,
    right=0.985,
    top=0.965,
    bottom=0.025,
    hspace=0.36,
    wspace=0.12,
)

for extension, kwargs in [
    ("pdf", {}),
    ("svg", {}),
    ("png", {"dpi": 600}),
]:
    fig.savefig(
        OUTDIR / f"SuppFig_chr_correspondence_B_pairs.{extension}",
        bbox_inches="tight",
        **kwargs,
    )

plt.close(fig)

print(
    TABLEDIR
    / "reciprocal_best_chromosome_pairs.tsv"
)

print(
    OUTDIR
    / "SuppFig_chr_correspondence_A_heatmap.pdf"
)

print(
    OUTDIR
    / "SuppFig_chr_correspondence_B_pairs.pdf"
)
