#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


BASE = Path(
    "project_primary_FASTA/"
    "12_suppfig2_read_support"
)

STRICT_FILE = (
    BASE
    / "02_tables/raw_snapshot/strict_spanning_by_manifest.tsv"
)

GAP_FILE = (
    BASE
    / "02_tables/raw_snapshot/strict_spanning_GH1_gap_sweep.tsv"
)

OUTDIR = BASE / "04_charts"
TABLEDIR = BASE / "02_tables/figure_data"

OUTDIR.mkdir(parents=True, exist_ok=True)
TABLEDIR.mkdir(parents=True, exist_ok=True)

FLANK = 1000
MIN_MAPQ = 20
PRIMARY_ONLY = "True"

ADOPTED_COLOR = "#0072B2"
ALTERNATIVE_COLOR = "#B8B8B8"
NO_GAP_COLOR = "#0072B2"
GAP_COLOR = "#B8B8B8"
TEXT_COLOR = "#222222"


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_unique(rows, test, junction):
    matches = [
        row
        for row in rows
        if row["test"] == test
        and row["junction"] == junction
        and int(row["flank"]) == FLANK
        and int(row["min_mapq"]) == MIN_MAPQ
        and row["primary_only"] == PRIMARY_ONLY
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Expected one row but found {len(matches)}: "
            f"{test}, {junction}"
        )

    return int(matches[0]["spanning_reads"])


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


strict_rows = read_tsv(STRICT_FILE)
gap_rows = read_tsv(GAP_FILE)

# ------------------------------------------------------------------
# Panel B:
# test03 is the selected ptg000020l-containing no-gap topology.
# test04 provides the corresponding ptg000010l-containing comparison
# while retaining ptg000020l downstream.
# ------------------------------------------------------------------

panel_b_rows = [
    {
        "comparison_group": "Testing ptg000010l",
        "display_label": (
            "ptg000027l–\n"
            "ptg000015l-derived"
        ),
        "model_role": "adopted",
        "test": "test03_G_noGap",
        "junction": "J2_Chr14_A2_to_Chr14_E",
        "spanning_reads": select_unique(
            strict_rows,
            "test03_G_noGap",
            "J2_Chr14_A2_to_Chr14_E",
        ),
    },
    {
        "comparison_group": "Testing ptg000010l",
        "display_label": (
            "ptg000027l–\n"
            "ptg000010l"
        ),
        "model_role": "alternative",
        "test": "test04_BG_noGap",
        "junction": "J2_Chr14_A2_to_Chr14_B",
        "spanning_reads": select_unique(
            strict_rows,
            "test04_BG_noGap",
            "J2_Chr14_A2_to_Chr14_B",
        ),
    },
    {
        "comparison_group": "Testing ptg000010l",
        "display_label": (
            "ptg000010l–\n"
            "ptg000015l-derived"
        ),
        "model_role": "alternative",
        "test": "test04_BG_noGap",
        "junction": "J3_Chr14_B_to_Chr14_E",
        "spanning_reads": select_unique(
            strict_rows,
            "test04_BG_noGap",
            "J3_Chr14_B_to_Chr14_E",
        ),
    },
    {
        "comparison_group": "Testing ptg000020l",
        "display_label": (
            "ptg000015l-derived–\n"
            "ptg000003l"
        ),
        "model_role": "alternative",
        "test": "test01_direct_noGap",
        "junction": "J3_Chr14_E_to_Chr14_H1",
        "spanning_reads": select_unique(
            strict_rows,
            "test01_direct_noGap",
            "J3_Chr14_E_to_Chr14_H1",
        ),
    },
    {
        "comparison_group": "Testing ptg000020l",
        "display_label": (
            "ptg000015l-derived–\n"
            "ptg000020l"
        ),
        "model_role": "adopted",
        "test": "test03_G_noGap",
        "junction": "J3_Chr14_E_to_Chr14_G",
        "spanning_reads": select_unique(
            strict_rows,
            "test03_G_noGap",
            "J3_Chr14_E_to_Chr14_G",
        ),
    },
    {
        "comparison_group": "Testing ptg000020l",
        "display_label": (
            "ptg000020l–\n"
            "ptg000003l"
        ),
        "model_role": "adopted",
        "test": "test03_G_noGap",
        "junction": "J4_Chr14_G_to_Chr14_H1",
        "spanning_reads": select_unique(
            strict_rows,
            "test03_G_noGap",
            "J4_Chr14_G_to_Chr14_H1",
        ),
    },
]

write_tsv(
    TABLEDIR / "SuppFig2B_strict_spanning_data.tsv",
    panel_b_rows,
    [
        "comparison_group",
        "display_label",
        "model_role",
        "test",
        "junction",
        "spanning_reads",
    ],
)

# ------------------------------------------------------------------
# Panel C:
# For inserted-gap models, both flanks of the N block must yield the
# same strict-spanning count. Stop if they differ.
# ------------------------------------------------------------------

def gap_support(test, left_junction, right_junction):
    left_value = select_unique(
        gap_rows,
        test,
        left_junction,
    )

    right_value = select_unique(
        gap_rows,
        test,
        right_junction,
    )

    if left_value != right_value:
        raise ValueError(
            f"Gap-edge counts differ for {test}: "
            f"{left_value} versus {right_value}"
        )

    return left_value


panel_c_rows = [
    {
        "gap_kb": 0,
        "display_label": "No gap",
        "test": "test09_GH1_noGap",
        "spanning_reads": select_unique(
            gap_rows,
            "test09_GH1_noGap",
            "J4_Chr14_G_to_Chr14_H1",
        ),
    },
    {
        "gap_kb": 1,
        "display_label": "1 kb",
        "test": "test10_GH1_gap1k",
        "spanning_reads": gap_support(
            "test10_GH1_gap1k",
            "J4_Chr14_G_to_GAP",
            "J5_GAP_to_Chr14_H1",
        ),
    },
    {
        "gap_kb": 5,
        "display_label": "5 kb",
        "test": "test11_GH1_gap5k",
        "spanning_reads": gap_support(
            "test11_GH1_gap5k",
            "J4_Chr14_G_to_GAP",
            "J5_GAP_to_Chr14_H1",
        ),
    },
    {
        "gap_kb": 10,
        "display_label": "10 kb",
        "test": "test12_GH1_gap10k",
        "spanning_reads": gap_support(
            "test12_GH1_gap10k",
            "J4_Chr14_G_to_GAP",
            "J5_GAP_to_Chr14_H1",
        ),
    },
]

write_tsv(
    TABLEDIR / "SuppFig2C_GH1_gap_sweep_data.tsv",
    panel_c_rows,
    [
        "gap_kb",
        "display_label",
        "test",
        "spanning_reads",
    ],
)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#555555")
    ax.spines["bottom"].set_color("#555555")
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(
        axis="both",
        width=0.6,
        color="#555555",
    )
    ax.grid(
        axis="x",
        linewidth=0.4,
        color="#D9D9D9",
        zorder=0,
    )


def add_bar_values(ax, bars, values, offset=0.4):
    for bar, value in zip(bars, values):
        ax.text(
            value + offset,
            bar.get_y() + bar.get_height() / 2,
            str(value),
            ha="left",
            va="center",
            fontsize=8,
            color=TEXT_COLOR,
        )


# ------------------------------------------------------------------
# Plot Panel B
# ------------------------------------------------------------------

fig_b, axes_b = plt.subplots(
    1,
    2,
    figsize=(8.2, 3.6),
    sharex=True,
)

for ax, group in zip(
    axes_b,
    ["Testing ptg000010l", "Testing ptg000020l"],
):
    subset = [
        row for row in panel_b_rows
        if row["comparison_group"] == group
    ]

    labels = [row["display_label"] for row in subset]
    values = [row["spanning_reads"] for row in subset]
    colors = [
        (
            ADOPTED_COLOR
            if row["model_role"] == "adopted"
            else ALTERNATIVE_COLOR
        )
        for row in subset
    ]

    y_positions = list(range(len(subset)))

    bars = ax.barh(
        y_positions,
        values,
        color=colors,
        edgecolor="none",
        height=0.62,
        zorder=2,
    )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_title(group, fontsize=9, fontweight="bold")
    ax.set_xlim(0, 26)

    style_axis(ax)
    add_bar_values(ax, bars, values)

axes_b[0].set_ylabel("Candidate junction")
axes_b[0].set_xlabel("Strictly spanning HiFi reads")
axes_b[1].set_xlabel("Strictly spanning HiFi reads")

fig_b.text(
    0.5,
    0.012,
    "Exact boundary; primary alignments only; MAPQ ≥20; 1-kb flank criterion",
    ha="center",
    va="bottom",
    fontsize=7,
    color="#555555",
)

fig_b.subplots_adjust(
    left=0.20,
    right=0.98,
    top=0.88,
    bottom=0.25,
    wspace=0.55,
)

for suffix, kwargs in [
    ("pdf", {}),
    ("svg", {}),
    ("png", {"dpi": 600}),
]:
    fig_b.savefig(
        OUTDIR / f"SuppFig2B_strict_spanning.{suffix}",
        bbox_inches="tight",
        **kwargs,
    )

plt.close(fig_b)

# ------------------------------------------------------------------
# Plot Panel C
# ------------------------------------------------------------------

fig_c, ax_c = plt.subplots(figsize=(4.2, 3.6))

labels = [row["display_label"] for row in panel_c_rows]
values = [row["spanning_reads"] for row in panel_c_rows]

colors = [
    NO_GAP_COLOR,
    GAP_COLOR,
    GAP_COLOR,
    GAP_COLOR,
]

bars = ax_c.bar(
    labels,
    values,
    color=colors,
    edgecolor="none",
    width=0.66,
    zorder=2,
)

for bar, value in zip(bars, values):
    ax_c.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.45,
        str(value),
        ha="center",
        va="bottom",
        fontsize=8,
        color=TEXT_COLOR,
    )

ax_c.set_ylabel("Strictly spanning HiFi reads")
ax_c.set_xlabel("Inserted gap between ptg000020l and ptg000003l")
ax_c.set_ylim(0, 26)
ax_c.set_title(
    "Gap-length comparison",
    fontsize=9,
    fontweight="bold",
)

ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.spines["left"].set_color("#555555")
ax_c.spines["bottom"].set_color("#555555")
ax_c.spines["left"].set_linewidth(0.6)
ax_c.spines["bottom"].set_linewidth(0.6)
ax_c.tick_params(
    axis="both",
    width=0.6,
    color="#555555",
)
ax_c.grid(
    axis="y",
    linewidth=0.4,
    color="#D9D9D9",
    zorder=0,
)

fig_c.text(
    0.5,
    0.012,
    "Exact boundary; primary alignments only; MAPQ ≥20; 1-kb flank criterion",
    ha="center",
    va="bottom",
    fontsize=7,
    color="#555555",
)

fig_c.subplots_adjust(
    left=0.18,
    right=0.97,
    top=0.86,
    bottom=0.28,
)

for suffix, kwargs in [
    ("pdf", {}),
    ("svg", {}),
    ("png", {"dpi": 600}),
]:
    fig_c.savefig(
        OUTDIR / f"SuppFig2C_GH1_gap_sweep.{suffix}",
        bbox_inches="tight",
        **kwargs,
    )

plt.close(fig_c)

print("Panel B data:")
for row in panel_b_rows:
    print(
        row["comparison_group"],
        row["display_label"].replace("\n", " "),
        row["spanning_reads"],
    )

print("\nPanel C data:")
for row in panel_c_rows:
    print(
        row["display_label"],
        row["spanning_reads"],
    )

print(f"\nWrote charts to: {OUTDIR}")
print(f"Wrote source data to: {TABLEDIR}")
