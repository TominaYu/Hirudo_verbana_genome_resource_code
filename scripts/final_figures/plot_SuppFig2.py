#!/usr/bin/env python3
"""Render the code-generated base for Supplementary Figure 2.

Presentation changes after inspecting Supplementary_Figure2_v33_refresh.pdf:
- remove text/box collisions in panel a;
- keep Model 3 emphasis behind the topology rather than through the model label;
- wrap the two '-derived' block labels inside their boxes;
- enlarge the small topology schematics below panel b;
- read the validated support counts from the audited source TSVs rather than hard-coding them;
- retain the full 'Boundary-spanning primary HiFi alignments' terminology.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from common import COLORS, set_pub_style, save_triplet_fixed

CONTIG_COLORS = {
    "ptg000025l": "#CCB7E8",
    "ptg000027l": "#AFC9E8",
    "ptg000010l": "#F3B37B",
    "ptg000015l-derived": "#CBE2B3",
    "ptg000020l": "#9BD9DC",
    "ptg000003l": "#F4D99A",
    "ptg000004l-derived": "#EFB6CF",
    "N-spacer": "#EEEEEE",
}
MODELS = [
    ("Model 1", ["ptg000025l","ptg000027l","ptg000015l-derived","ptg000003l","ptg000004l-derived"]),
    ("Model 2", ["ptg000025l","ptg000027l","ptg000010l","ptg000015l-derived","ptg000003l","ptg000004l-derived"]),
    ("Model 3", ["ptg000025l","ptg000027l","ptg000015l-derived","ptg000020l","ptg000003l","ptg000004l-derived"]),
    ("Model 4", ["ptg000025l","ptg000027l","ptg000010l","ptg000015l-derived","ptg000020l","ptg000003l","ptg000004l-derived"]),
]

DISPLAY_TEXT = {
    "ptg000015l-derived": "ptg000015l-\nderived",
    "ptg000004l-derived": "ptg000004l-\nderived",
}

def draw_boxes(ax, items, x0, y, available=0.80, gap=0.008, h=0.115, fs=5.0):
    n = len(items)
    width = (available - gap * (n - 1)) / n
    x = x0
    bounds = []
    for item in items:
        patch = FancyBboxPatch(
            (x, y - h/2), width, h,
            boxstyle="round,pad=0.003,rounding_size=0.006",
            facecolor=CONTIG_COLORS[item], edgecolor="none",
        )
        ax.add_patch(patch)
        ax.text(
            x + width/2, y, DISPLAY_TEXT.get(item, item),
            ha="center", va="center", fontsize=fs, linespacing=0.90,
        )
        bounds.append((x, x + width))
        x += width + gap
    for (a0, a1), (b0, b1) in zip(bounds[:-1], bounds[1:]):
        ax.plot([a1, b0], [y, y], color="#222222", lw=0.65, zorder=0)
    return bounds

def load_counts(strict_spanning_path, gap_sweep_path):
    bpath = Path(strict_spanning_path)
    cpath = Path(gap_sweep_path)
    print("Using:", bpath)
    print("Using:", cpath)

    b = pd.read_csv(bpath, sep="\t")
    required = {"test","junction","spanning_reads"}
    if not required.issubset(b.columns):
        raise SystemExit(f"{bpath} lacks {sorted(required)}; observed {list(b.columns)}")

    q = b.copy()
    if "flank" in q.columns:
        q = q[pd.to_numeric(q["flank"], errors="coerce") == 1000]
    if "min_mapq" in q.columns:
        q = q[pd.to_numeric(q["min_mapq"], errors="coerce") == 20]

    def get(test, junction):
        m = q[(q["test"].astype(str) == test) & (q["junction"].astype(str) == junction)]
        if len(m) != 1:
            raise SystemExit(f"Expected exactly one row for {test} / {junction} in {bpath}; found {len(m)}")
        return int(m.iloc[0]["spanning_reads"])

    vals1 = [
        get("test03_G_noGap", "J2_Chr14_A2_to_Chr14_E"),
        get("test04_BG_noGap", "J2_Chr14_A2_to_Chr14_B"),
        get("test04_BG_noGap", "J3_Chr14_B_to_Chr14_E"),
    ]
    vals2 = [
        get("test01_direct_noGap", "J3_Chr14_E_to_Chr14_H1"),
        get("test03_G_noGap", "J3_Chr14_E_to_Chr14_G"),
        get("test03_G_noGap", "J4_Chr14_G_to_Chr14_H1"),
    ]

    c = pd.read_csv(cpath, sep="\t")
    if not {"gap_kb","spanning_reads"}.issubset(c.columns):
        raise SystemExit(f"{cpath} must contain gap_kb and spanning_reads; observed {list(c.columns)}")
    c["gap_kb"] = pd.to_numeric(c["gap_kb"], errors="raise").astype(int)
    c["spanning_reads"] = pd.to_numeric(c["spanning_reads"], errors="raise").astype(int)
    vals_gap = []
    for gap in (0,1,5,10):
        m = c[c["gap_kb"] == gap]
        if len(m) != 1:
            raise SystemExit(f"Expected one row for gap_kb={gap} in {cpath}; found {len(m)}")
        vals_gap.append(int(m.iloc[0]["spanning_reads"]))

    print("Panel b counts:", vals1, vals2)
    print("Panel c counts:", vals_gap)
    return vals1, vals2, vals_gap

def bar_panel(ax, title, labels, values, selected_indices):
    y = np.arange(len(labels))
    colors = ["#B8B8B8"] * len(labels)
    for i in selected_indices:
        colors[i] = "#333333"
    bars = ax.barh(y, values, color=colors, height=0.44)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 25)
    ax.set_xticks([0,5,10,15,20,25])
    ax.set_xlabel("Boundary-spanning primary HiFi alignments", fontsize=6.2)
    ax.set_title(title, fontweight="bold", pad=5)
    ax.spines[["top","right","left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=5.3, pad=3)
    for yi, (bar, v) in enumerate(zip(bars, values)):
        ax.text(max(v, 0) + 0.45, yi, str(v), va="center", ha="left", fontsize=6.3)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strict-spanning", required=True, help="SuppFig2B_strict_spanning_data.tsv")
    ap.add_argument("--gap-sweep", required=True, help="SuppFig2C_GH1_gap_sweep_data.tsv")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--stem", default="Supplementary_Figure2_code_base")
    args = ap.parse_args()

    set_pub_style(7.0)
    vals1, vals2, gapvals = load_counts(args.strict_spanning, args.gap_sweep)

    # ~180 mm wide.
    fig = plt.figure(figsize=(7.08, 6.25))

    # a
    axa = fig.add_axes([0.055, 0.715, 0.89, 0.255])
    axa.set_xlim(0, 1); axa.set_ylim(0, 1); axa.axis("off")
    axa.text(-0.030, 1.00, "a", fontweight="bold", fontsize=9, va="top")
    axa.text(0.53, 0.995, "Four alternative backbone topologies",
             ha="center", va="top", fontweight="bold")
    ys = [0.78, 0.56, 0.34, 0.12]
    for (name, items), y in zip(MODELS, ys):
        if name == "Model 3":
            # Highlight only topology area, not the left-side model label.
            axa.add_patch(FancyBboxPatch(
                (0.135, y - 0.075), 0.80, 0.15,
                boxstyle="round,pad=0.008,rounding_size=0.010",
                facecolor="#F4F8FF", edgecolor=COLORS["blue"], lw=0.9,
            ))
        axa.text(0.015, y, name, ha="left", va="center", fontweight="bold",
                 color=COLORS["blue"] if name == "Model 3" else "black")
        # Reserve a narrow right margin for the adopted-model badge.
        draw_boxes(axa, items, 0.145, y, available=0.73, gap=0.006, h=0.108, fs=4.65)
        if name == "Model 3":
            axa.text(0.900, y, "Adopted\nmodel", ha="center", va="center",
                     color=COLORS["blue"], fontweight="bold", fontsize=5.8, linespacing=0.90)

    # b
    axb1 = fig.add_axes([0.12, 0.445, 0.31, 0.205])
    axb2 = fig.add_axes([0.635, 0.445, 0.31, 0.205])
    axb1.text(-0.18, 1.12, "b", transform=axb1.transAxes, fontweight="bold", fontsize=9)
    labs1 = [
        "ptg000027l–\nptg000015l-derived",
        "ptg000027l–\nptg000010l",
        "ptg000010l–\nptg000015l-derived",
    ]
    labs2 = [
        "ptg000015l-derived–\nptg000003l",
        "ptg000015l-derived–\nptg000020l",
        "ptg000020l–\nptg000003l",
    ]
    bar_panel(axb1, "Testing ptg000010l", labs1, vals1, [0])
    bar_panel(axb2, "Testing ptg000020l", labs2, vals2, [1,2])

    # Larger, clearer schematics under b.
    axsc = fig.add_axes([0.075, 0.275, 0.87, 0.115])
    axsc.axis("off"); axsc.set_xlim(0,1); axsc.set_ylim(0,1)

    axsc.text(0.01,0.82,"Without ptg000010l",fontweight="bold",fontsize=6.2)
    draw_boxes(axsc,["ptg000027l","ptg000015l-derived"],0.01,0.53,available=0.37,gap=0.012,h=0.19,fs=5.3)
    axsc.text(0.01,0.20,"With ptg000010l",fontweight="bold",fontsize=6.2)
    draw_boxes(axsc,["ptg000027l","ptg000010l","ptg000015l-derived"],0.01,-0.01,available=0.37,gap=0.010,h=0.19,fs=5.1)

    axsc.text(0.55,0.82,"Without ptg000020l",fontweight="bold",fontsize=6.2)
    draw_boxes(axsc,["ptg000015l-derived","ptg000003l"],0.55,0.53,available=0.37,gap=0.012,h=0.19,fs=5.3)
    axsc.text(0.55,0.20,"With ptg000020l",fontweight="bold",fontsize=6.2)
    draw_boxes(axsc,["ptg000015l-derived","ptg000020l","ptg000003l"],0.55,-0.01,available=0.37,gap=0.010,h=0.19,fs=5.1)

    # c
    axc = fig.add_axes([0.095, 0.065, 0.43, 0.155])
    gaps = [0,1,5,10]
    bars = axc.bar(range(4), gapvals, color=["#333333","#AFAFAF","#BBBBBB","#C9C9C9"], width=0.55)
    axc.set_xticks(range(4), ["0","1","5","10"])
    axc.set_xlabel("Inserted N-spacer (kb)")
    axc.set_ylabel("Boundary-spanning\nprimary HiFi alignments")
    axc.spines[["top","right"]].set_visible(False)
    axc.text(-0.19,1.12,"c",transform=axc.transAxes,fontweight="bold",fontsize=9)
    for b,v in zip(bars,gapvals):
        axc.text(b.get_x()+b.get_width()/2, v+0.65, str(v), ha="center", fontsize=6.4)

    axcs = fig.add_axes([0.615, 0.055, 0.33, 0.165])
    axcs.axis("off"); axcs.set_xlim(0,1); axcs.set_ylim(0,1)
    axcs.text(0,0.82,"No N-spacer",fontweight="bold")
    draw_boxes(axcs,["ptg000020l","ptg000003l"],0,0.60,available=0.88,gap=0.04,h=0.21,fs=5.5)
    axcs.text(0,0.31,"With N-spacer",fontweight="bold")
    draw_boxes(axcs,["ptg000020l","N-spacer","ptg000003l"],0,0.08,available=0.88,gap=0.035,h=0.21,fs=5.2)

    save_triplet_fixed(fig, args.outdir, args.stem, dpi=450)
    plt.close(fig)

if __name__ == "__main__":
    main()
