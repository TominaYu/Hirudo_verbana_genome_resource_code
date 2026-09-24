#!/usr/bin/env python3
"""Render the code-generated base for main Figure 2.

Reads the audited 100-kb source table. Scientific display ranges match the manuscript figure caption.
Key changes after visual review of Figure2_v33_refresh.pdf:
- output physical width fixed at ~180 mm;
- chromosome names moved inside the pale outer ring;
- local 0/10-Mb coordinate labels moved outside the ring and made smaller,
  avoiding collisions such as Chr01/10 and Chr02/10;
- chromosome labels are kept upright around the circle;
- compact right-side legend; same a-g order as the manuscript;
- color-blind-aware palette, with position providing non-color redundancy;
- line widths kept >= ~0.35 pt.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Rectangle
from common import COLORS, set_pub_style, save_triplet_fixed

DISPLAY = {
    "gc_skew": (-0.25, 0.25),
    "gc_fraction": (0.25, 0.60),
    "gene_density_overlapping_gene_spans": (0, 80),
    "LINE_fraction": (0, 0.35),
    "LTR_fraction": (0, 0.20),
    "interspersed_TE_fraction": (0, 0.75),
}
TRACKS = [
    ("gc_skew", "GC skew", COLORS["green"], "−0.25 to 0.25", "line"),
    ("gc_fraction", "GC fraction", COLORS["gray"], "0.25 to 0.60", "fill"),
    ("gene_density_overlapping_gene_spans", "Gene density", COLORS["blue"], "0 to 80 genes per 100 kb", "fill"),
    ("LINE_fraction", "LINE fraction", COLORS["vermillion"], "0 to 0.35", "fill"),
    ("LTR_fraction", "LTR fraction", COLORS["orange"], "0 to 0.20", "fill"),
    ("interspersed_TE_fraction", "Interspersed TE fraction", COLORS["purple"], "0 to 0.75", "fill"),
]

def norm_clip(v, lo, hi):
    a = np.asarray(v, float)
    return np.clip((a - lo) / (hi - lo), 0, 1)

def upright_rotation(theta):
    """Tangential-ish text rotation that never leaves text upside down."""
    rot = -np.degrees(theta)
    while rot <= -180:
        rot += 360
    while rot > 180:
        rot -= 360
    if rot < -90:
        rot += 180
    elif rot > 90:
        rot -= 180
    return rot

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Figure 2 100-kb source-data TSV")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--stem", default="Figure2_code_base")
    args = ap.parse_args()

    set_pub_style(7.0)
    src = Path(args.input)
    print("Using:", src)
    df = pd.read_csv(src, sep="\t")
    required = {"chromosome", "window_start_0based", "window_end", *DISPLAY.keys()}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {sorted(missing)}")

    chroms = list(dict.fromkeys(df["chromosome"].astype(str)))
    lengths = df.groupby("chromosome", sort=False)["window_end"].max().reindex(chroms).astype(int)
    offsets = lengths.cumsum().shift(fill_value=0)
    total = int(lengths.sum())

    mid = (df["window_start_0based"].to_numpy(float) + df["window_end"].to_numpy(float)) / 2
    global_mid = mid + df["chromosome"].map(offsets.to_dict()).to_numpy(float)
    theta = 2 * np.pi * global_mid / total

    # Exact ~180-mm width. Keep every element inside the figure to avoid tight-bbox expansion.
    fig = plt.figure(figsize=(7.08, 4.85))
    ax = fig.add_axes([0.015, 0.035, 0.705, 0.93], projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_axis_off()
    ax.set_ylim(0, 1.19)

    outer_r0, outer_r1 = 0.985, 1.075
    for i, chrom in enumerate(chroms):
        start = float(offsets[chrom])
        end = start + float(lengths[chrom])
        t1 = 360 * start / total
        t2 = 360 * end / total
        color = "#B9C7D2" if i % 2 == 0 else "#E3E8EC"
        ax.add_patch(Wedge(
            (0, 0), outer_r1, 90 - t2, 90 - t1,
            width=outer_r1 - outer_r0,
            facecolor=color, edgecolor="white", linewidth=0.7,
            transform=ax.transData._b,
        ))

        # Chromosome ID inside the neutral ring, separated from coordinate labels.
        tm = 2 * np.pi * ((start + end) / 2) / total
        ax.text(
            tm, (outer_r0 + outer_r1) / 2, chrom,
            ha="center", va="center", fontsize=6.2,
            rotation=upright_rotation(tm), rotation_mode="anchor",
        )

        # Local coordinate labels outside the ring.
        # Put 0- and 10-Mb labels on two radial levels.  This prevents
        # labels on opposite sides of a chromosome boundary from visually
        # merging into strings such as "100".
        for bp, lab, label_r in [
            (0, "0", outer_r1 + 0.032),
            (10_000_000, "10", outer_r1 + 0.066),
        ]:
            if bp <= lengths[chrom]:
                tt = 2 * np.pi * (start + bp) / total
                tick_end = outer_r1 + (0.013 if bp == 0 else 0.021)
                ax.plot([tt, tt], [outer_r1, tick_end], color="#333333", lw=0.45)
                ax.text(
                    tt, label_r, lab, fontsize=5.1,
                    ha="center", va="center",
                    rotation=upright_rotation(tt), rotation_mode="anchor",
                )

    centers = [0.915, 0.805, 0.695, 0.585, 0.475, 0.365]
    half = 0.044
    for (col, label, color, rng, mode), c in zip(TRACKS, centers):
        lo, hi = DISPLAY[col]
        vals = norm_clip(df[col].to_numpy(float), lo, hi)
        base = c - half
        rr = base + vals * (2 * half)
        grid_theta = np.linspace(0, 2 * np.pi, 721)
        ax.plot(grid_theta, np.full(721, base), color="#D5D5D5", lw=0.35)
        ax.plot(grid_theta, np.full(721, c + half), color="#ECECEC", lw=0.35)
        if mode == "line":
            ax.plot(theta, rr, color=color, lw=0.48)
        else:
            ax.fill_between(theta, base, rr, color=color, alpha=0.96, linewidth=0)

    for chrom in chroms:
        tt = 2 * np.pi * float(offsets[chrom]) / total
        ax.plot([tt, tt], [0.29, outer_r0], color="#C7C7C7", lw=0.35, zorder=0)

    # Compact legend fully inside the 180-mm canvas.
    lax = fig.add_axes([0.725, 0.10, 0.265, 0.80])
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    lax.text(0.00, 1.00, "Tracks (outer to inner)", fontsize=7.2, fontweight="bold", va="top")

    rows = [("a", "Chromosome-level records", "sequence coordinate (Mb)", None, "chrom")]
    for letter, tr in zip("bcdefg", TRACKS):
        col, label, color, rng, mode = tr
        rows.append((letter, label, rng, color, mode))

    y = 0.89
    for letter, label, detail, color, mode in rows:
        lax.text(0.00, y, letter, fontweight="bold", va="center")
        gx0, gx1 = 0.10, 0.26
        if mode == "chrom":
            lax.add_patch(Rectangle((gx0, y - 0.020), gx1 - gx0, 0.040,
                                    facecolor="#C5CFD8", edgecolor="#444", lw=0.45))
        elif mode == "line":
            x = np.linspace(gx0, gx1, 35)
            yy = y + 0.010 * np.sin(np.linspace(0, 4 * np.pi, 35))
            lax.plot(x, yy, color=color, lw=1.0, clip_on=False)
        else:
            lax.add_patch(Rectangle((gx0, y - 0.020), gx1 - gx0, 0.040,
                                    facecolor=color, edgecolor="none"))
        lax.text(0.32, y + 0.010, label, va="center", fontsize=6.4)
        lax.text(0.32, y - 0.026, detail, va="center", fontsize=5.6)
        y -= 0.125

    save_triplet_fixed(fig, args.outdir, args.stem, dpi=450)
    plt.close(fig)

if __name__ == "__main__":
    main()
