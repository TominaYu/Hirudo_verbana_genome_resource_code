#!/usr/bin/env python3
"""Targeted finalization of Supplementary Figure 4 from the final mitochondrial annotation TSV.

Changes after inspecting Supplementary_Figure4_v33_refresh.pdf:
- larger ~180-mm-wide canvas and more vertical room for 37 direct labels;
- two-stage elbow leader lines rather than direct long fan lines;
- label positions are spread monotonically on each side to avoid crossings;
- derived dense lower-right labels get more separation;
- retain direct labels, trnL1(TAG)/trnL2(TAA), 5.30-kb unannotated interval,
  and the 13 PCG / 22 tRNA / 2 rRNA summary.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Rectangle
from common import COLORS, set_pub_style, save_triplet_fixed

GENOME_LEN = 19644
TYPE_COLORS = {"PCG": "#4D4D4D", "tRNA": COLORS["orange"], "rRNA": COLORS["blue"]}

def load_annotation(path: Path):
    raw = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str,
                      keep_default_na=False, engine="python")
    raw = raw.dropna(how="all")
    if raw.shape[1] < 6:
        raise SystemExit(f"Expected >=6 tab-separated columns in {path}; got {raw.shape[1]}")
    s = pd.to_numeric(raw.iloc[:,0].astype(str).str.replace(",","",regex=False), errors="coerce")
    e = pd.to_numeric(raw.iloc[:,1].astype(str).str.replace(",","",regex=False), errors="coerce")
    keep = s.notna() & e.notna()
    out = pd.DataFrame({
        "start": s.loc[keep].astype(int).to_numpy(),
        "end": e.loc[keep].astype(int).to_numpy(),
        "type": raw.loc[keep, raw.columns[3]].astype(str).str.strip().to_numpy(),
        "label": raw.loc[keep, raw.columns[5]].astype(str).str.strip().to_numpy(),
    })
    t = out["type"].str.lower()
    out = out[t.isin(["pcg","trna","rrna"])].copy()
    out = out[(out.start.between(1,GENOME_LEN)) & (out.end.between(1,GENOME_LEN))].reset_index(drop=True)
    counts = out["type"].str.lower().value_counts().to_dict()
    print("Parsed mtDNA features:", len(out),
          f"(PCG={counts.get('pcg',0)}, tRNA={counts.get('trna',0)}, rRNA={counts.get('rrna',0)})")
    if (counts.get("pcg",0),counts.get("trna",0),counts.get("rrna",0)) != (13,22,2):
        raise SystemExit("Unexpected mitochondrial feature complement; expected 13 PCGs, 22 tRNAs and 2 rRNAs.")
    return out

def theta_for(pos):
    return 2*np.pi*(float(pos) % GENOME_LEN)/GENOME_LEN

def spread_side(items, low=-1.10, high=1.10, minsep=0.082):
    """Preserve vertical order while enforcing a minimum label separation."""
    items = sorted(items, key=lambda z: z[0])
    if not items:
        return {}
    ys = [max(low, min(high, y)) for y,_ in items]
    for i in range(1,len(ys)):
        ys[i] = max(ys[i], ys[i-1] + minsep)
    if ys[-1] > high:
        shift = ys[-1] - high
        ys = [y-shift for y in ys]
    for i in range(len(ys)-2,-1,-1):
        ys[i] = min(ys[i], ys[i+1] - minsep)
    if ys[0] < low:
        shift = low - ys[0]
        ys = [y+shift for y in ys]
    return {idx: y for y,(_,idx) in zip(ys, items)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Curated 19,644-bp mitochondrial annotation TSV")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--stem", default="Supplementary_Figure4_code_base")
    args = ap.parse_args()

    set_pub_style(7.0)
    src = Path(args.input)
    print("Using:", src)

    df = load_annotation(src)
    def cls(x):
        s = x.lower()
        if "trna" in s: return "tRNA"
        if "rrna" in s: return "rRNA"
        return "PCG"
    df["class"] = df["type"].map(cls)
    df["mid"] = (df.start + df.end) / 2

    # Wider/taller than v33 to make dense direct labels publication-readable.
    fig, ax = plt.subplots(figsize=(7.08, 5.65))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.62,1.62)
    ax.set_ylim(-1.38,1.38)

    for r in [0.72,0.82,0.95,1.02]:
        ax.add_patch(plt.Circle((0,0), r, fill=False, edgecolor="#D0D0D0", lw=0.42))

    for _, row in df.iterrows():
        a1 = 90 - 360*row.end/GENOME_LEN
        a2 = 90 - 360*row.start/GENOME_LEN
        if row["class"] == "tRNA":
            r0,r1 = 0.77,0.82
        else:
            r0,r1 = 0.86,1.00
        ax.add_patch(Wedge(
            (0,0), r1, a1, a2, width=r1-r0,
            facecolor=TYPE_COLORS[row["class"]],
            edgecolor="white", lw=0.35,
        ))

    def row_for(lbl):
        m = df[df.label.astype(str).str.lower() == lbl.lower()]
        return None if m.empty else m.iloc[0]

    trnr, trnh = row_for("trnR"), row_for("trnH")
    if trnr is not None and trnh is not None:
        start = int(trnr.end)+1
        end = int(trnh.start)-1
        arcs = [(start,GENOME_LEN),(0,end)] if start > end else [(start,end)]
        for s,e in arcs:
            a1 = 90 - 360*e/GENOME_LEN
            a2 = 90 - 360*s/GENOME_LEN
            ax.add_patch(Wedge((0,0),1.115,a1,a2,width=0.052,
                               facecolor="#D8CEE9",edgecolor="none",alpha=0.95))
        ax.text(0,1.235,"Unannotated interval (5.30 kb)",
                ha="center",va="center",fontweight="bold",fontsize=7.4)

    for kb in [0,5,10,15]:
        th = theta_for(kb*1000)
        x1,y1 = 0.67*np.sin(th),0.67*np.cos(th)
        x2,y2 = 0.71*np.sin(th),0.71*np.cos(th)
        ax.plot([x1,x2],[y1,y2],color="#666666",lw=0.48)
        ax.text(0.61*np.sin(th),0.61*np.cos(th),str(kb),ha="center",va="center",fontsize=6.5)

    mids = [theta_for(v) for v in df.mid]
    left, right = [], []
    for i, th in enumerate(mids):
        x, y = np.sin(th), np.cos(th)
        # Near the bottom, split by x sign; this preserves genomic order while balancing sides.
        (right if x >= 0 else left).append((y, i))
    ly = spread_side(left, minsep=0.086)
    ry = spread_side(right, minsep=0.086)

    # Ordered elbow leaders: short radial segment -> lane -> label.
    for i, row in df.iterrows():
        th = mids[i]
        x, y = np.sin(th), np.cos(th)
        side = 1 if x >= 0 else -1
        target_y = (ry if side == 1 else ly)[i]

        xa, ya = 1.005*x, 1.005*y
        xr, yr = 1.075*x, 1.075*y

        # Use an independent elbow for every feature.  The previous shared
        # near-vertical lane produced a "comb" of overlapping leaders in the
        # dense lower-right part of the map.  Because target_y preserves the
        # anchor order on each side, these diagonal elbows remain ordered.
        elbow_x = side * 1.27
        label_x = side * 1.48

        ax.plot([xa,xr],[ya,yr],color="#555555",lw=0.40)
        ax.plot([xr,elbow_x],[yr,target_y],color="#666666",lw=0.38)
        ax.plot([elbow_x,label_x],[target_y,target_y],color="#555555",lw=0.40)

        ax.text(label_x + 0.025*side, target_y, str(row.label),
                ha="left" if side==1 else "right",va="center",fontsize=6.2)

    ax.text(0,0.15,r"$\it{Hirudo\ verbana}$",ha="center",va="center",fontsize=8.2,fontweight="bold")
    ax.text(0,0.03,"mitochondrial genome",ha="center",fontsize=7.0)
    ax.text(0,-0.10,"19,644 bp",ha="center",fontsize=8.0)
    ax.text(0,-0.235,"13 protein-coding genes · 22 tRNAs · 2 rRNAs",
            ha="center",fontsize=6.3)

    yleg = -1.25
    entries = [
        (-0.72,"Protein-coding genes","PCG"),
        (0.08,"tRNA","tRNA"),
        (0.50,"rRNA","rRNA"),
    ]
    for x,label,cl in entries:
        ax.add_patch(Rectangle((x,yleg-0.024),0.075,0.048,
                               facecolor=TYPE_COLORS[cl],edgecolor="none"))
        ax.text(x+0.095,yleg,label,va="center",ha="left",fontsize=6.4)

    save_triplet_fixed(fig, args.outdir, args.stem, dpi=450)
    plt.close(fig)

if __name__ == "__main__":
    main()
