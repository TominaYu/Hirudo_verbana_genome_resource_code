#!/usr/bin/env python3
"""Render code-generated bases for Supplementary Figure 3."""
from __future__ import annotations
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from common import COLORS, set_pub_style


def chr_num(value):
    m=re.search(r"(\d+)$",str(value))
    if not m:
        raise ValueError(f"Could not parse chromosome number: {value}")
    return int(m.group(1))


def save(fig,out,stem,dpi=450):
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    for ext in ("pdf","svg"):
        p=out/f"{stem}.{ext}"
        fig.savefig(p,facecolor="white",bbox_inches="tight",pad_inches=0.02)
        print("Wrote:",p)
    p=out/f"{stem}.png"
    fig.savefig(p,dpi=dpi,facecolor="white",bbox_inches="tight",pad_inches=0.02)
    print("Wrote:",p)


def as_bool(value):
    if isinstance(value,(bool,np.bool_)):
        return bool(value)
    x=str(value).strip().lower()
    if x in {"true","1","yes"}:
        return True
    if x in {"false","0","no"}:
        return False
    raise ValueError(f"Could not parse boolean value: {value}")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pair-stats",required=True)
    ap.add_argument("--reciprocal-pairs",required=True)
    ap.add_argument("--visual-blocks",required=True)
    ap.add_argument("--outdir",default=".")
    args=ap.parse_args()
    set_pub_style(7.0)

    pairs=pd.read_csv(args.pair_stats,sep="\t")
    recip=pd.read_csv(args.reciprocal_pairs,sep="\t")
    recip=recip.assign(_order=recip["our_chr"].map(chr_num)).sort_values("_order").drop(columns="_order")
    blocks=pd.read_csv(args.visual_blocks,sep="\t")

    required_pair={"eu_chr","our_chr","aligned_bp"}
    required_recip={"eu_chr","our_chr","eu_length","our_length","dominant_orientation"}
    required_blocks={"eu_chr","our_chr","eu_start","eu_end","our_start","our_end","aligned_bp","same_orientation"}
    for name,df,req in [
        ("pair statistics",pairs,required_pair),
        ("reciprocal pairs",recip,required_recip),
        ("visualized blocks",blocks,required_blocks),
    ]:
        missing=req-set(df.columns)
        if missing:
            raise SystemExit(f"Missing columns in {name}: {sorted(missing)}")

    if len(recip)!=14:
        raise SystemExit(f"Expected 14 reciprocal-best pairs, found {len(recip)}")

    eu_order=sorted(pairs["eu_chr"].astype(str).unique(),key=chr_num)
    our_order=[f"Chr{x:02d}" for x in range(1,15)]
    mat=np.zeros((len(eu_order),len(our_order)),float)
    eu_i={x:i for i,x in enumerate(eu_order)}
    our_i={x:i for i,x in enumerate(our_order)}
    for _,r in pairs.iterrows():
        mat[eu_i[str(r.eu_chr)],our_i[str(r.our_chr)]]=float(r.aligned_bp)/1e6

    fig,ax=plt.subplots(figsize=(5.9,5.2))
    im=ax.imshow(mat,cmap="Blues",aspect="equal")
    ax.set_xticks(range(14),[x.replace("Chr","") for x in our_order])
    ax.set_yticks(range(14),[f"{chr_num(x):02d}" for x in eu_order])
    ax.set_xlabel("Hverb_1.0 chromosome")
    ax.set_ylabel("wcHirVerb1 chromosome")
    cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04)
    cb.set_label("Summed one-to-one aligned sequence (Mb)")
    rb={(str(r.eu_chr),str(r.our_chr)) for _,r in recip.iterrows()}
    for i,eu in enumerate(eu_order):
        for j,our in enumerate(our_order):
            if (eu,our) in rb:
                ax.add_patch(Rectangle((j-0.5,i-0.5),1,1,fill=False,edgecolor="black",linewidth=0.8))
    fig.tight_layout()
    save(fig,args.outdir,"SuppFig3_A_correspondence_matrix")
    plt.close(fig)

    fig,axes=plt.subplots(7,2,figsize=(7.2,10.2))
    axes=axes.ravel()
    for ax,(_,p) in zip(axes,recip.iterrows()):
        eu=str(p.eu_chr)
        our=str(p.our_chr)
        eu_len=int(p.eu_length)
        our_len=int(p.our_length)
        flip=(str(p.dominant_orientation).lower()=="reverse")
        sub=blocks[
            (blocks.eu_chr.astype(str)==eu)
            & (blocks.our_chr.astype(str)==our)
        ].copy()
        maxlen=max(eu_len,our_len)
        ax.set_xlim(0,maxlen)
        ax.set_ylim(0,1)
        ax.axis("off")
        y1,y2,h=0.72,0.18,0.075
        ax.add_patch(Rectangle((0,y1),eu_len,h,facecolor="#C8C8C8",edgecolor="none"))
        ax.add_patch(Rectangle((0,y2),our_len,h,facecolor="#C8C8C8",edgecolor="none"))
        sub=sub.sort_values("aligned_bp",ascending=False)
        for _,b in sub.iterrows():
            q1,q2=int(b.our_start),int(b.our_end)
            if flip:
                dq1=our_len-q2+1
                dq2=our_len-q1+1
            else:
                dq1,dq2=q1,q2
            same=as_bool(b.same_orientation)
            collinear=(same != flip)
            color=COLORS["blue"] if collinear else COLORS["orange"]
            pts=[
                (int(b.eu_start),y1),
                (int(b.eu_end),y1),
                (dq2,y2+h),
                (dq1,y2+h),
            ]
            ax.add_patch(
                Polygon(
                    pts,
                    closed=True,
                    facecolor=color,
                    edgecolor="none",
                    alpha=0.58,
                )
            )
        ax.text(0,0.96,f"wc {chr_num(eu):02d}",ha="left",va="top",fontsize=6.2)
        ax.text(0,0.04,our.replace("Chr",""),ha="left",va="bottom",fontsize=6.2)
        ax.text(maxlen,0.96,f"{eu_len/1e6:.2f} Mb",ha="right",va="top",fontsize=5.5)
        ax.text(maxlen,0.04,f"{our_len/1e6:.2f} Mb",ha="right",va="bottom",fontsize=5.5)
    fig.text(
        0.5,0.995,
        "wc = wcHirVerb1; Hverb_1.0 chromosome numbers are shown below each pair",
        ha="center",va="top",fontsize=6.5,
    )
    fig.tight_layout(rect=(0,0,1,0.985),h_pad=0.35,w_pad=0.35)
    save(fig,args.outdir,"SuppFig3_B_reciprocal_pairs")
    plt.close(fig)


if __name__=="__main__":
    main()
