#!/usr/bin/env python3
from pathlib import Path
import matplotlib as mpl
from matplotlib import font_manager

COLORS = {"green":"#009E73","blue":"#0072B2","vermillion":"#D55E00","orange":"#E69F00","purple":"#CC79A7","gray":"#666666"}

def _font_family():
    installed={f.name for f in font_manager.fontManager.ttflist}
    for x in ["Arial","Helvetica","Liberation Sans","DejaVu Sans"]:
        if x in installed: return x
    return "DejaVu Sans"

def set_pub_style(font_size=7.0):
    fam=_font_family()
    mpl.rcParams.update({"font.family":fam,"font.size":font_size,"axes.titlesize":font_size,"axes.labelsize":font_size,"xtick.labelsize":font_size-0.5,"ytick.labelsize":font_size-0.5,"legend.fontsize":font_size-0.5,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","axes.linewidth":0.6})
    return fam

def save_triplet_fixed(fig,outdir,stem,dpi=450):
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    for ext in ("pdf","svg"):
        p=out/f"{stem}.{ext}"; fig.savefig(p,facecolor="white",transparent=False,bbox_inches=None,pad_inches=0); print("Wrote:",p)
    p=out/f"{stem}.png"; fig.savefig(p,dpi=dpi,facecolor="white",transparent=False,bbox_inches=None,pad_inches=0); print("Wrote:",p)
