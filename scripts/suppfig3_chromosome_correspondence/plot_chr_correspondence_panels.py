#!/usr/bin/env python3
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

WORK = Path("09_suppfig_chr_correspondence")
coords_file = WORK / "02_tables/our_vs_wcHirVerb1.1to1.coords.tsv"
best_file = WORK / "02_tables/chr_correspondence_best_pairs.tsv"
out_png = WORK / "03_figures/SuppFig_chr_correspondence_one_to_one.png"
out_pdf = WORK / "03_figures/SuppFig_chr_correspondence_one_to_one.pdf"

# read best mapping
best_map = []
with best_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        best_map.append(row)

best_lookup = {r["eu_chr"]: r["best_our_chr"] for r in best_map}

# load coords
blocks = defaultdict(list)
eu_len = {}
our_len = {}

with coords_file.open() as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if not row or row[0].startswith("/"):
            continue
        if len(row) < 13:
            continue
        s1, e1, s2, e2 = map(int, row[:4])
        len1 = int(row[4])
        pid = float(row[6])
        ref_len = int(row[7])
        qry_len = int(row[8])
        ref_id = row[11]
        qry_id = row[12]

        if ref_id not in best_lookup:
            continue
        if qry_id != best_lookup[ref_id]:
            continue

        eu_len[ref_id] = ref_len
        our_len[qry_id] = qry_len

        orient = "forward" if s2 <= e2 else "reverse"
        blocks[ref_id].append({
            "s1": s1, "e1": e1,
            "s2": s2, "e2": e2,
            "len1": len1,
            "pid": pid,
            "qry_id": qry_id,
            "orient": orient
        })

# sort chromosomes by EU numeric order if possible
def chr_key(x):
    import re
    m = re.search(r'chromosome:\s*(\d+)', x)
    if m:
        return int(m.group(1))
    m = re.search(r'chr(?:omosome)?[_ ]?(\d+)', x, re.IGNORECASE)
    if m:
        return int(m.group(1))
    m = re.search(r'(\d+)$', x)
    if m:
        return int(m.group(1))
    return x

eu_order = sorted(blocks.keys(), key=chr_key)

n = len(eu_order)
ncol = 2
nrow = math.ceil(n / ncol)

fig, axes = plt.subplots(nrow, ncol, figsize=(14, 2.3 * nrow))
axes = axes.flatten()

for ax in axes[n:]:
    ax.axis("off")

for ax, eu_chr in zip(axes, eu_order):
    our_chr = best_lookup[eu_chr]
    elen = eu_len[eu_chr]
    qlen = our_len[our_chr]

    ax.set_xlim(0, max(elen, qlen))
    ax.set_ylim(0, 1)
    ax.axis("off")

    # bar positions
    y_ref = 0.72
    y_qry = 0.22
    h = 0.08

    ax.add_patch(Rectangle((0, y_ref), elen, h, color="#D0D0D0", lw=0))
    ax.add_patch(Rectangle((0, y_qry), qlen, h, color="#D0D0D0", lw=0))

    # blocks
    for b in sorted(blocks[eu_chr], key=lambda x: x["s1"]):
        x1a, x1b = b["s1"], b["e1"]
        x2a, x2b = sorted([b["s2"], b["e2"]])
        color = "#4C78A8" if b["orient"] == "forward" else "#F58518"
        poly = Polygon(
            [(x1a, y_ref), (x1b, y_ref), (x2b, y_qry + h), (x2a, y_qry + h)],
            closed=True, facecolor=color, edgecolor="none", alpha=0.75
        )
        ax.add_patch(poly)

    # labels
    ax.text(0, 0.95, eu_chr, fontsize=9, ha="left", va="top")
    ax.text(0, 0.05, our_chr, fontsize=9, ha="left", va="bottom")

    ax.text(elen, 0.95, f"{elen/1e6:.2f} Mb", fontsize=8, ha="right", va="top")
    ax.text(qlen, 0.05, f"{qlen/1e6:.2f} Mb", fontsize=8, ha="right", va="bottom")

fig.suptitle(
    "One-to-one NUCmer alignment blocks between wcHirVerb1 chromosomes and our primary chromosomes",
    fontsize=14, y=0.995
)

plt.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(out_png, dpi=400)
fig.savefig(out_pdf)
print(f"Wrote: {out_png}")
print(f"Wrote: {out_pdf}")
