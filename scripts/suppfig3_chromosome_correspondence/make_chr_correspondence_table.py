#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path

coords = Path("09_suppfig_chr_correspondence/02_tables/our_vs_wcHirVerb1.1to1.coords.tsv")
out1 = Path("09_suppfig_chr_correspondence/02_tables/chr_pair_alignment_summary.tsv")
out2 = Path("09_suppfig_chr_correspondence/02_tables/chr_correspondence_best_pairs.tsv")

rows = []
with coords.open() as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if not row or row[0].startswith("/"):
            continue
        # show-coords -T columns
        # [S1] [E1] [S2] [E2] [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]
        if len(row) < 13:
            continue
        s1, e1, s2, e2 = map(int, row[:4])
        len1, len2 = map(int, row[4:6])
        pid = float(row[6])
        len_r = int(row[7])
        len_q = int(row[8])
        cov_r = float(row[9])
        cov_q = float(row[10])
        ref_id = row[11]
        qry_id = row[12]
        rows.append({
            "ref_id": ref_id,
            "qry_id": qry_id,
            "len1": len1,
            "len2": len2,
            "pid": pid,
            "ref_len": len_r,
            "qry_len": len_q,
        })

pair_bp = defaultdict(int)
pair_blocks = defaultdict(int)
ref_lengths = {}
qry_lengths = {}

for r in rows:
    pair = (r["ref_id"], r["qry_id"])
    pair_bp[pair] += r["len1"]
    pair_blocks[pair] += 1
    ref_lengths[r["ref_id"]] = r["ref_len"]
    qry_lengths[r["qry_id"]] = r["qry_len"]

with out1.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow([
        "eu_chr", "our_chr", "aligned_bp", "block_count",
        "eu_chr_length", "our_chr_length",
        "eu_coverage_pct", "our_coverage_pct"
    ])
    for (eu_chr, our_chr), bp in sorted(pair_bp.items()):
        eu_len = ref_lengths[eu_chr]
        our_len = qry_lengths[our_chr]
        w.writerow([
            eu_chr,
            our_chr,
            bp,
            pair_blocks[(eu_chr, our_chr)],
            eu_len,
            our_len,
            100.0 * bp / eu_len if eu_len else 0.0,
            100.0 * bp / our_len if our_len else 0.0,
        ])

# best pair per EU chromosome
best = {}
for (eu_chr, our_chr), bp in pair_bp.items():
    if eu_chr not in best or bp > best[eu_chr][1]:
        best[eu_chr] = (our_chr, bp)

with out2.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow([
        "eu_chr", "best_our_chr", "aligned_bp",
        "eu_chr_length", "our_chr_length",
        "eu_coverage_pct", "our_coverage_pct"
    ])
    for eu_chr in sorted(best):
        our_chr, bp = best[eu_chr]
        eu_len = ref_lengths[eu_chr]
        our_len = qry_lengths[our_chr]
        w.writerow([
            eu_chr,
            our_chr,
            bp,
            eu_len,
            our_len,
            100.0 * bp / eu_len if eu_len else 0.0,
            100.0 * bp / our_len if our_len else 0.0,
        ])

print(f"Wrote: {out1}")
print(f"Wrote: {out2}")
