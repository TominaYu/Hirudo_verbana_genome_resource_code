#!/usr/bin/env python3
"""Build Supplementary Figure 3 chromosome-correspondence tables.

The public Figshare package contains:
  - one-to-one NUCmer coordinates; and
  - Supplementary Table S4, which supplies the wcHirVerb1 chromosome
    accession/length mapping.

Reciprocal-best pairs are selected independently in both directions by
summed one-to-one aligned sequence. Coverage is the union of merged
one-to-one alignment intervals, not summed aligned length divided by
record length.
"""
from __future__ import annotations
import argparse, csv, re
from collections import defaultdict
from pathlib import Path

def chr_number(value):
    m = re.search(r"(\d+)(?!.*\d)", str(value))
    if not m:
        raise ValueError(f"Could not parse chromosome number: {value}")
    return int(m.group(1))

def union_length(intervals):
    merged=[]
    for a,b in sorted((min(a,b),max(a,b)) for a,b in intervals):
        if not merged or a > merged[-1][1] + 1:
            merged.append([a,b])
        else:
            merged[-1][1]=max(merged[-1][1],b)
    return sum(b-a+1 for a,b in merged)

def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No rows for {path}")
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def meta_from_s4(path):
    with open(path,encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    required = {
        "wcHirVerb1 record",
        "wcHirVerb1 accession",
        "wcHirVerb1 length (bp)",
    }
    if not rows:
        raise SystemExit("S4 is empty")
    missing=required-set(rows[0])
    if missing:
        raise SystemExit(f"S4 missing required columns: {sorted(missing)}")
    meta={}
    for r in rows:
        n=chr_number(r["wcHirVerb1 record"])
        key=f"EU_chr{n:02d}"
        meta[key]={
            "number":n,
            "accession":r["wcHirVerb1 accession"].strip(),
            "length":int(r["wcHirVerb1 length (bp)"].replace(",","")),
        }
    if sorted(x["number"] for x in meta.values()) != list(range(1,15)):
        raise SystemExit(f"Expected wcHirVerb1 chromosomes 1-14 from S4; got {sorted(x['number'] for x in meta.values())}")
    return meta

def meta_from_accessions(path):
    meta={}
    with open(path,encoding="utf-8",newline="") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            meta[r["eu_chr"]]={
                "number":int(r["eu_chr_number"]),
                "accession":r["accession"],
                "length":int(r["length_bp"]),
            }
    return meta

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--coords",required=True)
    src=ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--s4",help="Supplementary Table S4 from the Figshare package")
    src.add_argument("--accessions",help="Legacy accession-map TSV (QA/backward compatibility)")
    ap.add_argument("--outdir",required=True)
    ap.add_argument("--visual-min-block",type=int,default=20000)
    ap.add_argument("--visual-min-identity",type=float,default=95.0)
    args=ap.parse_args()

    meta=meta_from_s4(args.s4) if args.s4 else meta_from_accessions(args.accessions)
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    blocks=[]
    with open(args.coords,encoding="utf-8") as f:
        for lineno,raw in enumerate(f,1):
            if not raw.strip() or raw.startswith("#"):
                continue
            x=raw.rstrip("\n").split()
            if len(x)<13:
                raise ValueError(f"Expected >=13 fields at line {lineno}, found {len(x)}")
            s1,e1,s2,e2=map(int,x[:4])
            len1,len2=map(int,x[4:6])
            ident=float(x[6]); lenr,lenq=map(int,x[7:9])
            ref,qry=x[11],x[12]
            if ref not in meta:
                raise ValueError(f"Unexpected wcHirVerb1 record in coords: {ref}")
            if not re.fullmatch(r"Chr(?:0[1-9]|1[0-4])",qry):
                raise ValueError(f"Unexpected Hverb_1.0 record in coords: {qry}")
            if lenr != meta[ref]["length"]:
                raise ValueError(f"{ref} length mismatch: coords={lenr}, S4/map={meta[ref]['length']}")
            blocks.append({
                "eu_chr":ref,"our_chr":qry,
                "eu_start":min(s1,e1),"eu_end":max(s1,e1),
                "our_start":min(s2,e2),"our_end":max(s2,e2),
                "eu_length":lenr,"our_length":lenq,
                "aligned_bp":min(len1,len2),"identity_pct":ident,
                "same_orientation":((e1-s1)*(e2-s2)>=0),
            })

    grouped=defaultdict(list)
    for b in blocks:
        grouped[(b["eu_chr"],b["our_chr"])].append(b)

    pair=[]
    for (eu,our),rows in grouped.items():
        aligned=sum(r["aligned_bp"] for r in rows)
        same=sum(r["aligned_bp"] for r in rows if r["same_orientation"])
        eu_cov=union_length((r["eu_start"],r["eu_end"]) for r in rows)
        our_cov=union_length((r["our_start"],r["our_end"]) for r in rows)
        eu_len=rows[0]["eu_length"]; our_len=rows[0]["our_length"]
        pair.append({
            "eu_chr":eu,
            "eu_chr_number":meta[eu]["number"],
            "eu_accession":meta[eu]["accession"],
            "our_chr":our,
            "aligned_bp":aligned,
            "block_count":len(rows),
            "same_orientation_bp":same,
            "reverse_orientation_bp":aligned-same,
            "dominant_orientation":"same" if same>=aligned-same else "reverse",
            "eu_length":eu_len,
            "our_length":our_len,
            "eu_union_covered_bp":eu_cov,
            "our_union_covered_bp":our_cov,
            "eu_union_coverage_pct":100*eu_cov/eu_len,
            "our_union_coverage_pct":100*our_cov/our_len,
        })
    pair.sort(key=lambda r:(r["eu_chr_number"],chr_number(r["our_chr"])))

    best_our={}
    for eu in meta:
        cand=[r for r in pair if r["eu_chr"]==eu]
        if not cand:
            raise SystemExit(f"No alignment found for {eu}")
        best_our[eu]=max(cand,key=lambda r:r["aligned_bp"])

    best_eu={}
    for our in [f"Chr{x:02d}" for x in range(1,15)]:
        cand=[r for r in pair if r["our_chr"]==our]
        if not cand:
            raise SystemExit(f"No alignment found for {our}")
        best_eu[our]=max(cand,key=lambda r:r["aligned_bp"])

    reciprocal=[r for eu,r in best_our.items() if best_eu[r["our_chr"]]["eu_chr"]==eu]
    reciprocal.sort(key=lambda r:chr_number(r["our_chr"]))
    if len(reciprocal)!=14:
        raise SystemExit(f"Expected 14 reciprocal-best pairs, found {len(reciprocal)}")

    visual=[b for b in blocks
            if b["aligned_bp"]>=args.visual_min_block
            and b["identity_pct"]>=args.visual_min_identity]

    write_tsv(out/"all_chromosome_pair_statistics.tsv",pair)
    write_tsv(out/"reciprocal_best_chromosome_pairs.tsv",reciprocal)
    write_tsv(out/"visualized_one_to_one_blocks.tsv",visual)
    print("All one-to-one blocks:",len(blocks))
    print("Pair rows:",len(pair))
    print("Reciprocal-best pairs:",len(reciprocal))
    print("Visualized blocks:",len(visual))
    print("Coverage definition: union of merged one-to-one alignment intervals")

if __name__=="__main__":
    main()
