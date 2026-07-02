#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess
from collections import Counter
from pathlib import Path

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_OPS = {"M", "D", "N", "=", "X"}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--focus-start", type=int, required=True)
    p.add_argument("--focus-end", type=int, required=True)
    p.add_argument("--out-prefix", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    proc = subprocess.run(
        ["samtools", "view", args.bam],
        check=True,
        capture_output=True,
        text=True,
    )

    rows = []

    for line in proc.stdout.splitlines():
        if not line:
            continue

        f = line.split("\t")
        qname = f[0]
        flag = int(f[1])
        pos = int(f[3])
        mapq = int(f[4])
        cigar = f[5]

        ref_pos = pos
        deletions = []
        matched_segments = []

        for length_s, op in CIGAR_RE.findall(cigar):
            length = int(length_s)

            if op in {"M", "=", "X"}:
                matched_segments.append(
                    (ref_pos, ref_pos + length - 1)
                )
                ref_pos += length

            elif op in {"D", "N"}:
                deletions.append(
                    (ref_pos, ref_pos + length - 1, length, op)
                )
                ref_pos += length

            elif op in REF_OPS:
                ref_pos += length

        overlapping = [
            event for event in deletions
            if event[0] <= args.focus_end
            and event[1] >= args.focus_start
        ]

        covering = [
            event for event in deletions
            if event[0] <= args.focus_start
            and event[1] >= args.focus_end
        ]

        matched_through = any(
            start <= args.focus_start and end >= args.focus_end
            for start, end in matched_segments
        )

        largest = max(overlapping, key=lambda x: x[2]) if overlapping else None

        rows.append({
            "qname": qname,
            "flag": flag,
            "mapq": mapq,
            "alignment_start": pos,
            "alignment_end": ref_pos - 1,
            "cigar": cigar,
            "matched_through_focus": matched_through,
            "deletion_overlaps_focus": bool(overlapping),
            "deletion_covers_focus": bool(covering),
            "largest_focus_deletion_length": largest[2] if largest else 0,
            "largest_focus_deletion_start": largest[0] if largest else "",
            "largest_focus_deletion_end": largest[1] if largest else "",
            "largest_focus_deletion_op": largest[3] if largest else "",
        })

    fields = list(rows[0].keys())
    per_read = Path(f"{out_prefix}.per_read.tsv")

    with per_read.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    lengths = Counter(
        row["largest_focus_deletion_length"]
        for row in rows
        if row["largest_focus_deletion_length"] > 0
    )

    histogram = Path(f"{out_prefix}.deletion_lengths.tsv")

    with histogram.open("w") as out:
        out.write("deletion_length\tread_count\n")
        for length, count in sorted(lengths.items()):
            out.write(f"{length}\t{count}\n")

    summary = Path(f"{out_prefix}.summary.txt")

    with summary.open("w") as out:
        out.write(f"bam\t{args.bam}\n")
        out.write(f"focus_start\t{args.focus_start}\n")
        out.write(f"focus_end\t{args.focus_end}\n")
        out.write(f"total_records\t{len(rows)}\n")
        out.write(
            "matched_through_focus\t"
            f"{sum(r['matched_through_focus'] for r in rows)}\n"
        )
        out.write(
            "deletion_overlaps_focus\t"
            f"{sum(r['deletion_overlaps_focus'] for r in rows)}\n"
        )
        out.write(
            "deletion_covers_focus\t"
            f"{sum(r['deletion_covers_focus'] for r in rows)}\n"
        )

    print(summary.read_text())
    print(f"Wrote: {per_read}")
    print(f"Wrote: {histogram}")


if __name__ == "__main__":
    main()
