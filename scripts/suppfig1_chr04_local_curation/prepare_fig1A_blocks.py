#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv

PROJECT = Path("project_compare_EU_Hverb")
SPLIT_FASTA = PROJECT / "01_query/Hirudo_v6_Chr14_split.fasta"
COORDS = PROJECT / "08_suppfig1A_rebuild/01_tables/Chr14_vs_EUchr2.1delta.coords.tsv"
OUTDIR = PROJECT / "08_suppfig1A_rebuild/01_tables"

FINAL_FASTA = Path(
    "project_primary_FASTA/"
    "07_finalize_v8_20260327_232659/"
    "05_outputs/reconstruct_v8_reconstruct_compat/"
    "Hirudo_v8_nuclear_primary.fasta"
)

COMPONENTS = [
    ("A1", "Chr14_A1", 1, 6435254),
    ("A2", "Chr14_A2", 6435255, 7500755),
    ("E",  "Chr14_E",  7500756, 9013138),
    ("G",  "Chr14_G",  9013139, 9382960),
    ("H1", "Chr14_H1", 9382961, 12299010),
    ("H3", "Chr14_H3", 12299011, 15289010),
]

TARGET_ID = "OZ237780.1"
TARGET_LENGTH = 16479789


def read_fasta(path, wanted=None):
    records = {}
    name = None
    chunks = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    records[name] = "".join(chunks).upper()

                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

        if name is not None and (wanted is None or name in wanted):
            records[name] = "".join(chunks).upper()

    return records


def revcomp(seq):
    table = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return seq.translate(table)[::-1]


def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [list(intervals[0])]

    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])

    return [tuple(x) for x in merged]


def union_length(intervals):
    return sum(end - start + 1 for start, end in merge_intervals(intervals))


OUTDIR.mkdir(parents=True, exist_ok=True)

split_names = {record_id for _, record_id, _, _ in COMPONENTS}
split_records = read_fasta(SPLIT_FASTA, split_names)
final_records = read_fasta(FINAL_FASTA, {"Chr04"})

if "Chr04" not in final_records:
    raise SystemExit("Chr04 was not found in the finalized primary FASTA.")

chr04 = final_records["Chr04"]

if len(chr04) != 15289010:
    raise SystemExit(
        f"Unexpected Chr04 length: {len(chr04)}; expected 15289010."
    )

orientation_rows = []
orientation_by_record = {}

for label, record_id, final_start, final_end in COMPONENTS:
    if record_id not in split_records:
        raise SystemExit(f"Missing split FASTA record: {record_id}")

    query_seq = split_records[record_id]
    final_segment = chr04[final_start - 1:final_end]

    if len(query_seq) != final_end - final_start + 1:
        state = "LENGTH_MISMATCH"
    elif query_seq == final_segment:
        state = "direct"
    elif revcomp(query_seq) == final_segment:
        state = "reverse"
    else:
        state = "NO_EXACT_MATCH"

    orientation_by_record[record_id] = state

    orientation_rows.append({
        "label": label,
        "record_id": record_id,
        "record_length": len(query_seq),
        "chr04_start_1based": final_start,
        "chr04_end_1based": final_end,
        "orientation_in_final_chr04": state,
    })

orientation_file = OUTDIR / "fig1A_component_orientation_check.tsv"

with orientation_file.open("w", newline="") as out:
    fields = list(orientation_rows[0].keys())
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(orientation_rows)

bad = [
    row for row in orientation_rows
    if row["orientation_in_final_chr04"] not in {"direct", "reverse"}
]

if bad:
    print(orientation_file.read_text())
    raise SystemExit(
        "At least one component did not exactly match the corresponding "
        "final Chr04 segment. Plot preparation was stopped."
    )

component_info = {
    record_id: {
        "label": label,
        "final_start": final_start,
        "final_end": final_end,
        "query_length": final_end - final_start + 1,
        "orientation": orientation_by_record[record_id],
    }
    for label, record_id, final_start, final_end in COMPONENTS
}

blocks = []

with COORDS.open() as handle:
    for line in handle:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        fields = line.split()

        if len(fields) < 13:
            continue

        s1 = int(fields[0])
        e1 = int(fields[1])
        s2 = int(fields[2])
        e2 = int(fields[3])
        len1 = int(fields[4])
        len2 = int(fields[5])
        identity = float(fields[6])
        len_ref = int(fields[7])
        len_query = int(fields[8])
        cov_ref = float(fields[9])
        cov_query = float(fields[10])
        ref_id = fields[11]
        query_id = fields[12]

        if ref_id != TARGET_ID or query_id not in component_info:
            continue

        info = component_info[query_id]
        query_orientation = info["orientation"]
        query_length = info["query_length"]
        offset0 = info["final_start"] - 1

        def query_pos_to_chr04(pos1):
            if query_orientation == "direct":
                return offset0 + pos1
            return offset0 + (query_length - pos1 + 1)

        chr04_at_s2 = query_pos_to_chr04(s2)
        chr04_at_e2 = query_pos_to_chr04(e2)

        chr04_start = min(chr04_at_s2, chr04_at_e2)
        chr04_end = max(chr04_at_s2, chr04_at_e2)

        nucmer_query_strand = "+" if s2 <= e2 else "-"

        if query_orientation == "direct":
            final_strand_vs_target = nucmer_query_strand
        else:
            final_strand_vs_target = (
                "-" if nucmer_query_strand == "+" else "+"
            )

        blocks.append({
            "label": info["label"],
            "query_id": query_id,
            "query_length": query_length,
            "query_start_1based": min(s2, e2),
            "query_end_1based": max(s2, e2),
            "query_strand_vs_target": nucmer_query_strand,
            "orientation_in_final_chr04": query_orientation,
            "chr04_start_1based": chr04_start,
            "chr04_end_1based": chr04_end,
            "final_chr04_strand_vs_target": final_strand_vs_target,
            "target_id": ref_id,
            "target_length": len_ref,
            "target_start_1based": min(s1, e1),
            "target_end_1based": max(s1, e1),
            "reference_aligned_length": len1,
            "query_aligned_length": len2,
            "minimum_aligned_length": min(len1, len2),
            "identity_pct": identity,
            "reference_coverage_pct": cov_ref,
            "query_coverage_pct": cov_query,
        })

if not blocks:
    raise SystemExit("No selected NUCmer blocks were found.")

block_fields = list(blocks[0].keys())
all_blocks_file = OUTDIR / "fig1A_backbone6_nucmer_blocks_all.tsv"

with all_blocks_file.open("w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=block_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(
        sorted(
            blocks,
            key=lambda row: (
                row["chr04_start_1based"],
                row["target_start_1based"],
            ),
        )
    )

thresholds = [
    (1000, 90.0),
    (5000, 90.0),
    (10000, 90.0),
    (20000, 90.0),
    (50000, 90.0),
    (1000, 95.0),
    (5000, 95.0),
    (10000, 95.0),
    (20000, 95.0),
    (50000, 95.0),
    (1000, 97.0),
    (5000, 97.0),
    (10000, 97.0),
    (20000, 97.0),
    (50000, 97.0),
]

summary_rows = []

for min_length, min_identity in thresholds:
    selected = [
        row for row in blocks
        if row["minimum_aligned_length"] >= min_length
        and row["identity_pct"] >= min_identity
    ]

    selected_file = (
        OUTDIR
        / f"fig1A_blocks_min{min_length}_id{int(min_identity)}.tsv"
    )

    with selected_file.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=block_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            sorted(
                selected,
                key=lambda row: (
                    row["chr04_start_1based"],
                    row["target_start_1based"],
                ),
            )
        )

    for label, record_id, final_start, final_end in COMPONENTS:
        subset = [row for row in selected if row["query_id"] == record_id]

        query_intervals = [
            (
                row["query_start_1based"],
                row["query_end_1based"],
            )
            for row in subset
        ]

        target_intervals = [
            (
                row["target_start_1based"],
                row["target_end_1based"],
            )
            for row in subset
        ]

        query_length = final_end - final_start + 1
        query_covered = union_length(query_intervals)
        target_covered = union_length(target_intervals)

        summary_rows.append({
            "minimum_block_length": min_length,
            "minimum_identity_pct": min_identity,
            "label": label,
            "query_id": record_id,
            "query_length": query_length,
            "block_count": len(subset),
            "query_union_covered_bp": query_covered,
            "query_union_coverage_pct": (
                100.0 * query_covered / query_length
                if query_length else 0.0
            ),
            "target_union_covered_bp": target_covered,
            "forward_block_count": sum(
                row["final_chr04_strand_vs_target"] == "+"
                for row in subset
            ),
            "reverse_block_count": sum(
                row["final_chr04_strand_vs_target"] == "-"
                for row in subset
            ),
        })

summary_file = OUTDIR / "fig1A_nucmer_filter_summary.tsv"

with summary_file.open("w", newline="") as out:
    fields = list(summary_rows[0].keys())
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"Wrote: {orientation_file}")
print(f"Wrote: {all_blocks_file}")
print(f"Wrote: {summary_file}")
