#!/usr/bin/env python3

import argparse
import csv
import re
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
QUERY_ALIGNED_OPS = {"M", "I", "=", "X"}
REFERENCE_OPS = {"M", "D", "N", "=", "X"}
MATCH_OPS = {"M", "=", "X"}


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--bam", required=True)
    parser.add_argument("--boundary", type=int, required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--out-prefix", required=True)

    return parser.parse_args()


def parse_cigar(cigar, reference_start):
    operations = [
        (int(length), op)
        for length, op in CIGAR_RE.findall(cigar)
    ]

    left_clip = 0
    for length, op in operations:
        if op in {"S", "H"}:
            left_clip += length
        else:
            break

    right_clip = 0
    for length, op in reversed(operations):
        if op in {"S", "H"}:
            right_clip += length
        else:
            break

    aligned_query_length = sum(
        length
        for length, op in operations
        if op in QUERY_ALIGNED_OPS
    )

    total_query_length = (
        left_clip
        + aligned_query_length
        + right_clip
    )

    reference_position = reference_start
    match_segments = []
    deletion_segments = []

    for length, op in operations:
        if op in MATCH_OPS:
            match_segments.append((
                reference_position,
                reference_position + length - 1,
            ))
            reference_position += length

        elif op in {"D", "N"}:
            deletion_segments.append((
                reference_position,
                reference_position + length - 1,
                length,
                op,
            ))
            reference_position += length

        elif op in REFERENCE_OPS:
            reference_position += length

    return {
        "left_clip": left_clip,
        "right_clip": right_clip,
        "aligned_query_length": aligned_query_length,
        "total_query_length": total_query_length,
        "reference_end": reference_position - 1,
        "match_segments": match_segments,
        "deletion_segments": deletion_segments,
    }


def match_covers(segments, position):
    return any(
        start <= position <= end
        for start, end in segments
    )


def match_crosses_boundary(segments, boundary):
    return any(
        start <= boundary
        and end >= boundary + 1
        for start, end in segments
    )


def deletion_crosses_boundary(segments, boundary):
    events = [
        event
        for event in segments
        if event[0] <= boundary
        and event[1] >= boundary + 1
    ]

    return events


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()

    bam = Path(args.bam)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    process = subprocess.run(
        ["samtools", "view", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )

    alignments = []

    for record_index, line in enumerate(
        process.stdout.splitlines(),
        start=1,
    ):
        if not line:
            continue

        fields = line.split("\t")

        qname = fields[0]
        flag = int(fields[1])
        reference_start = int(fields[3])
        mapq = int(fields[4])
        cigar = fields[5]

        if flag & 4:
            continue

        # Exclude secondary, retain primary and supplementary.
        if flag & 256:
            continue

        if mapq < args.min_mapq:
            continue

        if cigar == "*":
            continue

        parsed = parse_cigar(
            cigar,
            reference_start,
        )

        strand = "-" if flag & 16 else "+"

        if strand == "+":
            query_start = parsed["left_clip"]
            query_end = (
                parsed["left_clip"]
                + parsed["aligned_query_length"]
            )
        else:
            query_start = parsed["right_clip"]
            query_end = (
                parsed["right_clip"]
                + parsed["aligned_query_length"]
            )

        match_segments = parsed["match_segments"]
        deletion_segments = parsed["deletion_segments"]

        deletion_events = deletion_crosses_boundary(
            deletion_segments,
            args.boundary,
        )

        alignment = {
            "alignment_index": record_index,
            "qname": qname,
            "flag": flag,
            "alignment_type": (
                "supplementary"
                if flag & 2048
                else "primary"
            ),
            "strand": strand,
            "mapq": mapq,
            "reference_start": reference_start,
            "reference_end": parsed["reference_end"],
            "query_start_0based": query_start,
            "query_end_0based": query_end,
            "query_aligned_length": parsed["aligned_query_length"],
            "cigar": cigar,
            "match_covers_left_boundary_base": match_covers(
                match_segments,
                args.boundary,
            ),
            "match_covers_right_boundary_base": match_covers(
                match_segments,
                args.boundary + 1,
            ),
            "continuous_match_crosses_boundary": (
                match_crosses_boundary(
                    match_segments,
                    args.boundary,
                )
            ),
            "deletion_crosses_boundary": bool(
                deletion_events
            ),
            "largest_boundary_deletion": (
                max(
                    event[2]
                    for event in deletion_events
                )
                if deletion_events
                else 0
            ),
        }

        alignments.append(alignment)

    if not alignments:
        raise SystemExit("No alignments passed the filters.")

    alignment_fields = list(alignments[0].keys())

    write_tsv(
        Path(f"{out_prefix}.alignments.tsv"),
        alignments,
        alignment_fields,
    )

    by_read = defaultdict(list)

    for alignment in alignments:
        by_read[alignment["qname"]].append(alignment)

    read_rows = []
    split_rows = []

    for qname, records in sorted(by_read.items()):
        continuous = [
            record
            for record in records
            if record["continuous_match_crosses_boundary"]
        ]

        deletion_crossing = [
            record
            for record in records
            if record["deletion_crosses_boundary"]
        ]

        left_records = [
            record
            for record in records
            if record["match_covers_left_boundary_base"]
        ]

        right_records = [
            record
            for record in records
            if record["match_covers_right_boundary_base"]
        ]

        best_pair = None

        for left in left_records:
            for right in right_records:
                if (
                    left["alignment_index"]
                    == right["alignment_index"]
                ):
                    continue

                same_strand = (
                    left["strand"] == right["strand"]
                )

                if not same_strand:
                    query_gap = None
                    score = (
                        1,
                        10**12,
                        -left["mapq"] - right["mapq"],
                    )
                elif left["strand"] == "+":
                    query_gap = (
                        right["query_start_0based"]
                        - left["query_end_0based"]
                    )
                    score = (
                        0,
                        abs(query_gap),
                        -left["mapq"] - right["mapq"],
                    )
                else:
                    query_gap = (
                        left["query_start_0based"]
                        - right["query_end_0based"]
                    )
                    score = (
                        0,
                        abs(query_gap),
                        -left["mapq"] - right["mapq"],
                    )

                candidate = {
                    "qname": qname,
                    "same_strand": same_strand,
                    "strand": (
                        left["strand"]
                        if same_strand
                        else "mixed"
                    ),
                    "query_gap_bp": (
                        query_gap
                        if query_gap is not None
                        else ""
                    ),
                    "left_alignment_type": (
                        left["alignment_type"]
                    ),
                    "left_reference_start": (
                        left["reference_start"]
                    ),
                    "left_reference_end": (
                        left["reference_end"]
                    ),
                    "left_query_start_0based": (
                        left["query_start_0based"]
                    ),
                    "left_query_end_0based": (
                        left["query_end_0based"]
                    ),
                    "left_mapq": left["mapq"],
                    "left_cigar": left["cigar"],
                    "right_alignment_type": (
                        right["alignment_type"]
                    ),
                    "right_reference_start": (
                        right["reference_start"]
                    ),
                    "right_reference_end": (
                        right["reference_end"]
                    ),
                    "right_query_start_0based": (
                        right["query_start_0based"]
                    ),
                    "right_query_end_0based": (
                        right["query_end_0based"]
                    ),
                    "right_mapq": right["mapq"],
                    "right_cigar": right["cigar"],
                    "_score": score,
                }

                if (
                    best_pair is None
                    or candidate["_score"]
                    < best_pair["_score"]
                ):
                    best_pair = candidate

        if continuous:
            classification = "continuous_single_alignment"

        elif deletion_crossing:
            classification = "deletion_crossing_single_alignment"

        elif best_pair is not None and best_pair["same_strand"]:
            classification = "split_alignment_pair"

        elif best_pair is not None:
            classification = "split_pair_mixed_orientation"

        else:
            classification = "one_sided_or_unresolved"

        largest_deletion = max(
            [
                record["largest_boundary_deletion"]
                for record in deletion_crossing
            ],
            default=0,
        )

        read_rows.append({
            "qname": qname,
            "classification": classification,
            "retained_alignment_count": len(records),
            "primary_alignment_count": sum(
                record["alignment_type"] == "primary"
                for record in records
            ),
            "supplementary_alignment_count": sum(
                record["alignment_type"] == "supplementary"
                for record in records
            ),
            "continuous_alignment_count": len(continuous),
            "deletion_crossing_alignment_count": len(
                deletion_crossing
            ),
            "largest_boundary_deletion": largest_deletion,
            "has_best_split_pair": best_pair is not None,
            "best_split_same_strand": (
                best_pair["same_strand"]
                if best_pair is not None
                else ""
            ),
            "best_split_strand": (
                best_pair["strand"]
                if best_pair is not None
                else ""
            ),
            "best_split_query_gap_bp": (
                best_pair["query_gap_bp"]
                if best_pair is not None
                else ""
            ),
        })

        if best_pair is not None:
            best_pair.pop("_score")
            split_rows.append(best_pair)

    read_fields = list(read_rows[0].keys())

    write_tsv(
        Path(f"{out_prefix}.per_read.tsv"),
        read_rows,
        read_fields,
    )

    if split_rows:
        split_fields = list(split_rows[0].keys())

        write_tsv(
            Path(f"{out_prefix}.best_split_pairs.tsv"),
            split_rows,
            split_fields,
        )

    classification_counts = defaultdict(int)

    for row in read_rows:
        classification_counts[row["classification"]] += 1

    same_strand_gaps = [
        int(row["best_split_query_gap_bp"])
        for row in read_rows
        if row["classification"] == "split_alignment_pair"
        and row["best_split_query_gap_bp"] != ""
    ]

    summary_path = Path(f"{out_prefix}.summary.txt")

    with summary_path.open("w") as out:
        out.write(f"bam\t{bam}\n")
        out.write(f"boundary\t{args.boundary}\n")
        out.write(f"min_mapq\t{args.min_mapq}\n")
        out.write(
            "secondary_alignments_excluded\tYES\n"
        )
        out.write(
            "supplementary_alignments_retained\tYES\n"
        )
        out.write(
            f"retained_alignment_records\t"
            f"{len(alignments)}\n"
        )
        out.write(
            f"unique_reads\t{len(read_rows)}\n"
        )

        for classification in [
            "continuous_single_alignment",
            "deletion_crossing_single_alignment",
            "split_alignment_pair",
            "split_pair_mixed_orientation",
            "one_sided_or_unresolved",
        ]:
            out.write(
                f"{classification}\t"
                f"{classification_counts[classification]}\n"
            )

        if same_strand_gaps:
            out.write(
                f"split_query_gap_min_bp\t"
                f"{min(same_strand_gaps)}\n"
            )
            out.write(
                f"split_query_gap_median_bp\t"
                f"{statistics.median(same_strand_gaps)}\n"
            )
            out.write(
                f"split_query_gap_max_bp\t"
                f"{max(same_strand_gaps)}\n"
            )

    qname_groups = defaultdict(list)

    for row in read_rows:
        qname_groups[row["classification"]].append(
            row["qname"]
        )

    for classification, qnames in qname_groups.items():
        qname_path = Path(
            f"{out_prefix}.{classification}.qnames.txt"
        )

        with qname_path.open("w") as out:
            for qname in qnames:
                out.write(qname + "\n")

    print(summary_path.read_text())

    print(
        f"Wrote: {out_prefix}.alignments.tsv"
    )
    print(
        f"Wrote: {out_prefix}.per_read.tsv"
    )

    if split_rows:
        print(
            f"Wrote: {out_prefix}.best_split_pairs.tsv"
        )


if __name__ == "__main__":
    main()
