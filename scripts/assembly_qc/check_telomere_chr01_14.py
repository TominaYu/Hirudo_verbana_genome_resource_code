#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import csv
from Bio import SeqIO

TEL_FWD = "TTAGGG"
TEL_REV = "CCCTAA"

# 端部サマリー用の長さ
END_WINDOW = 50000

# 端部スキャン用
SCAN_TOTAL = 200000     # 端から 200 kb を見る
SCAN_STEP = 10000       # 10 kb 窓
SCAN_WINDOW = 10000     # 10 kb 窓


def count_motif(seq: str, motif: str) -> int:
    return seq.count(motif)


def scan_left_end(seq: str, record_id: str):
    """
    左端から内側へ向かって SCAN_TOTAL bp を SCAN_WINDOW ごとにスキャン
    """
    rows = []
    n = len(seq)
    max_len = min(SCAN_TOTAL, n)

    start = 0
    while start < max_len:
        end = min(start + SCAN_WINDOW, max_len)
        frag = seq[start:end]
        rows.append({
            "record_id": record_id,
            "examined_end": "left",
            "window_index_from_end": start // SCAN_STEP + 1,
            "window_start_1based": start + 1,
            "window_end_1based": end,
            "TTAGGG_count": count_motif(frag, TEL_FWD),
            "CCCTAA_count": count_motif(frag, TEL_REV),
        })
        start += SCAN_STEP

    return rows


def scan_right_end(seq: str, record_id: str):
    """
    右端から内側へ向かって SCAN_TOTAL bp を SCAN_WINDOW ごとにスキャン
    """
    rows = []
    n = len(seq)
    max_len = min(SCAN_TOTAL, n)

    offset = 0
    while offset < max_len:
        end = n - offset
        start = max(0, end - SCAN_WINDOW)
        frag = seq[start:end]
        rows.append({
            "record_id": record_id,
            "examined_end": "right",
            "window_index_from_end": offset // SCAN_STEP + 1,
            "window_start_1based": start + 1,
            "window_end_1based": end,
            "TTAGGG_count": count_motif(frag, TEL_FWD),
            "CCCTAA_count": count_motif(frag, TEL_REV),
        })
        offset += SCAN_STEP

    return rows


def main():
    fasta = Path("05_output_v6/Hirudo_v6_nuclear_primary.fasta")
    outdir = Path("06_qc_v6")
    outdir.mkdir(parents=True, exist_ok=True)

    summary_tsv = outdir / "telomere_chr01_14_end_summary.tsv"
    scan_tsv = outdir / "telomere_chr01_14_end_scan.tsv"

    # 対象:
    # Chr01-13 -> 左右端
    # Chr14 -> Chr14_A1 の左端、Chr14_H3 の右端のみ
    target_both = {f"Chr{i:02d}" for i in range(1, 14)}
    target_special_left = {"Chr14_A1"}
    target_special_right = {"Chr14_H3"}

    summary_rows = []
    scan_rows = []

    found_ids = set()

    for rec in SeqIO.parse(str(fasta), "fasta"):
        rec_id = rec.id
        seq = str(rec.seq).upper()
        n = len(seq)

        if rec_id in target_both:
            found_ids.add(rec_id)

            left = seq[: min(END_WINDOW, n)]
            right = seq[max(0, n - END_WINDOW):]

            summary_rows.append({
                "record_id": rec_id,
                "examined_end": "left",
                "length_bp": n,
                "window_bp": min(END_WINDOW, n),
                "TTAGGG_count": count_motif(left, TEL_FWD),
                "CCCTAA_count": count_motif(left, TEL_REV),
            })
            summary_rows.append({
                "record_id": rec_id,
                "examined_end": "right",
                "length_bp": n,
                "window_bp": min(END_WINDOW, n),
                "TTAGGG_count": count_motif(right, TEL_FWD),
                "CCCTAA_count": count_motif(right, TEL_REV),
            })

            scan_rows.extend(scan_left_end(seq, rec_id))
            scan_rows.extend(scan_right_end(seq, rec_id))

        elif rec_id in target_special_left:
            found_ids.add(rec_id)

            left = seq[: min(END_WINDOW, n)]
            summary_rows.append({
                "record_id": rec_id,
                "examined_end": "left",
                "length_bp": n,
                "window_bp": min(END_WINDOW, n),
                "TTAGGG_count": count_motif(left, TEL_FWD),
                "CCCTAA_count": count_motif(left, TEL_REV),
            })
            scan_rows.extend(scan_left_end(seq, rec_id))

        elif rec_id in target_special_right:
            found_ids.add(rec_id)

            right = seq[max(0, n - END_WINDOW):]
            summary_rows.append({
                "record_id": rec_id,
                "examined_end": "right",
                "length_bp": n,
                "window_bp": min(END_WINDOW, n),
                "TTAGGG_count": count_motif(right, TEL_FWD),
                "CCCTAA_count": count_motif(right, TEL_REV),
            })
            scan_rows.extend(scan_right_end(seq, rec_id))

    # 出力
    with summary_tsv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "record_id",
                "examined_end",
                "length_bp",
                "window_bp",
                "TTAGGG_count",
                "CCCTAA_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    with scan_tsv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "record_id",
                "examined_end",
                "window_index_from_end",
                "window_start_1based",
                "window_end_1based",
                "TTAGGG_count",
                "CCCTAA_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(scan_rows)

    print(f"Wrote summary: {summary_tsv}")
    print(f"Wrote scan:    {scan_tsv}")

    expected_ids = target_both | target_special_left | target_special_right
    missing = sorted(expected_ids - found_ids)
    if missing:
        print("Warning: the following expected records were not found in FASTA:")
        for x in missing:
            print(f"  {x}")


if __name__ == "__main__":
    main()
