#!/usr/bin/env python3
import csv
from collections import Counter, defaultdict
from pathlib import Path

targets = {
    "our_primary_lophotrochozoa": Path("07_finalize_v8_20260327_232659/09_busco6_setup_20260328_002021/04_busco_runs/busco6_primary_lophotrochozoa/run_lophotrochozoa_odb12/full_table.tsv"),
    "eu_lophotrochozoa": Path("07_finalize_v8_20260327_232659/09_busco6_setup_20260328_002021/04_busco_runs/busco6_eu_lophotrochozoa/run_lophotrochozoa_odb12/full_table.tsv"),
}

outdir = Path("07_finalize_v8_20260327_232659/09_busco6_setup_20260328_002021/05_compare_tables")
outdir.mkdir(parents=True, exist_ok=True)

def parse_full_table(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            rows.append({
                "busco_id": fields[0] if len(fields) > 0 else "",
                "status": fields[1] if len(fields) > 1 else "",
                "sequence": fields[2] if len(fields) > 2 else "",
                "start": fields[3] if len(fields) > 3 else "",
                "end": fields[4] if len(fields) > 4 else "",
                "score": fields[5] if len(fields) > 5 else "",
                "length": fields[6] if len(fields) > 6 else "",
                "description": fields[8] if len(fields) > 8 else "",
            })
    return rows

for label, path in targets.items():
    if not path.exists():
        raise FileNotFoundError(path)

    rows = parse_full_table(path)

    non_complete = [r for r in rows if r["status"] != "Complete"]
    out_noncomplete = outdir / f"{label}.non_complete.tsv"
    with open(out_noncomplete, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["busco_id", "status", "sequence", "start", "end", "score", "length", "description"],
            delimiter="\t"
        )
        w.writeheader()
        w.writerows(non_complete)

    per_seq = defaultdict(Counter)
    for r in rows:
        seq = r["sequence"] if r["sequence"] else "NA"
        per_seq[seq][r["status"]] += 1

    out_perseq = outdir / f"{label}.status_by_sequence.tsv"
    with open(out_perseq, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sequence", "Complete", "Duplicated", "Fragmented", "Missing", "Total_non_complete"])
        for seq, c in sorted(per_seq.items(), key=lambda x: (x[0] == "NA", x[0])):
            total_non_complete = c["Duplicated"] + c["Fragmented"] + c["Missing"]
            w.writerow([seq, c["Complete"], c["Duplicated"], c["Fragmented"], c["Missing"], total_non_complete])

    print(f"[{label}]")
    print(f"  full_table: {path}")
    print(f"  non_complete_tsv: {out_noncomplete}")
    print(f"  status_by_sequence_tsv: {out_perseq}")
    print(f"  rows_total: {len(rows)}")
    print(f"  rows_non_complete: {len(non_complete)}")
