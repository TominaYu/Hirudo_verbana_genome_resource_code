#!/usr/bin/env python3
import json
import csv
from pathlib import Path

rundir = Path("07_finalize_v8_20260327_232659/09_busco6_setup_20260328_002021/04_busco_runs")
out_tsv = Path("07_finalize_v8_20260327_232659/09_busco6_setup_20260328_002021/05_compare_tables/busco_summary_comparison.tsv")

targets = {
    ("our_primary", "eukaryota_odb10"): rundir / "busco6_primary_eukaryota_odb10" / "short_summary.specific.eukaryota_odb10.busco6_primary_eukaryota_odb10.json",
    ("our_primary", "lophotrochozoa_odb12"): rundir / "busco6_primary_lophotrochozoa" / "short_summary.specific.lophotrochozoa_odb12.busco6_primary_lophotrochozoa.json",
    ("our_combined", "lophotrochozoa_odb12"): rundir / "busco6_combined_lophotrochozoa" / "short_summary.specific.lophotrochozoa_odb12.busco6_combined_lophotrochozoa.json",
    ("eu_assembly", "eukaryota_odb10"): rundir / "busco6_eu_eukaryota_odb10" / "short_summary.specific.eukaryota_odb10.busco6_eu_eukaryota_odb10.json",
    ("eu_assembly", "lophotrochozoa_odb12"): rundir / "busco6_eu_lophotrochozoa" / "short_summary.specific.lophotrochozoa_odb12.busco6_eu_lophotrochozoa.json",
}

rows = []
for (assembly, lineage), path in targets.items():
    if not path.exists():
        print(f"WARNING: missing {path}")
        continue
    with open(path) as f:
        d = json.load(f)
    r = d["results"]
    m = d["metrics"]
    rows.append({
        "assembly": assembly,
        "lineage": lineage,
        "C_pct": r["Complete percentage"],
        "S_pct": r["Single copy percentage"],
        "D_pct": r["Multi copy percentage"],
        "F_pct": r["Fragmented percentage"],
        "M_pct": r["Missing percentage"],
        "n_markers": r["n_markers"],
        "complete_n": r["Complete BUSCOs"],
        "single_n": r["Single copy BUSCOs"],
        "dup_n": r["Multi copy BUSCOs"],
        "frag_n": r["Fragmented BUSCOs"],
        "miss_n": r["Missing BUSCOs"],
        "internal_stop_pct": r.get("internal_stop_codon_percent", "NA"),
        "n_scaffolds": m["Number of scaffolds"],
        "total_length": m["Total length"],
        "scaffold_N50": m["Scaffold N50"],
    })

if not rows:
    raise SystemExit("No BUSCO summary JSON files found.")

with open(out_tsv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"OUTPUT: {out_tsv}")
print(out_tsv.read_text())
