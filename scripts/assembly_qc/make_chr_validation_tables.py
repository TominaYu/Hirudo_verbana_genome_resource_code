#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path
from collections import defaultdict

import pandas as pd


def normalize_component_name(name: str) -> str:
    name = str(name).strip()
    if "_length_" in name:
        name = name.split("_length_")[0]
    return name


def main():
    project_dir = Path(".")
    recon_tsv = project_dir / "05_output" / "Hirudo_reconstruction_map.tsv"
    outdir = project_dir / "06_qc"
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(recon_tsv, sep="\t")

    # Chr01–Chr13 のみ
    chr_df = df[
        (df["set_name"] == "primary") &
        (df["assembly_group"].astype(str).str.match(r"^Chr([1-9]|0[1-9]|1[0-3])$|^Chr0?[1-9]$|^Chr1[0-3]$"))
    ].copy()

    # Chr名を自然順に並べるためのキー
    def chr_num(x: str) -> int:
        return int(str(x).replace("Chr", "").lstrip("0") or "0")

    chr_df["chr_num"] = chr_df["assembly_group"].apply(chr_num)
    chr_df = chr_df.sort_values(["chr_num", "order_in_group", "component_beg", "component_end"])

    # 1) Chr構成表
    chr_components_tsv = outdir / "chr01_13_component_table.tsv"
    chr_df.to_csv(chr_components_tsv, sep="\t", index=False)

    # 2) split contig 抽出
    split_rows = []
    for chr_name, sub in chr_df.groupby("assembly_group"):
        counts = sub["component_normalized"].value_counts()
        duplicated = counts[counts > 1].index.tolist()
        if duplicated:
            for comp in duplicated:
                sub2 = sub[sub["component_normalized"] == comp].sort_values(
                    ["order_in_group", "component_beg", "component_end"]
                )
                for _, row in sub2.iterrows():
                    split_rows.append({
                        "chromosome": chr_name,
                        "component_original": row["component_original"],
                        "component_normalized": row["component_normalized"],
                        "component_beg": row["component_beg"],
                        "component_end": row["component_end"],
                        "orientation_curated": row["orientation_curated"],
                        "order_in_group": row["order_in_group"],
                        "length_bp": row["length_bp"],
                    })

    split_tsv = outdir / "chr01_13_split_components.tsv"
    pd.DataFrame(split_rows).to_csv(split_tsv, sep="\t", index=False)

    # 3) Chr summary
    summary_rows = []
    for chr_name, sub in chr_df.groupby("assembly_group"):
        n_components = len(sub)
        n_unique = sub["component_normalized"].nunique()
        n_reverse = (sub["orientation_curated"] == "-").sum()
        has_split = n_unique < n_components
        summary_rows.append({
            "chromosome": chr_name,
            "n_components": n_components,
            "n_unique_components": n_unique,
            "n_reverse_orientation": int(n_reverse),
            "has_split_component": has_split,
        })

    summary_tsv = outdir / "chr01_13_validation_summary.tsv"
    pd.DataFrame(summary_rows).sort_values(
        by="chromosome",
        key=lambda s: s.str.replace("Chr", "").astype(int)
    ).to_csv(summary_tsv, sep="\t", index=False)

    print(f"Wrote: {chr_components_tsv}")
    print(f"Wrote: {split_tsv}")
    print(f"Wrote: {summary_tsv}")


if __name__ == "__main__":
    main()
