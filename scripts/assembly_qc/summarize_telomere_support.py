#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import pandas as pd


def classify_support(ttaggg: int, ccctaa: int) -> tuple[str, str, float]:
    """
    Return:
      dominant_motif, support_level, bias_ratio

    bias_ratio = dominant / (minor + 1)
    """
    if ttaggg >= ccctaa:
        dominant = "TTAGGG"
        major = ttaggg
        minor = ccctaa
    else:
        dominant = "CCCTAA"
        major = ccctaa
        minor = ttaggg

    bias_ratio = major / (minor + 1)

    # Heuristic thresholds
    if major >= 300 and bias_ratio >= 5:
        support = "strong"
    elif major >= 150 and bias_ratio >= 3:
        support = "moderate"
    elif major >= 50 and bias_ratio >= 1.5:
        support = "weak"
    else:
        support = "ambiguous"

    return dominant, support, bias_ratio


def expected_end_direction(record_id: str, examined_end: str) -> str:
    """
    Expected pattern under the observed genome-wide tendency:
      left end  -> CCCTAA-rich
      right end -> TTAGGG-rich

    Chr14 special case is already encoded in the input summary
    (only Chr14_A1 left and Chr14_H3 right are present).
    """
    if examined_end == "left":
        return "CCCTAA"
    elif examined_end == "right":
        return "TTAGGG"
    return "unknown"


def main():
    inp = Path("06_qc_v6/telomere_chr01_14_end_summary.tsv")
    out = Path("06_qc_v6/telomere_chr01_14_support_summary.tsv")

    df = pd.read_csv(inp, sep="\t")

    rows = []
    for _, row in df.iterrows():
        record_id = str(row["record_id"])
        examined_end = str(row["examined_end"])
        ttaggg = int(row["TTAGGG_count"])
        ccctaa = int(row["CCCTAA_count"])

        dominant, support, bias_ratio = classify_support(ttaggg, ccctaa)
        expected = expected_end_direction(record_id, examined_end)
        matches_expected = (dominant == expected)

        rows.append({
            "record_id": record_id,
            "examined_end": examined_end,
            "TTAGGG_count": ttaggg,
            "CCCTAA_count": ccctaa,
            "dominant_motif": dominant,
            "expected_motif_for_end": expected,
            "matches_expected": matches_expected,
            "bias_ratio": round(bias_ratio, 2),
            "support_level": support,
        })

    out_df = pd.DataFrame(rows)

    # sort nicely
    def sort_key(rec: str):
        if rec.startswith("Chr14_A1"):
            return (14, "A1")
        if rec.startswith("Chr14_H3"):
            return (14, "H3")
        if rec.startswith("Chr"):
            num = rec.replace("Chr", "")
            try:
                return (int(num), "")
            except ValueError:
                return (999, rec)
        return (999, rec)

    out_df["_sort"] = out_df["record_id"].map(sort_key)
    out_df = out_df.sort_values(by=["_sort", "examined_end"]).drop(columns=["_sort"])

    out_df.to_csv(out, sep="\t", index=False)

    print(f"Wrote: {out}")

    print("\n=== Summary by support_level ===")
    print(out_df["support_level"].value_counts())

    print("\n=== Ends not matching expected motif direction ===")
    mismatch = out_df[~out_df["matches_expected"]]
    if len(mismatch) == 0:
        print("None")
    else:
        print(mismatch.to_string(index=False))


if __name__ == "__main__":
    main()
