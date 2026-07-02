#!/usr/bin/env python3

import csv
import statistics
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = Path(
    "project_primary_FASTA/"
    "12_suppfig2_read_support/06_panel_E"
)

NO_GAP_TSV = Path(
    "project_primary_FASTA/"
    "12_suppfig2_read_support/05_igv/"
    "00_diagnostics/D4_noGap_focus.per_read.tsv"
)

GAP1K_TSV = Path(
    "project_primary_FASTA/"
    "12_suppfig2_read_support/05_igv/"
    "00_diagnostics/D5_gap1k_focus.per_read.tsv"
)

TABLEDIR = BASE / "01_tables"
CHARTDIR = BASE / "04_charts"

TABLEDIR.mkdir(parents=True, exist_ok=True)
CHARTDIR.mkdir(parents=True, exist_ok=True)


def parse_bool(value):
    return str(value).strip().lower() in {
        "true", "1", "yes"
    }


def read_unique_rows(path):
    rows = {}

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            qname = row["qname"]

            if qname in rows:
                raise ValueError(
                    f"Duplicate qname in {path}: {qname}"
                )

            row["deletion_covers_focus"] = parse_bool(
                row["deletion_covers_focus"]
            )

            row["deletion_overlaps_focus"] = parse_bool(
                row["deletion_overlaps_focus"]
            )

            row["largest_focus_deletion_length"] = int(
                row["largest_focus_deletion_length"]
            )

            rows[qname] = row

    return rows


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


no_gap_all = read_unique_rows(NO_GAP_TSV)
gap1k_all = read_unique_rows(GAP1K_TSV)

no_gap_covering = {
    qname: row
    for qname, row in no_gap_all.items()
    if row["deletion_covers_focus"]
}

gap1k_covering = {
    qname: row
    for qname, row in gap1k_all.items()
    if row["deletion_covers_focus"]
}

paired_qnames = sorted(
    set(no_gap_covering)
    & set(gap1k_covering)
)

no_gap_only = sorted(
    set(no_gap_covering)
    - set(gap1k_covering)
)

gap1k_only = sorted(
    set(gap1k_covering)
    - set(no_gap_covering)
)

if not paired_qnames:
    raise SystemExit(
        "No read identities were shared between the two models."
    )

paired_rows = []

for qname in paired_qnames:
    no_length = no_gap_covering[qname][
        "largest_focus_deletion_length"
    ]

    gap_length = gap1k_covering[qname][
        "largest_focus_deletion_length"
    ]

    paired_rows.append({
        "qname": qname,
        "no_spacer_deletion_length_bp": no_length,
        "gap1k_deletion_length_bp": gap_length,
        "difference_bp": gap_length - no_length,
    })

paired_tsv = TABLEDIR / "PanelE_paired_deletion_lengths.tsv"

write_tsv(
    paired_tsv,
    paired_rows,
    [
        "qname",
        "no_spacer_deletion_length_bp",
        "gap1k_deletion_length_bp",
        "difference_bp",
    ],
)

paired_qname_file = TABLEDIR / "PanelE_paired_qnames.txt"

with paired_qname_file.open("w") as out:
    for qname in paired_qnames:
        out.write(qname + "\n")

with (
    TABLEDIR / "PanelE_no_spacer_only_qnames.txt"
).open("w") as out:
    for qname in no_gap_only:
        out.write(qname + "\n")

with (
    TABLEDIR / "PanelE_gap1k_only_qnames.txt"
).open("w") as out:
    for qname in gap1k_only:
        out.write(qname + "\n")

differences = [
    row["difference_bp"]
    for row in paired_rows
]

no_lengths = [
    row["no_spacer_deletion_length_bp"]
    for row in paired_rows
]

gap_lengths = [
    row["gap1k_deletion_length_bp"]
    for row in paired_rows
]

difference_counts = Counter(differences)

summary_file = TABLEDIR / "PanelE_paired_summary.txt"

with summary_file.open("w") as out:
    out.write(
        f"no_spacer_deletion_covers_focus\t"
        f"{len(no_gap_covering)}\n"
    )

    out.write(
        f"gap1k_deletion_covers_focus\t"
        f"{len(gap1k_covering)}\n"
    )

    out.write(
        f"paired_read_identities\t"
        f"{len(paired_qnames)}\n"
    )

    out.write(
        f"no_spacer_only\t{len(no_gap_only)}\n"
    )

    out.write(
        f"gap1k_only\t{len(gap1k_only)}\n"
    )

    out.write(
        f"no_spacer_length_min_bp\t"
        f"{min(no_lengths)}\n"
    )

    out.write(
        f"no_spacer_length_median_bp\t"
        f"{statistics.median(no_lengths)}\n"
    )

    out.write(
        f"no_spacer_length_max_bp\t"
        f"{max(no_lengths)}\n"
    )

    out.write(
        f"gap1k_length_min_bp\t"
        f"{min(gap_lengths)}\n"
    )

    out.write(
        f"gap1k_length_median_bp\t"
        f"{statistics.median(gap_lengths)}\n"
    )

    out.write(
        f"gap1k_length_max_bp\t"
        f"{max(gap_lengths)}\n"
    )

    out.write(
        f"difference_min_bp\t"
        f"{min(differences)}\n"
    )

    out.write(
        f"difference_median_bp\t"
        f"{statistics.median(differences)}\n"
    )

    out.write(
        f"difference_max_bp\t"
        f"{max(differences)}\n"
    )

    out.write(
        "difference_distribution\t"
        + ",".join(
            f"{length}:{count}"
            for length, count
            in sorted(difference_counts.items())
        )
        + "\n"
    )

# Paired plot
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

fig, ax = plt.subplots(figsize=(3.8, 4.2))

for row in paired_rows:
    y0 = row["no_spacer_deletion_length_bp"]
    y1 = row["gap1k_deletion_length_bp"]

    ax.plot(
        [0, 1],
        [y0, y1],
        color="#B8B8B8",
        linewidth=0.8,
        alpha=0.75,
        zorder=1,
    )

ax.scatter(
    [0] * len(no_lengths),
    no_lengths,
    color="#0072B2",
    edgecolor="white",
    linewidth=0.4,
    s=30,
    zorder=3,
)

ax.scatter(
    [1] * len(gap_lengths),
    gap_lengths,
    color="#8F8F8F",
    edgecolor="white",
    linewidth=0.4,
    s=30,
    zorder=3,
)

ax.set_xlim(-0.35, 1.35)
ax.set_xticks([0, 1])
ax.set_xticklabels([
    "No N spacer",
    "1-kb N spacer",
])

ax.set_ylabel(
    "Deletion-like CIGAR length (bp)"
)

ax.set_title(
    "Same-read comparison",
    fontsize=9,
    fontweight="bold",
)

ax.grid(
    axis="y",
    color="#D9D9D9",
    linewidth=0.45,
    zorder=0,
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(0.6)
ax.spines["bottom"].set_linewidth(0.6)

median_difference = statistics.median(differences)

ax.text(
    0.5,
    0.97,
    (
        f"n={len(paired_rows)} paired reads\n"
        f"median Δ={median_difference:g} bp"
    ),
    transform=ax.transAxes,
    ha="center",
    va="top",
    fontsize=7,
    color="#444444",
)

fig.tight_layout()

for suffix, kwargs in [
    ("pdf", {}),
    ("svg", {}),
    ("png", {"dpi": 600}),
]:
    fig.savefig(
        CHARTDIR
        / f"PanelE_same_read_deletion_shift.{suffix}",
        bbox_inches="tight",
        **kwargs,
    )

plt.close(fig)

print(summary_file.read_text())
print(f"Wrote: {paired_tsv}")
print(f"Wrote: {paired_qname_file}")
print(
    "Wrote: "
    f"{CHARTDIR}/PanelE_same_read_deletion_shift.pdf"
)
