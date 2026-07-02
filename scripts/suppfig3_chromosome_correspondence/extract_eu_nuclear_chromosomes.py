#!/usr/bin/env python3

import re
from pathlib import Path

input_fasta = Path(
    "09_suppfig_chr_correspondence/00_input/wcHirVerb1.fasta"
)

output_fasta = Path(
    "09_suppfig_chr_correspondence/v2_chr_only/"
    "00_input/wcHirVerb1_chr01_14.fasta"
)

mapping_tsv = Path(
    "09_suppfig_chr_correspondence/v2_chr_only/"
    "02_tables/wcHirVerb1_chromosome_accessions.tsv"
)


def read_fasta(path):
    header = None
    seq = []

    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)

                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())

        if header is not None:
            yield header, "".join(seq)


records = []

for header, sequence in read_fasta(input_fasta):
    match = re.search(
        r"chromosome:\s*([0-9]+)(?:\D|$)",
        header,
        flags=re.IGNORECASE,
    )

    if not match:
        continue

    chromosome_number = int(match.group(1))

    if not 1 <= chromosome_number <= 14:
        continue

    accession = header.split()[0]
    short_id = f"EU_chr{chromosome_number:02d}"

    records.append({
        "chromosome_number": chromosome_number,
        "short_id": short_id,
        "accession": accession,
        "original_header": header,
        "sequence": sequence,
        "length": len(sequence),
    })

records.sort(key=lambda row: row["chromosome_number"])

observed = [row["chromosome_number"] for row in records]
expected = list(range(1, 15))

if observed != expected:
    raise SystemExit(
        f"Expected chromosomes 1–14 but observed: {observed}"
    )

output_fasta.parent.mkdir(parents=True, exist_ok=True)
mapping_tsv.parent.mkdir(parents=True, exist_ok=True)

with output_fasta.open("w") as out:
    for row in records:
        out.write(
            f">{row['short_id']} "
            f"accession={row['accession']} "
            f"length={row['length']}\n"
        )

        sequence = row["sequence"]

        for index in range(0, len(sequence), 80):
            out.write(sequence[index:index + 80] + "\n")

with mapping_tsv.open("w") as out:
    out.write(
        "eu_chr\teu_chr_number\taccession\tlength_bp\toriginal_header\n"
    )

    for row in records:
        out.write(
            f"{row['short_id']}\t"
            f"{row['chromosome_number']}\t"
            f"{row['accession']}\t"
            f"{row['length']}\t"
            f"{row['original_header']}\n"
        )

print(f"Wrote: {output_fasta}")
print(f"Wrote: {mapping_tsv}")
print(f"Records: {len(records)}")
