# Hirudo verbana genome resource code

Custom scripts and documentation associated with:

**Chromosome-level assembly and annotation resource for the medicinal leech *Hirudo verbana***

## Resource and identifiers

**Hverb_1.0** is the submitted nuclear assembly name. The default nuclear reference comprises 14 chromosome-level records (Chr01–Chr14; 177,659,782 bp). A separate retained non-primary set contains 40 records (3,135,724 bp); the combined nuclear FASTA contains 54 records (180,795,506 bp).

The full nuclear annotation contains 17,789 genes and 30,131 transcripts across primary and retained non-primary records. Of these, 17,639 genes and 29,969 transcripts are on the 14 primary records. The separate mitochondrial genome model is 19,644 bp and contains 13 protein-coding genes, 22 tRNAs and two rRNAs.

Use `hirudo_augmented_release_v1.gff3` with `Hirudo_verbana_v8_primary_plus_retained.fasta` for the complete reusable nuclear annotation. For primary-only analyses, use `Hirudo_verbana_v8_primary_nuclear.fasta` and restrict annotation features to Chr01–Chr14. The `v8` strings in these existing filenames are retained for file traceability; they are not an alternative assembly name.

## Data access

- **Sequencing data:** DDBJ/DRA BioProject PRJDB42640; experiments DRX1045546–DRX1045557 and runs DRR1069955–DRR1069966. These sequencing records and the seven source BioSamples SAMD01936330–SAMD01936336 were released on 15 September 2026.
- **Supporting data:** associated Figshare record with reserved DOI **10.6084/m9.figshare.32880206**. The supporting-data release package contains sequence and annotation files, QC summaries, figure source data and manuscript-facing table exports.
- **Nuclear sequence submission:** derived BioSample SAMD02106527 and locus-tag prefix LGRP. Nuclear INSDC accession assignment is pending DDBJ Mass Submission System curation.
- **Mitochondrial sequence submission:** source BioSample SAMD01936330 and HiFi run DRR1069955. Hi-C run DRR1069956 was not used for mitochondrial recovery. Mitochondrial INSDC accession assignment is pending DDBJ curation.

See `docs/data_package.md` for data/code roles and reference-file selection.

## Repository structure

```text
scripts/
  final_figures/                    manuscript figure code bases
  fig2_genome_landscape/            Figure 2 source-track preparation
  suppfig1_chr04_local_curation/    Supplementary Figure 1 plotting
  suppfig2_chr04_read_support/      Chr04 read-support analysis
  suppfig3_chromosome_correspondence/ chromosome correspondence analysis
  assembly_qc/                      selected assembly-QC utilities
docs/
scripts_manifest.tsv
requirements.txt
LICENSE
```

The repository contains selected custom analysis, validation and plotting scripts. Large sequence/alignment inputs and deposited numerical source data are kept in DDBJ/DRA or the associated Figshare supporting-data record rather than duplicated here.

## Figure code and reproduction scope

The scripts under `scripts/final_figures/` reproduce the code-generated bases for Figure 2, Figure 3, Supplementary Figure 2, Supplementary Figure 3 and Supplementary Figure 4. The final manuscript artwork received subsequent typography/layout adjustments in Adobe Illustrator, so byte-identical reproduction of the Illustrator-edited PDFs is not claimed.

Figure 3 uses the validated raw, unbalanced 250-kb Cooler matrix, `log10(raw contact count + 1)`, deposited display limits and true cumulative sequence-coordinate bin widths. Figure 2 uses the display ranges stated in the manuscript caption.

Supplementary Figure 3 is reproducible from the deposited one-to-one NUCmer coordinate table,
`05_figure_source_data/Supplementary_Figure_S3/wcHirVerb1_vs_final_v8.current.1to1.coords.tsv`,
together with Supplementary Table S4. Reciprocal-best chromosome pairs are selected independently in both directions by summed one-to-one aligned sequence. Reported chromosome coverage is the union of merged one-to-one alignment intervals.

See `scripts/final_figures/README.md` for source-data inputs and example commands.

## Dependencies

`requirements.txt` lists the direct Python packages used by the selected scripts. The tested working versions recorded during release preparation are listed in `docs/tested_environment_20260924.md`.

## License and citation

Code in this repository is distributed under the MIT License; see `LICENSE`. Third-party software and data retain their own licenses and citation requirements. The code license does not define the license of the associated Figshare dataset.

Please cite the accompanying manuscript and the versioned supporting-data record when reusing this resource.
