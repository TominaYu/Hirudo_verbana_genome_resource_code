# Hirudo verbana genome resource code

This repository contains selected scripts used to generate, validate, and document tables and figures associated with the manuscript:

**Chromosome-level genome assembly and annotation of the medicinal leech *Hirudo verbana***

The primary data package will be deposited separately in a DOI-bearing repository such as Figshare or Zenodo. This GitHub repository contains scripts and documentation only. It does not contain raw sequencing reads, FASTA/GFF3 data files, BAM/SAM files, Hi-C `.hic` files, or other large intermediate files.

## Associated data package

Current local data package candidate:

`Hirudo_verbana_SciData_v07_public_upload_candidate_with_SuppFig5_and_Fig3HiC_20260629_214727.tar.gz`

MD5:

`ad0b623bf5fbb24f3bc31f8e9c160c4f`

The data package contains the primary nuclear assembly, retained non-primary sequences, annotation files, repeat annotation files, QC summaries, figure source data, and submission-ready tables. Raw sequencing reads are being deposited separately in DDBJ/DRA.

After repository deposition, replace this section with the final DOI.

## Repository structure

```text
scripts/
  fig2_genome_landscape/
  suppfig1_chr04_local_curation/
  suppfig2_chr04_read_support/
  suppfig3_chromosome_correspondence/
  assembly_qc/
  utilities/

scripts_manifest.tsv
requirements.txt
LICENSE
Script categories
fig2_genome_landscape/: scripts for genome feature tracks used in the genome landscape figure.
suppfig1_chr04_local_curation/: scripts used to prepare and plot local chromosome curation panels.
suppfig2_chr04_read_support/: scripts used to summarize and plot read-backed support for alternative local backbone models.
suppfig3_chromosome_correspondence/: scripts used to summarize chromosome correspondence with a published H. verbana assembly.
assembly_qc/: scripts used for basic assembly statistics, BUSCO parsing, telomere checks, and validation tables.
utilities/: helper scripts used during public release staging.
Data inputs

The scripts assume that the DOI data package has been downloaded and unpacked. Paths may need to be adjusted depending on the local working directory.

Large inputs such as FASTA, GFF3, BAM, .hic, and raw read files are intentionally not stored in this GitHub repository.

Scope

These scripts are provided as provenance and reproducibility aids. They are not intended to reproduce the complete genome assembly and annotation pipeline from raw reads. The main analysis workflow used standard tools described in the manuscript, while the scripts here primarily support figure generation, table preparation, and QC summarization.

Requirements

Python dependencies are listed in requirements.txt.

pip install -r requirements.txt
License

Code in this repository is released under the MIT License.
