# Associated data record and file selection

## Repository roles

- **DDBJ/DRA:** PacBio HiFi genomic reads, Hi-C reads, five short-read RNA-seq datasets and five tissue-specific Kinnex FLNC BAM datasets under BioProject PRJDB42640.
- **DDBJ nucleotide sequence archive:** Hverb_1.0 nuclear assembly/annotation and the mitochondrial sequence submission. Accession assignment is pending curator completion at the time of this code-release preparation.
- **Figshare:** associated supporting-data record, reserved DOI 10.6084/m9.figshare.32880206, containing sequence/annotation release files, QC summaries, figure source data and table exports.
- **GitHub:** selected custom analysis, validation and plotting code.

## Reference-file selection

For chromosome-focused analyses, use the 14-record `Hirudo_verbana_v8_primary_nuclear.fasta`.

For analyses using the complete 17,789-gene, 30,131-transcript reusable annotation, use `hirudo_augmented_release_v1.gff3` with `Hirudo_verbana_v8_primary_plus_retained.fasta`. The combined FASTA contains the 14 primary records and 40 retained non-primary records. Analyses restricted to the primary FASTA must likewise restrict GFF3 features to Chr01–Chr14.

The retained non-primary sequences are provided for reuse but are not validated biological haplotypes or structural variants. The mitochondrial sequence and annotation are separate resources.

The DDBJ representation contains fewer mRNA/CDS feature pairs because alternative transcripts with identical CDS locations within a locus are collapsed for submission compatibility. Use the augmented GFF3 rather than the reduced DDBJ representation for multi-isoform analyses.

## Figures and source data

The Figshare package contains numerical source data for Figure 2, Figure 3 and Supplementary Figures 1–4. The GitHub plotting scripts accept deposited source files as explicit inputs. Final manuscript PDFs may include subsequent Adobe Illustrator typography/layout edits; these do not replace the deposited numerical source tables.

### Supplementary Figure 3

The primary deposited analysis input is:

- `05_figure_source_data/Supplementary_Figure_S3/wcHirVerb1_vs_final_v8.current.1to1.coords.tsv`

Supplementary Table S4 provides the wcHirVerb1 chromosome/accession/length mapping needed to regenerate the chromosome-pair summaries:

- `07_tables_submission/Supplementary_Table_S4_reciprocal_best_chromosome_correspondence_wcHirVerb1.tsv`

The package also contains compact publication-facing S3 source tables:

- `01_S3a_heatmap_matrix_Mb.tsv`
- `02_S3_reciprocal_pairs_publication.tsv`
- `03_S3b_visualized_block_summary.tsv`

Chromosome-pair assignment is based on reciprocal-best summed one-to-one aligned sequence, selected independently in both directions. Reported chromosome coverage is the percentage of each chromosome covered by the union of merged one-to-one alignment intervals, not summed aligned length divided by chromosome length.
