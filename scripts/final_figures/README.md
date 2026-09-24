# Final manuscript figure code bases

These scripts reproduce the code-generated bases for the manuscript figures listed below. The manuscript artwork was subsequently adjusted in Adobe Illustrator for typography and layout. Therefore, the scripts reproduce the numerical content, coordinates, display ranges and code-generated layout bases; they are not claimed to reproduce the Illustrator-edited PDF byte for byte.

The examples below assume the Figshare package root is stored in `$PKG`.

## Figure 2 — genome-wide feature landscape

Input:

- `05_figure_source_data/Figure2/Fig2_final_authoritative_source_data_100kb.tsv`

Display ranges:

- GC skew: -0.25 to 0.25
- GC fraction: 0.25 to 0.60
- gene density: 0 to 80 overlapping gene spans per 100-kb window
- LINE fraction: 0 to 0.35
- LTR fraction: 0 to 0.20
- interspersed TE fraction: 0 to 0.75

```bash
python scripts/final_figures/plot_Figure2.py \
  --input "$PKG/05_figure_source_data/Figure2/Fig2_final_authoritative_source_data_100kb.tsv" \
  --outdir figure_outputs
```

## Figure 3 — Hi-C contact map

Inputs:

- `05_figure_source_data/Figure3/final_v8_primary.250kb.raw.cool`
- `05_figure_source_data/Figure3/Fig3_final_v8_primary_250kb_raw_true_coordinates_display_parameters.tsv`

The script uses raw, unbalanced 250-kb contacts, applies `log10(raw contact count + 1)`, and draws terminal partial bins at their true cumulative sequence widths. If the display-parameter TSV includes the expected Cooler SHA256, the script checks it before rendering.

```bash
python scripts/final_figures/plot_Figure3.py \
  --cooler "$PKG/05_figure_source_data/Figure3/final_v8_primary.250kb.raw.cool" \
  --display-params "$PKG/05_figure_source_data/Figure3/Fig3_final_v8_primary_250kb_raw_true_coordinates_display_parameters.tsv" \
  --outdir figure_outputs
```

## Supplementary Figure 2 — alternative Chr04 models and HiFi support

Inputs:

- `05_figure_source_data/Supplementary_Figures_S1_S2/S2/SuppFig2B_strict_spanning_data.tsv`
- `05_figure_source_data/Supplementary_Figures_S1_S2/S2/SuppFig2C_GH1_gap_sweep_data.tsv`

```bash
python scripts/final_figures/plot_SuppFig2.py \
  --strict-spanning "$PKG/05_figure_source_data/Supplementary_Figures_S1_S2/S2/SuppFig2B_strict_spanning_data.tsv" \
  --gap-sweep "$PKG/05_figure_source_data/Supplementary_Figures_S1_S2/S2/SuppFig2C_GH1_gap_sweep_data.tsv" \
  --outdir figure_outputs
```

## Supplementary Figure 3 — chromosome correspondence

Public reproduction uses the deposited one-to-one NUCmer coordinates plus Supplementary Table S4.

Inputs:

- `05_figure_source_data/Supplementary_Figure_S3/wcHirVerb1_vs_final_v8.current.1to1.coords.tsv`
- `07_tables_submission/Supplementary_Table_S4_reciprocal_best_chromosome_correspondence_wcHirVerb1.tsv`

Regenerate chromosome-pair statistics, reciprocal-best pairs and visualization blocks:

```bash
python scripts/suppfig3_chromosome_correspondence/build_correspondence_tables.py \
  --coords "$PKG/05_figure_source_data/Supplementary_Figure_S3/wcHirVerb1_vs_final_v8.current.1to1.coords.tsv" \
  --s4 "$PKG/07_tables_submission/Supplementary_Table_S4_reciprocal_best_chromosome_correspondence_wcHirVerb1.tsv" \
  --outdir s3_tables
```

Render panels A and B from the regenerated tables:

```bash
python scripts/final_figures/plot_SuppFig3.py \
  --pair-stats s3_tables/all_chromosome_pair_statistics.tsv \
  --reciprocal-pairs s3_tables/reciprocal_best_chromosome_pairs.tsv \
  --visual-blocks s3_tables/visualized_one_to_one_blocks.tsv \
  --outdir figure_outputs
```

The Figshare package also contains compact publication-facing S3 source tables:

- `05_figure_source_data/Supplementary_Figure_S3/01_S3a_heatmap_matrix_Mb.tsv`
- `05_figure_source_data/Supplementary_Figure_S3/02_S3_reciprocal_pairs_publication.tsv`
- `05_figure_source_data/Supplementary_Figure_S3/03_S3b_visualized_block_summary.tsv`

Reciprocal-best pairs are selected independently in both directions using the largest summed one-to-one aligned length. Union coverage is calculated from merged alignment intervals; it is not summed aligned length divided by chromosome length.

The public reproduction workflow was validated against the RC2 package on 24 September 2026: the deposited coordinate table was byte-identical to the original analysis coordinate table; regenerated pair statistics, reciprocal-best pairs and visualization blocks matched the original processed tables; 14 unique reciprocal-best pairs were recovered; chromosome mapping, lengths, summed aligned bp, dominant orientation and union coverage matched Supplementary Table S4; and both S3 panels rendered successfully.

## Supplementary Figure 4 — mitochondrial genome map

Input:

- `05_figure_source_data/Supplementary_Figure_S4/01_final_mtDNA_annotation.tsv`

```bash
python scripts/final_figures/plot_SuppFig4.py \
  --input "$PKG/05_figure_source_data/Supplementary_Figure_S4/01_final_mtDNA_annotation.tsv" \
  --outdir figure_outputs
```

## Fonts and output formats

The shared helper selects Arial, Helvetica, Liberation Sans or DejaVu Sans in that order, depending on what is installed. PDF and SVG use editable text settings; PNG is exported at the requested DPI. Font substitution can cause small typographic differences without changing plotted data.
