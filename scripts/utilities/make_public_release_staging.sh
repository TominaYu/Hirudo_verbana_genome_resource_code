#!/usr/bin/env bash
set -euo pipefail

# =========================
# 0) paths
# =========================
BASE="${HOME}/hirudo/project_primary_FASTA"
SRC_FINAL="${BASE}/07_finalize_v8_20260327_232659"

# 公表用 staging フォルダ名（必要に応じて変更）
DEST="${BASE}/90_public_release_staging_$(date +%Y%m%d)"

echo "[INFO] Creating staging folder: ${DEST}"

# =========================
# 1) create folder hierarchy
# =========================
mkdir -p "${DEST}"/{00_README,01_reference_assembly/{primary,retained_non_primary,combined},02_qc_summaries/{assembly_stats,busco,busco_comparison},03_tables_working/{main,supplementary},04_annotation_release/{gff3,fasta,tsv},05_raw_reads/{hifi,hic,rna_short,rna_kinnex},06_mtDNA/{assembly,validation},07_figures_source,99_internal_archive_optional}

# =========================
# 2) copy assembly FASTA files already identified
# =========================
cp -av \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/00_inputs/Hirudo_v8_nuclear_primary.fasta" \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/00_inputs/Hirudo_v8_nuclear_primary.fasta.fai" \
  "${DEST}/01_reference_assembly/primary/"

cp -av \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/00_inputs/Hirudo_v8_nuclear_provisional.fasta" \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/00_inputs/Hirudo_v8_nuclear_provisional.fasta.fai" \
  "${DEST}/01_reference_assembly/retained_non_primary/"

cp -av \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/00_inputs/Hirudo_v8_nuclear_primary_plus_unplaced.fasta" \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/00_inputs/Hirudo_v8_nuclear_primary_plus_unplaced.fasta.fai" \
  "${DEST}/01_reference_assembly/combined/"

# =========================
# 3) copy assembly/QC summary files used for Table 2/3/4
# =========================
cp -av \
  "${SRC_FINAL}/02_qc_reconstruct_compat/03_summary.txt" \
  "${SRC_FINAL}/02_qc_reconstruct_compat/05_group_summary.tsv" \
  "${SRC_FINAL}/02_qc_reconstruct_compat/01_primary_lengths.tsv" \
  "${SRC_FINAL}/02_qc_reconstruct_compat/02_provisional_lengths.tsv" \
  "${SRC_FINAL}/02_qc_reconstruct_compat/06_reconstruction_map.tsv" \
  "${DEST}/02_qc_summaries/assembly_stats/"

cp -av \
  "${SRC_FINAL}/04_tables/Hirudo_v8_unplaced_component_table.tsv" \
  "${SRC_FINAL}/04_tables/Hirudo_v8_unplaced_component_table.txt" \
  "${DEST}/02_qc_summaries/assembly_stats/"

cp -av \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/02_stats/primary_summary.tsv" \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/02_stats/unplaced_summary.tsv" \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/02_stats/combined_summary.tsv" \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/02_stats/primary_lengths.tsv" \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/02_stats/unplaced_lengths.tsv" \
  "${SRC_FINAL}/08_genome_qc_20260328_000200/02_stats/combined_lengths.tsv" \
  "${DEST}/02_qc_summaries/assembly_stats/"

# =========================
# 4) copy BUSCO summaries used for Table 7
# =========================
cp -av \
  "${SRC_FINAL}/09_busco6_setup_20260328_002021/04_busco_runs/busco6_primary/short_summary.specific.metazoa_odb12.busco6_primary.txt" \
  "${SRC_FINAL}/09_busco6_setup_20260328_002021/04_busco_runs/busco6_primary_lophotrochozoa/short_summary.specific.lophotrochozoa_odb12.busco6_primary_lophotrochozoa.txt" \
  "${SRC_FINAL}/09_busco6_setup_20260328_002021/04_busco_runs/busco6_primary_eukaryota_odb10/short_summary.specific.eukaryota_odb10.busco6_primary_eukaryota_odb10.txt" \
  "${SRC_FINAL}/09_busco6_setup_20260328_002021/04_busco_runs/busco6_combined_lophotrochozoa/short_summary.specific.lophotrochozoa_odb12.busco6_combined_lophotrochozoa.txt" \
  "${DEST}/02_qc_summaries/busco/"

# JSON も残したい場合
cp -av \
  "${SRC_FINAL}/09_busco6_setup_20260328_002021/04_busco_runs/busco6_primary/short_summary.specific.metazoa_odb12.busco6_primary.json" \
  "${SRC_FINAL}/09_busco6_setup_20260328_002021/04_busco_runs/busco6_primary_lophotrochozoa/short_summary.specific.lophotrochozoa_odb12.busco6_primary_lophotrochozoa.json" \
  "${SRC_FINAL}/09_busco6_setup_20260328_002021/04_busco_runs/busco6_primary_eukaryota_odb10/short_summary.specific.eukaryota_odb10.busco6_primary_eukaryota_odb10.json" \
  "${SRC_FINAL}/09_busco6_setup_20260328_002021/04_busco_runs/busco6_combined_lophotrochozoa/short_summary.specific.lophotrochozoa_odb12.busco6_combined_lophotrochozoa.json" \
  "${DEST}/02_qc_summaries/busco/"

# BUSCO比較表
if compgen -G "${SRC_FINAL}/09_busco6_setup_20260328_002021/05_compare_tables/*" > /dev/null; then
  cp -av \
    "${SRC_FINAL}/09_busco6_setup_20260328_002021/05_compare_tables/"* \
    "${DEST}/02_qc_summaries/busco_comparison/"
fi

# =========================
# 5) create placeholder working files for tables
# =========================
cat > "${DEST}/03_tables_working/README_tables.txt" <<'EOF'
This directory is for manuscript-facing table preparation.

Suggested mapping:
- Table 2 / 3 / 4: derive from 02_qc_summaries/assembly_stats
- Table 7: derive from 02_qc_summaries/busco and busco_comparison
- Table 5 / 6 / Supplementary Table S1: annotation-side tables to be added later
EOF

# 空のテンプレート置き場
touch "${DEST}/03_tables_working/main/Table2_nuclear_genome_assembly_statistics.tsv"
touch "${DEST}/03_tables_working/main/Table3_chromosome_level_record_summary.tsv"
touch "${DEST}/03_tables_working/main/Table4_retained_non_primary_sequence_summary.tsv"
touch "${DEST}/03_tables_working/main/Table7_technical_validation_summary.tsv"
touch "${DEST}/03_tables_working/supplementary/SupplementaryTableS2_retained_non_primary_records.tsv"

# =========================
# 6) create placeholders for files not yet located in this chat
# =========================
cat > "${DEST}/05_raw_reads/README_raw_reads.txt" <<'EOF'
Raw read directories created, but exact source file paths have not yet been identified in this chat.

To fill later:
- hifi/
- hic/
- rna_short/
- rna_kinnex/

Recommended action:
1. Identify exact raw FASTQ/BAM source paths.
2. Copy or symlink them here.
3. Record original source locations in a manifest.
EOF

cat > "${DEST}/04_annotation_release/README_annotation_release.txt" <<'EOF'
Place annotation release files here once finalized:
- slim GFF3
- public GFF3
- protein FASTA
- CDS FASTA
- transcript FASTA
- functional annotation TSV
- long-read support summary TSV
EOF

cat > "${DEST}/06_mtDNA/README_mtDNA.txt" <<'EOF'
Place finalized mitochondrial-genome files here once selected:
- mitochondrial FASTA
- annotation map
- validation summary tables
EOF

# =========================
# 7) create top-level README and manifest skeleton
# =========================
cat > "${DEST}/00_README/README_public_release_staging.txt" <<EOF
Public-release staging folder for the Hirudo verbana genome manuscript.

Created: $(date)
Source project root: ${BASE}
Primary source finalize folder: ${SRC_FINAL}

Included now:
- final primary / provisional / combined FASTA files
- assembly summary TSV/TXT files used for Tables 2/3/4
- BUSCO summary TXT/JSON files used for Table 7
- unplaced component table

Not yet added:
- raw HiFi / Hi-C / RNA-seq / Kinnex read files
- finalized annotation release files
- finalized mtDNA release files
- accession / DOI / checksum manifests
EOF

# file manifest skeleton
find "${DEST}" -type f | sort > "${DEST}/00_README/file_manifest_initial.txt"

echo "[INFO] Done."
echo "[INFO] Staging folder created at: ${DEST}"
echo "[INFO] Review the following first:"
echo "  ${DEST}/00_README/README_public_release_staging.txt"
echo "  ${DEST}/00_README/file_manifest_initial.txt"
