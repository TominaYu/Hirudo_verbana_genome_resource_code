#!/usr/bin/env bash
set -euo pipefail

FA=${FA:?}
GTF=${GTF:?}
RMGFF=${RMGFF:?}

echo "=== chrom.sizes ==="
cut -f1 data/chrom.sizes

echo
echo "=== FASTA headers (first 30) ==="
grep '^>' "$FA" | sed 's/^>//' | cut -d' ' -f1 | head -n 30

echo
echo "=== GTF chromosome names (unique, first 50) ==="
grep -v '^#' "$GTF" | cut -f1 | sort -u | head -n 50

echo
echo "=== RepeatMasker GFF chromosome names (unique, first 50) ==="
grep -v '^#' "$RMGFF" | cut -f1 | sort -u | head -n 50

echo
echo "=== first 5 non-comment lines of GTF ==="
grep -v '^#' "$GTF" | head -n 5

echo
echo "=== first GTF gene line ==="
grep -v '^#' "$GTF" | awk '$3=="gene"{print; exit}'

echo
echo "=== first 10 non-comment lines of RepeatMasker GFF ==="
grep -v '^#' "$RMGFF" | head -n 10

echo
echo "=== GTF names not in chrom.sizes ==="
comm -23 \
  <(grep -v '^#' "$GTF" | cut -f1 | sort -u) \
  <(cut -f1 data/chrom.sizes | sort -u) || true

echo
echo "=== GFF names not in chrom.sizes ==="
comm -23 \
  <(grep -v '^#' "$RMGFF" | cut -f1 | sort -u) \
  <(cut -f1 data/chrom.sizes | sort -u) || true
