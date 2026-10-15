#!/usr/bin/env bash
# Download NCBI mouse gene_info and write a GeneID<TAB>Symbol table used by the bulk loader's
# `entrez_symbol_map` option (needed for datasets keyed only by Entrez ID, e.g. GSE171921).
# Usage: bash setup_scripts/download_gene_annotation.sh /exports/eddie/scratch/s2906787/cochlear_datasets
# Run on the Eddie login node (needs internet).
set -euo pipefail
ROOT="${1:?usage: $0 <DATA_ROOT>}"
OUT="$ROOT/gene_annotation"; mkdir -p "$OUT"
URL="https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"
curl -fL --retry 3 -o "$OUT/Mus_musculus.gene_info.gz" "$URL"
# gene_info columns: #tax_id GeneID Symbol ...  (tax 10090 only)
zcat "$OUT/Mus_musculus.gene_info.gz" | awk -F'\t' 'BEGIN{OFS="\t"; print "GeneID","Symbol"} !/^#/ && $1==10090 {print $2,$3}' > "$OUT/mouse_entrez_symbol.tsv"
N=$(($(wc -l < "$OUT/mouse_entrez_symbol.tsv")-1))
echo "wrote $OUT/mouse_entrez_symbol.tsv ($N genes)"
[ "$N" -gt 20000 ] || { echo "ERROR: too few genes mapped ($N)"; exit 1; }
( cd "$OUT" && md5sum Mus_musculus.gene_info.gz mouse_entrez_symbol.tsv > MANIFEST.md5 && date -u +"downloaded %FT%TZ from $URL" >> MANIFEST.md5 )
