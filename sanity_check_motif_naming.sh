#!/usr/bin/env bash
# sanity_check_motif_naming.sh — proves the consensus_peaks_named.bed fix
# actually makes HOMER emit real chr:start-end PositionIDs, on a 10-peak
# slice, before committing to the full 224,719-peak / 1,466-motif rerun.
#
# Run this directly on the login node — it's small and fast (seconds).
set -euo pipefail

SCRATCH=/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed
HOMER_BIN=/exports/eddie/scratch/s2906787/linger_pipeline/homer/bin
GRN_DIR=/exports/eddie/scratch/s2906787/linger_pipeline/provide_data/provide_data
GENOME=mm10
HOMER_GENOME_DIR="${HOMER_BIN%/*}/data/genomes/${GENOME}/"

TMPDIR=$(mktemp -d)
echo "Working in: $TMPDIR"
export PATH="$HOMER_BIN:$PATH"

echo ""
echo "=== Step A: take a 10-peak slice of the real consensus_peaks.bed ==="
head -10 "$SCRATCH/consensus_peaks.bed" > "$TMPDIR/slice.bed"
cat "$TMPDIR/slice.bed"

echo ""
echo "=== Step B (the fix): add a 4th name column = chr:start-end ==="
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1":"$2"-"$3}' "$TMPDIR/slice.bed" > "$TMPDIR/slice_named.bed"
cat "$TMPDIR/slice_named.bed"

echo ""
echo "=== Step C: homerTools extract on the NAMED slice ==="
homerTools extract "$TMPDIR/slice_named.bed" "$HOMER_GENOME_DIR" -fa > "$TMPDIR/seqs.fa"
echo "--- resulting FASTA headers (should now be chr:start-end, not default-N) ---"
grep "^>" "$TMPDIR/seqs.fa"

echo ""
echo "=== Step D: homer2 find on the extracted sequences ==="
homer2 find -i "$TMPDIR/seqs.fa" -m "$GRN_DIR/all_motif_rmdup_Mammal" -offset 0 > "$TMPDIR/motif_raw.txt"
echo "--- first few PositionID values (column 1) — the actual thing that was broken ---"
cut -f1 "$TMPDIR/motif_raw.txt" | sort -u | head -10

echo ""
echo "=== VERDICT ==="
if [ "$(cut -f1 "$TMPDIR/motif_raw.txt" | grep -c "^chr")" -gt 0 ]; then
    echo "PASS: PositionID values are real chr:start-end coordinates."
    echo "Safe to deploy the fix and run the full pipeline."
else
    echo "FAIL: PositionID values are still not coordinate-shaped — do NOT"
    echo "proceed to the full rerun. Paste this script's full output back for"
    echo "further diagnosis instead of guessing further."
fi

echo ""
echo "(Leaving $TMPDIR in place in case you want to inspect it further;"
echo " safe to 'rm -rf $TMPDIR' once you're satisfied.)"
