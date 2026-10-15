#!/usr/bin/env bash
# =============================================================================
# download_gold_standard.sh  —  LINGER Pipeline: raw gold-standard tables
# =============================================================================
#
# Fetches the three public tables the frozen `directional_gold` table is built
# from (Track 4 / Round C):
#   GenAge human            genomics.senescence.info/genes/human_genes.zip
#   GenAge model organisms  genomics.senescence.info/genes/models_genes.zip
#   CellAge                 genomics.senescence.info/cells/cellAge.zip
#
# ADDED 2026-10-15 — two MSigDB MOUSE gene-set collections for the pathway-level metrics (Stage 6; read by
# compute_metrics.py, never used to build the gold table):
#   Hallmark  mh.all.v2023.2.Mm.symbols.gmt
#   GO:BP     m5.go.bp.v2023.2.Mm.symbols.gmt
#   from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Mm/
# They go into <DATA_ROOT>/gold_standard/msigdb/ with md5 + gene-set counts in the manifest. If the Broad server is
# not reachable from Eddie the script does NOT abort the three tables above: it prints the two URLs and the exact
# target path so the files can be downloaded on a laptop and copied over (scp), and exits non-zero at the very end
# so the missing files are not silently ignored. (A GMT that downloads but is HTML / empty / has < 40 sets stops
# the script immediately.)
#
# These are PROVENANCE files only. The pipeline never downloads anything at run
# time: a build script reads them once and writes a versioned table into the
# repository. That is why this script is separate from download_datasets.sh
# (which fetches raw GEO matrices, takes 12-48 h and makes them read-only).
#
# USAGE (login node, wget only, a few MB, about a minute):
#   bash setup_scripts/download_gold_standard.sh \
#        /exports/eddie/scratch/s2906787/cochlear_datasets
#
# Output (under <DATA_ROOT>/gold_standard/):
#   *.zip                 the downloads
#   genage_human/  genage_models/  cellage/   the extracted files
#   manifest.tsv          file, url, bytes, md5, downloaded_utc
# Re-running is safe: the three files are small, so each run re-downloads them and
# rewrites the manifest. The script stops (non-zero exit) if a download is missing, empty or
# not a valid zip. It finishes by printing each extracted table's line count,
# header and first rows, so the structure can be read without opening the files.
#
# Test hook: the three URLs can be overridden through the environment
# (URL_GENAGE_HUMAN, URL_GENAGE_MODELS, URL_CELLAGE).
# =============================================================================

set -euo pipefail

DATA_ROOT="${1:-/exports/eddie/scratch/s2906787/cochlear_datasets}"
GOLD_DIR="${DATA_ROOT}/gold_standard"

URL_GENAGE_HUMAN="${URL_GENAGE_HUMAN:-https://genomics.senescence.info/genes/human_genes.zip}"
URL_GENAGE_MODELS="${URL_GENAGE_MODELS:-https://genomics.senescence.info/genes/models_genes.zip}"
URL_CELLAGE="${URL_CELLAGE:-https://genomics.senescence.info/cells/cellAge.zip}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; exit 1; }

for tool in wget unzip md5sum; do
    command -v "$tool" >/dev/null 2>&1 || die "'$tool' not found on PATH"
done

mkdir -p "$GOLD_DIR"
# an earlier run leaves the directory read-only (see the end of this script)
chmod -R u+w "$GOLD_DIR" 2>/dev/null || true
cd "$GOLD_DIR"

# fetch <url> <zip file name> <extract dir>
fetch() {
    local url="$1" zip="$2" dest="$3"
    log "fetching ${url}"
    wget --tries=5 --timeout=60 --quiet -O "$zip.part" "$url" \
        || die "download failed: ${url}"
    [ -s "$zip.part" ] || die "empty download: ${url}"
    unzip -tq "$zip.part" >/dev/null 2>&1 || die "not a valid zip: ${url}"
    mv -f "$zip.part" "$zip"
    rm -rf "$dest"; mkdir -p "$dest"
    unzip -oq "$zip" -d "$dest" || die "unzip failed: ${zip}"
    [ -n "$(ls -A "$dest")" ] || die "zip extracted nothing: ${zip}"
}

fetch "$URL_GENAGE_HUMAN"  genage_human.zip   genage_human
fetch "$URL_GENAGE_MODELS" genage_models.zip  genage_models
fetch "$URL_CELLAGE"       cellAge.zip        cellage

# manifest
printf 'file\turl\tbytes\tmd5\tdownloaded_utc\n' > manifest.tsv
now="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
for pair in "genage_human.zip|$URL_GENAGE_HUMAN" "genage_models.zip|$URL_GENAGE_MODELS" "cellAge.zip|$URL_CELLAGE"; do
    f="${pair%%|*}"; u="${pair#*|}"
    printf '%s\t%s\t%s\t%s\t%s\n' "$f" "$u" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -d' ' -f1)" "$now" >> manifest.tsv
done
log "manifest:"; column -t -s $'\t' manifest.tsv || cat manifest.tsv

# structure preview: line count, header and first rows of every extracted table
log "extracted tables (paste this block back so the parser can be written against the real columns):"
find genage_human genage_models cellage -type f \( -iname '*.csv' -o -iname '*.tsv' -o -iname '*.txt' \) | sort | while read -r t; do
    echo "----- ${t}  ($(wc -l < "$t") lines)"
    head -n 4 "$t" | cut -c1-400
done

# ---- MSigDB mouse GMTs (Hallmark, GO:BP) -------------------------------------------------------------
MSIG_BASE="${MSIG_BASE:-https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Mm}"
MSIG_DIR="${GOLD_DIR}/msigdb"
mkdir -p "$MSIG_DIR"
MSIG_MISSING=0
# name | minimum number of gene sets (Hallmark has 50, GO:BP thousands)
for spec in "mh.all.v2023.2.Mm.symbols.gmt|40" "m5.go.bp.v2023.2.Mm.symbols.gmt|1000"; do
    f="${spec%%|*}"; min="${spec#*|}"
    if [ ! -s "${MSIG_DIR}/${f}" ]; then
        log "fetching ${MSIG_BASE}/${f}"
        if wget --tries=3 --timeout=60 --quiet -O "${MSIG_DIR}/${f}.part" "${MSIG_BASE}/${f}" && [ -s "${MSIG_DIR}/${f}.part" ]; then
            mv -f "${MSIG_DIR}/${f}.part" "${MSIG_DIR}/${f}"
        else
            rm -f "${MSIG_DIR}/${f}.part"
            echo "[$(date '+%F %T')] WARNING: could not download ${f}." >&2
            echo "    Download it on a machine with internet access:  ${MSIG_BASE}/${f}" >&2
            echo "    then copy it to Eddie:  scp ${f} <you>@eddie.ecdf.ed.ac.uk:${MSIG_DIR}/${f}" >&2
            MSIG_MISSING=1
            continue
        fi
    fi
    n=$(wc -l < "${MSIG_DIR}/${f}")
    head -c 200 "${MSIG_DIR}/${f}" | grep -qi '<html' && die "${f} is an HTML page, not a GMT"
    [ "$n" -ge "$min" ] || die "${f} has only ${n} gene sets (expected >= ${min}) — truncated or wrong file"
    awk -F'\t' 'NF < 3 {bad++} END {exit bad > 0}' "${MSIG_DIR}/${f}" || die "${f}: lines with < 3 tab-separated fields"
    printf '%s\t%s\t%s\t%s\t%s\n' "msigdb/${f}" "${MSIG_BASE}/${f}" "$(stat -c %s "${MSIG_DIR}/${f}")" "$(md5sum "${MSIG_DIR}/${f}" | cut -d' ' -f1)" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> manifest.tsv
    log "ok ${f}: ${n} gene sets"
done

# provenance files: read-only from here on
chmod -R a-w "$GOLD_DIR"
if [ "$MSIG_MISSING" -ne 0 ]; then
    die "MSigDB GMT file(s) missing (see the WARNING lines above) — copy them to ${MSIG_DIR} and re-run this script"
fi
log "done: ${GOLD_DIR}"
