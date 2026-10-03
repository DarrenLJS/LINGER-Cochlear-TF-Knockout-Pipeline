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

# provenance files: read-only from here on
chmod -R a-w "$GOLD_DIR"
log "done: ${GOLD_DIR}"
