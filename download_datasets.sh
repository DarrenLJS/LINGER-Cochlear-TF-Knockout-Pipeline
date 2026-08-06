#!/usr/bin/env bash
# =============================================================================
# download_datasets.sh  —  LINGER Pipeline: Full GEO Dataset Download
# Tailored for University of Edinburgh Eddie HPC, scratch space
# =============================================================================
#
# USAGE (in tmux session 2 — fully independent from setup):
#   ssh s2906787@eddie.ecdf.ed.ac.uk
#   tmux new-session -s download
#   bash /exports/eddie/scratch/s2906787/download_datasets.sh \
#        /exports/eddie/scratch/s2906787/cochlear_datasets
#
# Runs on a LOGIN NODE — no qlogin/module needed (wget only).
# Re-running is safe — wget --continue resumes incomplete files.
#
# Watch progress:
#   tail -f /exports/eddie/scratch/s2906787/cochlear_datasets/download.log
#   grep "file(s)" /exports/eddie/scratch/s2906787/cochlear_datasets/download.log
#
# Estimated time  : 12-48 hours
# Estimated space : 150-500 GB
# =============================================================================

set -uo pipefail

DATA_ROOT="${1:-/exports/eddie/scratch/s2906787/cochlear_datasets}"
LOG="${DATA_ROOT}/download.log"
NTHREADS="${NTHREADS:-8}"
mkdir -p "$DATA_ROOT"

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN:  $*" | tee -a "$LOG"; }
sep()  { log "============================================================"; }

# Pre-flight disk check
AVAIL_KB=$(df -k "/exports/eddie/scratch/s2906787" 2>/dev/null | awk 'NR==2 {print $4}')
AVAIL_GB=$(( AVAIL_KB / 1024 / 1024 ))

sep
log "LINGER Dataset Download"
log "  DATA_ROOT = ${DATA_ROOT}"
log "  SCRATCH   = ${AVAIL_GB} GB available"
(( AVAIL_GB < 600 )) && warn "  Less than 600 GB free — monitor with: df -h /exports/eddie/scratch/s2906787"
log "  Follow with: tail -f ${LOG}"
log "  Count done: grep 'file(s)' ${LOG}"
sep

# =============================================================================
# HELPER: download_geo_suppl  <ACCESSION>  <DEST_DIR>
#
# FIX vs previous version:
#   - "--append-output=${LOG}" was writing ALL wget verbose output (including
#     the per-50KB progress bar lines) directly to the log file, bypassing
#     the grep filter pipe. This caused a 500,000+ line log for a single file.
#   - Now uses "--quiet" (suppresses progress bars entirely) with stderr
#     redirected to the log. Only genuine errors appear in the log.
#   - The pipe and grep filter are removed — no longer needed.
#
# GEO FTP path:  ftp.ncbi.nlm.nih.gov/geo/series/GSE{N÷1000}nnn/GSE{N}/suppl/
# =============================================================================
download_geo_suppl() {
    local accession="$1"
    local dest_dir="$2"
    mkdir -p "$dest_dir"

    local digits="${accession#GSE}"
    local prefix_num=$(( 10#$digits / 1000 ))
    local prefix="GSE${prefix_num}nnn"
    local suppl_url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${accession}/suppl/"

    log "  [${accession}] downloading → ${dest_dir}"

    wget \
        --quiet \
        --recursive \
        --no-parent \
        --no-directories \
        --continue \
        --tries=5 \
        --timeout=180 \
        --wait=2 \
        --random-wait \
        --reject "index.html,robots.txt" \
        --directory-prefix="${dest_dir}" \
        "${suppl_url}" \
        2>>"$LOG" \
    || warn "  [${accession}] wget returned non-zero — partial download or empty suppl dir"

    local n size
    n=$(find "$dest_dir" -maxdepth 1 -type f 2>/dev/null | wc -l)
    size=$(du -sh "$dest_dir" 2>/dev/null | awk '{print $1}')
    log "  [${accession}] done — ${n} file(s)  ${size}"
}

# =============================================================================
# HELPER: download_sra  <GSE_ACCESSION>  <DEST_DIR>
#
# Fallback for datasets with no processed matrices on GEO.
# Requires linger_tools (prefetch/fasterq-dump) and LINGER (pysradb).
# Only call this after confirming matrix files are absent from the GEO suppl.
#
# To use:
#   module load anaconda
#   source $(conda info --base)/etc/profile.d/conda.sh
#   conda activate linger_tools
#   source /exports/eddie/scratch/s2906787/download_datasets.sh   # defines the function
#   download_sra GSE157398 /exports/eddie/scratch/s2906787/cochlear_datasets/01_multiome/GSE157398_sra
# =============================================================================
download_sra() {
    local accession="$1"
    local dest_dir="$2"
    mkdir -p "$dest_dir"

    if ! command -v prefetch &>/dev/null; then
        warn "  prefetch not found — activate linger_tools first"
        return 1
    fi

    log "  [${accession}] Fetching SRR list via pysradb..."
    local srr_file="${dest_dir}/srr_list.txt"

    # pysradb is in the LINGER env; conda run avoids needing a second activate
    conda run -n LINGER --no-capture-output python3 << PYCODE > "$srr_file" 2>>"$LOG"
from pysradb.sraweb import SRAweb
import sys
db = SRAweb()
df = db.gse_to_srp('${accession}')
if df is None or df.empty:
    sys.exit(0)
srp = df['study_accession'].iloc[0]
runs = db.srp_to_srr(srp)
if runs is not None and not runs.empty:
    for r in runs['run_accession'].unique():
        print(r)
PYCODE

    if [ ! -s "$srr_file" ]; then
        warn "  [${accession}] No SRR IDs found — check SRA Run Selector manually"
        return 1
    fi

    while IFS= read -r srr; do
        [[ -z "$srr" ]] && continue
        if [ -f "${dest_dir}/${srr}_1.fastq.gz" ] || [ -f "${dest_dir}/${srr}.fastq.gz" ]; then
            log "    [skip] ${srr} already downloaded"
            continue
        fi
        log "    Prefetching ${srr}..."
        prefetch --output-directory "$dest_dir" "$srr" >> "$LOG" 2>&1 \
        && fasterq-dump \
               --outdir "$dest_dir" --threads "$NTHREADS" --split-files \
               "${dest_dir}/${srr}" >> "$LOG" 2>&1 \
        && pigz -p "$NTHREADS" "${dest_dir}/${srr}"*.fastq 2>>"$LOG" \
        || warn "    ${srr} failed — check ${LOG}"
    done < "$srr_file"
}

# =============================================================================
# SECTION 1 — PRIMARY MULTIOME DATA                        [CRITICAL for LINGER]
# After download, verify each folder contains:
#   RNA  : barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
#   ATAC : peaks.bed  fragments.tsv.gz  OR  peak_matrix.*
# If absent, the dataset only deposited FASTQs — use download_sra fallback.
# =============================================================================
log ""
log "=== SECTION 1: Primary Multiome Inputs ==="
MULTI="${DATA_ROOT}/01_multiome"

download_geo_suppl "GSE157398"  "${MULTI}/GSE157398"    # hair-cell scMultiome (P2, dev)
download_geo_suppl "GSE331301"  "${MULTI}/GSE331301"    # single-cell multiome
download_geo_suppl "GSE182202"  "${MULTI}/GSE182202"    # A/GA/GAP reprogramming (P8/P15)
download_geo_suppl "GSE224563"  "${MULTI}/GSE224563"    # mature SC barrier + GRN (P70)

# =============================================================================
# SECTION 2 — CHROMATIN PRIOR DATA
# Expected: .hic/.cool (Hi-C), .bed/.narrowPeak (CUT&RUN/ATAC)
# =============================================================================
log ""
log "=== SECTION 2: Chromatin Prior Datasets ==="
CHROM="${DATA_ROOT}/02_chromatin_priors"

download_geo_suppl "GSE305205"  "${CHROM}/GSE305205_hic"           # Hi-C 3D genome
download_geo_suppl "GSE150386"  "${CHROM}/GSE150386_cutrun_atoh1"  # ATOH1/POU4F3 CUT&RUN ★
download_geo_suppl "GSE150391"  "${CHROM}/GSE150391_cutrun_pou4f3" # POU4F3 pioneer ★
download_geo_suppl "GSE181307"  "${CHROM}/GSE181307_cutrun"        # additional CUT&RUN
download_geo_suppl "GSE181310"  "${CHROM}/GSE181310_atac"          # ATAC regulatory prior
download_geo_suppl "GSE288375"  "${CHROM}/GSE288375_ihc_ohc_rna"   # IHC/OHC subtype RNA
download_geo_suppl "GSE288376"  "${CHROM}/GSE288376_ihc_ohc_atac"  # IHC/OHC subtype ATAC ★

# =============================================================================
# SECTION 3 — MODEL CONSTRUCTION DATA
# =============================================================================
log ""
log "=== SECTION 3: Model Construction Datasets ==="
MODEL="${DATA_ROOT}/03_model_construction"

download_geo_suppl "GSE135703"  "${MODEL}/GSE135703_adult_SC"
download_geo_suppl "GSE283534"  "${MODEL}/GSE283534_mature_hc"
download_geo_suppl "GSE111349"  "${MODEL}/GSE111349_sorted_ihc_ohc"
download_geo_suppl "GSE202920"  "${MODEL}/GSE202920_P14_P28"
download_geo_suppl "GSE136196"  "${MODEL}/GSE136196_stria"
download_geo_suppl "GSE168041"  "${MODEL}/GSE168041_sgn_lateralwall"

# =============================================================================
# SECTION 4 — BULK RNA-SEQ TRAINING DATA  (~34 datasets)
# Run up to 5 in parallel to save time without overloading GEO FTP.
# =============================================================================
log ""
log "=== SECTION 4: Bulk RNA-seq Training Datasets ==="
BULK="${DATA_ROOT}/04_bulk_rnaseq"

BULK_ACCESSIONS=(
    GSE111603 GSE111604 GSE113719 GSE116702 GSE120462
    GSE129533 GSE139145 GSE149254 GSE149257 GSE163798
    GSE193158 GSE198167 GSE203228 GSE213784 GSE239361
    GSE240605 GSE242321 GSE242669 GSE256069 GSE262205
    GSE266157 GSE268504 GSE271610 GSE275243 GSE294736
    GSE296324 GSE299064 GSE300428 GSE312253 GSE312476
    GSE329565 GSE83599  GSE85519  GSE90821
)

log "  Downloading ${#BULK_ACCESSIONS[@]} datasets (max 5 parallel)..."
JOBS=0
for acc in "${BULK_ACCESSIONS[@]}"; do
    download_geo_suppl "$acc" "${BULK}/${acc}" &
    JOBS=$(( JOBS + 1 ))
    if (( JOBS >= 5 )); then
        wait
        JOBS=0
    fi
done
wait
log "  All bulk RNA-seq downloads complete"

# =============================================================================
# SECTION 5 — TUNING VALIDATION DATA
# =============================================================================
log ""
log "=== SECTION 5: Tuning Validation Datasets ==="
TUNE="${DATA_ROOT}/05_tuning_validation"

download_geo_suppl "GSE150000"  "${TUNE}/GSE150000_SC_maturation"    # SC maturation + enhancer decommissioning (RNA-seq + ATAC)
download_geo_suppl "GSE150002"  "${TUNE}/GSE150002_notch_dapt"       # Notch/DAPT regeneration validation (RNA-seq + ATAC)
download_geo_suppl "GSE182202"  "${TUNE}/GSE182202_gap_reprogramming" # A/GA/GAP reprogramming P8/P15 (also GRN prior — data already in 01_multiome)
download_geo_suppl "GSE224563"  "${TUNE}/GSE224563_mature_SC_barrier" # Mature SC barrier P70 (also GRN prior — data already in 01_multiome)
download_geo_suppl "GSE234926"  "${TUNE}/GSE234926_P90_sensory"
download_geo_suppl "GSE165662"  "${TUNE}/GSE165662_adult_stria"
download_geo_suppl "GSE189798"  "${TUNE}/GSE189798"
download_geo_suppl "GSE266947"  "${TUNE}/GSE266947_small_rna"
download_geo_suppl "GSE289514"  "${TUNE}/GSE289514_serpine2"
download_geo_suppl "GSE297020"  "${TUNE}/GSE297020_aging"
download_geo_suppl "GSE305163"  "${TUNE}/GSE305163_otoferlin"
download_geo_suppl "GSE310246"  "${TUNE}/GSE310246_ebf1"
download_geo_suppl "GSE316951"  "${TUNE}/GSE316951"

# =============================================================================
# SECTION 6 — HELD-OUT VALIDATION DATA         !! FROZEN — DO NOT OPEN !!
# chmod 555 applied immediately after download.
# =============================================================================
log ""
log "=== SECTION 6: Held-Out Validation Datasets [FROZEN] ==="
HELD="${DATA_ROOT}/06_held_out"

download_geo_suppl "GSE224627"  "${HELD}/GSE224627_GAP_reprogramming"  # ★ primary GAP test
download_geo_suppl "GSE233559"  "${HELD}/GSE233559_Tbx2_OHC_IHC"       # Tbx2 conversion
download_geo_suppl "GSE274279"  "${HELD}/GSE274279_aging_snRNAseq"      # ★ primary aging
download_geo_suppl "GSE153882"  "${HELD}/GSE153882_hc_aging"
download_geo_suppl "GSE154833"  "${HELD}/GSE154833_stria_aging"
download_geo_suppl "GSE196870"  "${HELD}/GSE196870_SNHL"
download_geo_suppl "GSE273939"  "${HELD}/GSE273939_TNF_organoids"

log "  Locking ${HELD} (chmod 555)..."
chmod -R 555 "${HELD}" 2>/dev/null \
    && log "  [LOCKED] 06_held_out — do not unlock until final evaluation" \
    || warn "  chmod failed — run manually: chmod -R 555 ${HELD}"

# =============================================================================
# SECTION 7 — NEGATIVE CONTROL
# =============================================================================
log ""
log "=== SECTION 7: Negative Control ==="
NEG="${DATA_ROOT}/07_negative_control"
download_geo_suppl "GSE281207"  "${NEG}/GSE281207_Tmie_control"

# =============================================================================
# SECTION 8 — REFERENCE-ONLY
# =============================================================================
log ""
log "=== SECTION 8: Reference-Only Datasets ==="
REF="${DATA_ROOT}/08_reference_only"

download_geo_suppl "GSE152551"  "${REF}/GSE152551_stria_spindle"
download_geo_suppl "GSE181057"  "${REF}/GSE181057_fibrocyte"
download_geo_suppl "GSE312224"  "${REF}/GSE312224_stria_endothelial"

# =============================================================================
# SECTION 9 — AUXILIARY LOW-WEIGHT
# =============================================================================
log ""
log "=== SECTION 9: Auxiliary Low-Weight Datasets ==="
AUX="${DATA_ROOT}/09_auxiliary"

download_geo_suppl "GSE122732"  "${AUX}/GSE122732"
download_geo_suppl "GSE127683"  "${AUX}/GSE127683"
download_geo_suppl "GSE209791"  "${AUX}/GSE209791"
download_geo_suppl "GSE240187"  "${AUX}/GSE240187"
download_geo_suppl "GSE283708"  "${AUX}/GSE283708_vestibular"

# =============================================================================
# SUMMARY
# =============================================================================
sep
log "Download complete."
log ""
log "Disk usage by section:"
du -sh "${DATA_ROOT}"/0*/  2>/dev/null | tee -a "$LOG"
log ""
log "Total:"
du -sh "${DATA_ROOT}" 2>/dev/null | tee -a "$LOG"
sep

log ""
log "POST-DOWNLOAD CHECKLIST"
log ""
log "1. Check multiome folders have matrix files (not just FASTQs):"
log "   for d in ${DATA_ROOT}/01_multiome/*/; do"
log "     echo \"=== \$(basename \$d) ===\""
log "     ls \"\$d\" | grep -E 'barcodes|features|matrix|peaks|fragments|h5'"
log "   done"
log ""
log "2. Check chromatin priors:"
log "   ls ${DATA_ROOT}/02_chromatin_priors/GSE305205_hic/    # expect .hic or .cool"
log "   ls ${DATA_ROOT}/02_chromatin_priors/GSE150386*/       # expect .bed/.narrowPeak"
log ""
log "3. Held-out lock:"
log "   ls -ld ${DATA_ROOT}/06_held_out/   # should show dr-xr-xr-x (555)"
log ""
log "4. SRA fallback (if matrix files missing from a multiome dataset):"
log "   module load anaconda"
log "   source \$(conda info --base)/etc/profile.d/conda.sh"
log "   conda activate linger_tools"
log "   source ${DATA_ROOT}/../download_datasets.sh   # loads download_sra function"
log "   download_sra GSE157398 ${DATA_ROOT}/01_multiome/GSE157398_sra"
log ""
log "Full log: ${LOG}"
log "Log size: $(du -sh "${LOG}" | cut -f1)  (kept small by --quiet flag)"
sep
