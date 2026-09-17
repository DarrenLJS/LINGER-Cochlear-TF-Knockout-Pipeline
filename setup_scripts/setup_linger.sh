#!/usr/bin/env bash
# =============================================================================
# setup_linger.sh  —  LINGER GRN Pipeline: Environment & Tool Setup
# Tailored for University of Edinburgh Eddie HPC, scratch space
# =============================================================================
#
# USAGE (in tmux session 1):
#   ssh s2906787@eddie.ecdf.ed.ac.uk
#   tmux new-session -s setup
#   qlogin -l h_vmem=8G -l h_rt=04:00:00        # interactive node
#   module load anaconda
#   bash /exports/eddie/scratch/s2906787/setup_linger.sh \
#        /exports/eddie/scratch/s2906787/linger_pipeline
#
# Re-running is safe — all steps are idempotent:
#   existing conda envs are kept; packages are upgraded/verified
#   existing reference files are skipped; provide_data is skipped if present
#
# Optional flags (set before running):
#   INSTALL_R=1 bash setup_linger.sh ...    # also install R + Seurat (~60 min)
#
# Watch progress:
#   tail -f /exports/eddie/scratch/s2906787/linger_pipeline/setup.log
# =============================================================================

set -uo pipefail     # -u: unset vars are errors   -o pipefail: pipe errors caught
                     # NO -e: conda activate/deactivate can return non-zero

WORK_DIR="${1:-/exports/eddie/scratch/s2906787/linger_pipeline}"
LOG="${WORK_DIR}/setup.log"
INSTALL_R="${INSTALL_R:-0}"
mkdir -p "$WORK_DIR"

# ---- Logging ----------------------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN:  $*" | tee -a "$LOG"; }
die()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" | tee -a "$LOG"; exit 1; }
sep()  { log "------------------------------------------------------------"; }

log "=================================================================="
log "  LINGER setup"
log "  WORK_DIR  = ${WORK_DIR}"
log "  USER      = $(whoami)   HOST = $(hostname)"
log "  INSTALL_R = ${INSTALL_R}"
log "=================================================================="
sep

# =============================================================================
# STEP 1 — VERIFY CONDA IS AVAILABLE
# User runs "module load anaconda" before calling this script.
# The subshell inherits the updated PATH so conda is reachable.
# We source conda.sh here so that conda activate/deactivate work throughout.
# =============================================================================
log "Step 1 — Verifying conda is in PATH"

if ! command -v conda &>/dev/null; then
    die "conda not found in PATH.
    Did you run:  module load anaconda
    before calling this script?"
fi
log "  conda found: $(conda --version)"

CONDA_BASE=$(conda info --base 2>/dev/null) \
    || die "conda info --base failed"
source "${CONDA_BASE}/etc/profile.d/conda.sh" \
    || die "Failed to source ${CONDA_BASE}/etc/profile.d/conda.sh"
log "  Sourced: ${CONDA_BASE}/etc/profile.d/conda.sh"

# =============================================================================
# STEP 2 — DISK SPACE CHECK
# =============================================================================
sep
log "Step 2 — Checking scratch disk space"

SCRATCH_ROOT="/exports/eddie/scratch/s2906787"
AVAIL_KB=$(df -k "$SCRATCH_ROOT" 2>/dev/null | awk 'NR==2 {print $4}')
AVAIL_GB=$(( AVAIL_KB / 1024 / 1024 ))
log "  Available on scratch: ${AVAIL_GB} GB"
(( AVAIL_GB < 50 )) && warn "  Less than 50 GB free — monitor with: df -h ${SCRATCH_ROOT}"

# =============================================================================
# STEP 3 — CONFIGURE CONDA ENV/PKG PATHS
# =============================================================================
sep
log "Step 3 — Configuring conda env/pkg paths"

SCRATCH_ENVS="${WORK_DIR}/conda_envs"
SCRATCH_PKGS="${WORK_DIR}/conda_pkgs"
mkdir -p "$SCRATCH_ENVS" "$SCRATCH_PKGS"

if grep -q "envs_dirs" "${HOME}/.condarc" 2>/dev/null; then
    EXISTING_ENV_DIR=$(grep -A3 "envs_dirs:" "${HOME}/.condarc" 2>/dev/null \
        | grep -E "^\s+-\s+" | head -1 | sed 's/^[[:space:]]*-[[:space:]]*//')
    log "  ~/.condarc envs_dirs: ${EXISTING_ENV_DIR}"
    if echo "${EXISTING_ENV_DIR}" | grep -qE "^${HOME}|^~"; then
        die "envs_dirs points to home directory — will fill quota.
    Edit ~/.condarc and change envs_dirs to:  - ${SCRATCH_ENVS}"
    fi
else
    cat >> "${HOME}/.condarc" << CONDARC

# Added by setup_linger.sh ($(date '+%Y-%m-%d'))
envs_dirs:
  - ${SCRATCH_ENVS}
pkgs_dirs:
  - ${SCRATCH_PKGS}
CONDARC
    log "  ~/.condarc updated — envs → ${SCRATCH_ENVS}"
fi

conda config --add channels defaults    2>>"$LOG" || true
conda config --add channels bioconda    2>>"$LOG" || true
conda config --add channels conda-forge 2>>"$LOG" || true
conda config --set channel_priority flexible 2>>"$LOG" || true
log "  Channels: conda-forge > bioconda > defaults (flexible)"

# =============================================================================
# STEP 4 — LINGER CONDA ENVIRONMENT
#
# FIX vs previous version:
#   - "conda run -n LINGER pip install" was silently failing because conda run
#     can pick up the wrong pip, especially for first installs into a fresh env.
#   - Now uses explicit "conda activate LINGER → pip install → conda deactivate"
#     which is reliable since we sourced conda.sh in Step 1.
#   - pip install runs UNCONDITIONALLY (not gated on env existence) so re-runs
#     catch any previously failed package installs.
#   - Failure of "import linger" is now a die(), not a warn().
# =============================================================================
sep
log "Step 4 — LINGER environment  (Python 3.10.0 + LingerGRN + scanpy stack)"
log "  Slowest step — typically 20-40 min"

if conda env list 2>/dev/null | grep -qE "^LINGER[[:space:]]"; then
    log "  LINGER env already exists — skipping conda create"
else
    conda create -y -n LINGER python==3.10.0 >> "$LOG" 2>&1 \
        && log "  [OK] Environment created" \
        || die "conda create LINGER failed — see ${LOG}"
fi

# bedtools + pybedtools 0.12.0 via conda (pre-compiled binary, no source build)
#
# WHY conda, and WHY --no-deps for LingerGRN (Phase 1b below):
#   LingerGRN==1.105 pins pybedtools==0.10.0. That version only exists on PyPI
#   as a source tarball — no pre-built wheel. Every pip-based approach fails:
#     - default:            isolated build env has no setuptools → error
#     - --no-build-isolation: pip subprocess still can't find setuptools
#   conda's bioconda provides pybedtools 0.12.0 as a pre-compiled binary.
#   We install it here, then install LingerGRN with --no-deps (Phase 1b) to
#   skip pip's pybedtools==0.10.0 version check entirely.
#   pybedtools 0.12.0 is fully API-compatible with 0.10.0 for all functions
#   LingerGRN calls — the strict pin is metadata-only.
log "  Installing bedtools + pybedtools 0.12.0 into LINGER (via conda)..."
conda install -y -n LINGER -c conda-forge -c bioconda bedtools pybedtools >> "$LOG" 2>&1 \
    && log "  [OK] bedtools + pybedtools (conda 0.12.0)" \
    || warn "  bedtools/pybedtools install had errors — check ${LOG}"

# pip install — three phases:
#   1a. All of LingerGRN's strict transitive deps EXCEPT pybedtools
#       (explicit list derived from pip resolver output across debug runs)
#   1b. LingerGRN itself with --no-deps, bypassing the pybedtools version check
#   2.  Minimal runtime extras that don't conflict with LingerGRN's pins
#
# Phase 1a — torch (~2 GB) + all other pinned deps. All have pre-built wheels.
LINGER_CORE="${WORK_DIR}/requirements_linger_core.txt"
cat > "$LINGER_CORE" << 'REQS'
torch
scipy==1.11.3
numpy==1.24.3
pandas==2.0.3
shap==0.42.0
scikit-learn==1.3.0
joblib==1.3.2
matplotlib==3.8.0
seaborn==0.13.0
statsmodels==0.14.1
umap-learn
scanpy==1.9.5
anndata==0.9.2
REQS

log "  Phase 1a/2 — installing LingerGRN deps (torch ~2 GB, ~15 min)..."
conda activate LINGER
pip install -r "$LINGER_CORE" 2>&1 | tee -a "$LOG"
PIP_STATUS=${PIPESTATUS[0]}
conda deactivate
[[ $PIP_STATUS -eq 0 ]] \
    && log "  [OK] LingerGRN deps installed" \
    || die "LingerGRN deps install failed — see ${LOG}"

# Phase 1b — LingerGRN package itself, skipping pybedtools==0.10.0 check
log "  Phase 1b/2 — installing LingerGRN==1.105 --no-deps..."
conda activate LINGER
pip install --no-deps LingerGRN==1.105 2>&1 | tee -a "$LOG"
PIP_STATUS=${PIPESTATUS[0]}
conda deactivate
[[ $PIP_STATUS -eq 0 ]] \
    && log "  [OK] LingerGRN==1.105 installed" \
    || die "LingerGRN install failed — see ${LOG}"

# Phase 2 — MINIMAL set: only packages LINGER actually calls at runtime.
#
# WHY SO FEW PACKAGES HERE:
#   Phase 2 previously included muon, episcanpy, scanorama, scrublet, pydeseq2,
#   harmonypy, pysradb, GEOparse. Those packages require anndata >= 0.10 and
#   scipy >= 1.13 — incompatible with LingerGRN's hard pins (anndata==0.9.2,
#   scipy==1.11.3). pip resolved the conflict by UPGRADING both, which then
#   caused 'import linger' to fail with ModuleNotFoundError because LingerGRN's
#   __init__.py hit API incompatibilities in the upgraded packages.
#
#   Fix: only install what LingerGRN actually imports at runtime here. All the
#   preprocessing packages (episcanpy, scanorama, muon, scrublet, etc.) are
#   already in linger_preproc, which has no version pinning constraints.
#
# rpy2 — LingerGRN uses R internally and requires rpy2==3.5.16. Install via
#   conda (not pip) so conda handles the R dependency automatically. pip would
#   require R to already be present and may need a source build.
log "  Installing rpy2 into LINGER (via conda — handles R dependency)..."
conda install -y -n LINGER -c conda-forge rpy2 >> "$LOG" 2>&1 \
    && log "  [OK] rpy2 installed" \
    || warn "  rpy2 install had errors — check ${LOG}"

LINGER_EXTRA="${WORK_DIR}/requirements_linger_extra.txt"
cat > "$LINGER_EXTRA" << 'REQS'
leidenalg
python-igraph
gdown
REQS

log "  Phase 2/2 — installing LINGER runtime extras (leidenalg, igraph, gdown)..."
conda activate LINGER
pip install -r "$LINGER_EXTRA" 2>&1 | tee -a "$LOG"
PIP_STATUS=${PIPESTATUS[0]}
conda deactivate
[[ $PIP_STATUS -eq 0 ]] \
    && log "  [OK] Extra packages installed" \
    || die "Extra package install failed for LINGER env — see ${LOG}"

# Verify import — use conda run (more reliable than activate+python after
# multiple activate/deactivate cycles in a non-interactive script)
log "  Verifying: import LingerGRN..."
conda run -n LINGER python -c \
    "import LingerGRN; print('  [OK] LingerGRN imported successfully')" \
    2>&1 | tee -a "$LOG"
IMPORT_STATUS=${PIPESTATUS[0]}
[[ $IMPORT_STATUS -eq 0 ]] \
    && log "  [OK] import LingerGRN verified" \
    || die "'import LingerGRN' failed — see output above for the actual error"

log "  LINGER environment ready"

# =============================================================================
# STEP 5 — PREPROCESSING ENVIRONMENT (linger_preproc)
# Same fix: activate → pip install → deactivate
# =============================================================================
sep
log "Step 5 — linger_preproc environment  (scanpy, snapatac2, QC tools)"

if conda env list 2>/dev/null | grep -qE "^linger_preproc[[:space:]]"; then
    log "  linger_preproc env already exists — skipping conda create"
else
    conda create -y -n linger_preproc python=3.10 >> "$LOG" 2>&1 \
        && log "  [OK] Environment created" \
        || die "conda create linger_preproc failed — see ${LOG}"
fi

PREPROC_REQS="${WORK_DIR}/requirements_preproc.txt"
cat > "$PREPROC_REQS" << 'REQS'
scanpy
anndata
muon
snapatac2
scrublet
scanorama
harmonypy
leidenalg
python-igraph
umap-learn
pysradb
GEOparse
pydeseq2
biopython
pandas
numpy
scipy
matplotlib
seaborn
REQS

log "  Installing pip packages (activate → pip → deactivate)..."
conda activate linger_preproc
pip install -r "$PREPROC_REQS" 2>&1 | tee -a "$LOG"
PIP_STATUS=${PIPESTATUS[0]}
conda deactivate
[[ $PIP_STATUS -eq 0 ]] \
    && log "  [OK] pip packages installed" \
    || die "pip install failed for linger_preproc — see ${LOG}"

log "  Verifying: import scanpy..."
conda activate linger_preproc
python -c "import scanpy; print('  scanpy version:', scanpy.__version__)" \
    2>&1 | tee -a "$LOG"
IMPORT_STATUS=${PIPESTATUS[0]}
conda deactivate
[[ $IMPORT_STATUS -eq 0 ]] \
    && log "  [OK] import scanpy verified" \
    || warn "  scanpy import failed — check ${LOG}"

# Optional R + Seurat
if [[ "${INSTALL_R}" == "1" ]]; then
    log "  INSTALL_R=1 — installing R 4.3 + Seurat (~60 min)..."
    conda install -y -n linger_preproc -c conda-forge \
        r-base=4.3 r-seurat r-harmony r-ggplot2 >> "$LOG" 2>&1 \
        && log "  [OK] R + Seurat" \
        || warn "  R/Seurat install failed — retry individually if needed"
else
    log "  Skipping R/Seurat (re-run with INSTALL_R=1 to include)"
fi

log "  linger_preproc environment ready"

# =============================================================================
# STEP 6 — BIOINFORMATICS TOOLS ENVIRONMENT (linger_tools)
# Dedicated env for alignment/peak-calling/SRA — never modifies base.
# =============================================================================
sep
log "Step 6 — linger_tools environment  (STAR, samtools, sra-tools, MACS2...)"

if conda env list 2>/dev/null | grep -qE "^linger_tools[[:space:]]"; then
    log "  linger_tools env already exists — skipping conda create"
else
    conda create -y -n linger_tools python=3.10 >> "$LOG" 2>&1 \
        && log "  [OK] Environment created" \
        || die "conda create linger_tools failed — see ${LOG}"
fi

log "  Installing bioinformatics tools..."
{
conda install -y -n linger_tools \
    -c conda-forge -c bioconda \
    sra-tools samtools star bowtie2 subread picard deeptools macs2 htslib pigz
} >> "$LOG" 2>&1 \
    && log "  [OK] Bioinformatics tools" \
    || warn "  Some tools had errors — check ${LOG}"

# =============================================================================
# STEP 7 — CELL RANGER ARC  (manual — requires 10x Genomics login)
# =============================================================================
sep
log "Step 7 — Cell Ranger ARC  [MANUAL — action on your laptop]"
CELLRANGER_TOOLS="${WORK_DIR}/tools"
mkdir -p "$CELLRANGER_TOOLS"
cat << MANUAL | tee -a "$LOG"

  Cell Ranger ARC needs a free 10x Genomics account. Cannot auto-download.

  ON YOUR LAPTOP:
    1. https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads
    2. Download: cellranger-arc-2.0.2.tar.gz

  TRANSFER TO EDDIE:
    scp cellranger-arc-2.0.2.tar.gz \\
        s2906787@eddie.ecdf.ed.ac.uk:${CELLRANGER_TOOLS}/

  ON EDDIE (after scp):
    tar -zxvf ${CELLRANGER_TOOLS}/cellranger-arc-2.0.2.tar.gz -C ${CELLRANGER_TOOLS}/
    echo 'export PATH="${CELLRANGER_TOOLS}/cellranger-arc-2.0.2:\$PATH"' >> ~/.bashrc
    source ~/.bashrc && cellranger-arc --version

  Only needed if GEO deposits lack processed count matrices.

MANUAL

# =============================================================================
# STEP 8 — LINGER provide_data (TSS locations + TF motifs for mm10)
#
# This replaces the previous "pretrained weights" step. LINGER does not use
# pretrained model checkpoints — it trains its own neural network on your data.
# What it DOES need is a reference data package (provide_data/) containing:
#   - TSS location files for mm10
#   - TF motif information (pre-formatted for LINGER's GRN construction)
# This becomes the GRNdir argument in all LINGER function calls.
#
# Downloaded from Google Drive via gdown (installed in LINGER env, Step 4).
# File ID: 1Dog5JTS_SNIoa5aohgZmOWXrTUuAKHXV
# =============================================================================
sep
log "Step 8 — Downloading LINGER provide_data (TSS + TF motifs for mm10)"

PROVIDE_DATA_DIR="${WORK_DIR}/provide_data"
mkdir -p "$PROVIDE_DATA_DIR"
PROVIDE_TARBALL="${PROVIDE_DATA_DIR}/provide_data.tar.gz"

if [ -d "${PROVIDE_DATA_DIR}/provide_data" ] && \
   [ "$(ls -A "${PROVIDE_DATA_DIR}/provide_data" 2>/dev/null)" ]; then
    log "  provide_data/ already present and non-empty — skipping download"
else
    log "  Downloading from Google Drive via gdown (in LINGER env)..."
    conda activate LINGER
    python -m gdown 1Dog5JTS_SNIoa5aohgZmOWXrTUuAKHXV \
        -O "$PROVIDE_TARBALL" \
        2>&1 | tee -a "$LOG"
    GDOWN_STATUS=${PIPESTATUS[0]}
    conda deactivate

    if [[ $GDOWN_STATUS -ne 0 ]]; then
        warn "  gdown failed — manual fallback:"
        cat << FALLBACK | tee -a "$LOG"

  Manual provide_data download (run on Eddie):
    conda activate LINGER
    cd ${PROVIDE_DATA_DIR}
    wget --load-cookies /tmp/cookies.txt \\
      "https://drive.usercontent.google.com/download?export=download&confirm=\$(wget \\
      --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies \\
      --no-check-certificate \\
      'https://drive.usercontent.google.com/download?id=1Dog5JTS_SNIoa5aohgZmOWXrTUuAKHXV' \\
      -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1Dog5JTS_SNIoa5aohgZmOWXrTUuAKHXV" \\
      -O provide_data.tar.gz && rm -rf /tmp/cookies.txt
    tar -xzf provide_data.tar.gz
    conda deactivate

FALLBACK
    else
        log "  Extracting provide_data.tar.gz..."
        tar -xzf "$PROVIDE_TARBALL" -C "$PROVIDE_DATA_DIR" \
            && log "  [OK] provide_data extracted to ${PROVIDE_DATA_DIR}" \
            || warn "  tar extraction failed — try manually: tar -xzf ${PROVIDE_TARBALL} -C ${PROVIDE_DATA_DIR}"
        # List top-level contents so user can confirm what landed
        log "  Contents of ${PROVIDE_DATA_DIR}/:"
        ls "${PROVIDE_DATA_DIR}/" 2>&1 | tee -a "$LOG"
    fi
fi

log "  GRNdir for LINGER calls = ${PROVIDE_DATA_DIR}/provide_data"
log "  (use this path as the 'GRNdir' argument in LINGER_tr.training() etc.)"

# =============================================================================
# STEP 9 — HOMER (motif scanning, required by LINGER pipeline)
#
# LINGER's scNN.md tutorial calls findMotifsGenome.pl from HOMER to scan
# ATAC peaks for TF motif enrichment. HOMER is not on Eddie as a module
# and is not conda-installable — it uses its own Perl installer.
#
# Two sub-steps:
#   9a. Install HOMER core (Perl scripts + databases, ~500 MB)
#   9b. Install mm10 genome package (~3-4 GB) — submitted as a qsub job
#       because it is too large to run on a login node.
# =============================================================================
sep
log "Step 9 — Installing HOMER (motif scanner)"

HOMER_DIR="${WORK_DIR}/homer"
mkdir -p "$HOMER_DIR"

# 9a. Core HOMER install
if [ -f "${HOMER_DIR}/bin/homer" ]; then
    log "  HOMER core already installed — skipping"
else
    log "  Downloading HOMER installer..."
    wget -q http://homer.ucsd.edu/homer/configureHomer.pl \
        -O "${HOMER_DIR}/configureHomer.pl" \
        || die "Failed to download configureHomer.pl — check network from this node"

    log "  Installing HOMER core (~500 MB, 5-10 min)..."
    perl "${HOMER_DIR}/configureHomer.pl" -install \
        2>&1 | tee -a "$LOG" \
        && log "  [OK] HOMER core installed" \
        || die "HOMER core install failed — see ${LOG}"
fi

# Add HOMER to PATH for this session
export PATH="${HOMER_DIR}/bin:${PATH}"

# Add to ~/.bashrc if not already there (persists across sessions)
if ! grep -qF "${HOMER_DIR}/bin" "${HOME}/.bashrc" 2>/dev/null; then
    echo "export PATH=\"${HOMER_DIR}/bin:\$PATH\"" >> "${HOME}/.bashrc"
    log "  Added ${HOMER_DIR}/bin to ~/.bashrc"
fi

# 9b. mm10 genome package — write as a qsub job (too large for login node)
HOMER_MM10_JOB="${WORK_DIR}/homer_mm10_install.sh"
cat > "$HOMER_MM10_JOB" << HOMERJOB
#!/bin/bash
#$ -N homer_mm10
#$ -l h_vmem=8G
#$ -l h_rt=01:30:00
#$ -cwd
#$ -o ${WORK_DIR}/homer_mm10.o
#$ -e ${WORK_DIR}/homer_mm10.e

export PATH="${HOMER_DIR}/bin:\${PATH}"
perl ${HOMER_DIR}/configureHomer.pl -install mm10
echo "HOMER mm10 genome package install complete"
HOMERJOB

log "  HOMER mm10 genome package job written: ${HOMER_MM10_JOB}"
log "  Submit it now with:"
log "    qsub ${HOMER_MM10_JOB}"
log "  (requires ~3-4 GB download, ~1 hour)"
log "  Check status: qstat -j homer_mm10"

# Verify HOMER core is reachable
if command -v findMotifsGenome.pl &>/dev/null; then
    log "  [OK] findMotifsGenome.pl is in PATH"
else
    warn "  findMotifsGenome.pl not found — source ~/.bashrc or log out and back in"
fi

# =============================================================================
# STEP 10 — LINGER REPOSITORY
# =============================================================================
sep
log "Step 10 — Cloning LINGER repository"

LINGER_REPO="${WORK_DIR}/LINGER_repo"
if [ -d "${LINGER_REPO}/.git" ]; then
    log "  Already cloned — pulling latest"
    git -C "$LINGER_REPO" pull --ff-only >> "$LOG" 2>&1 || warn "  git pull failed"
else
    git clone https://github.com/Durenlab/LINGER.git "$LINGER_REPO" >> "$LOG" 2>&1 \
        && log "  [OK] Cloned to ${LINGER_REPO}" \
        || die "git clone failed — check outbound network from this node"
fi

# =============================================================================
# STEP 11 — MOUSE mm10 REFERENCE FILES
# =============================================================================
sep
log "Step 11 — Downloading mm10 (GRCm38) reference files"

REF_DIR="${WORK_DIR}/references/mm10"
mkdir -p "$REF_DIR"

if [ ! -f "${REF_DIR}/mm10.fa.gz" ]; then
    log "  Downloading mm10 genome FASTA (~850 MB)..."
    wget -q --show-progress -c \
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz" \
        -O "${REF_DIR}/mm10.fa.gz" \
        && log "  [OK] mm10.fa.gz" || warn "  mm10 FASTA download failed"
else
    log "  mm10.fa.gz already present ($(du -sh "${REF_DIR}/mm10.fa.gz" | cut -f1))"
fi

if [ ! -f "${REF_DIR}/GRCm38.102.gtf.gz" ]; then
    log "  Downloading Ensembl GRCm38 r102 GTF..."
    wget -q --show-progress -c \
        "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz" \
        -O "${REF_DIR}/GRCm38.102.gtf.gz" \
        && log "  [OK] GRCm38.102.gtf.gz" || warn "  GTF download failed"
else
    log "  GRCm38.102.gtf.gz already present"
fi

if [ ! -f "${REF_DIR}/mm10.chrom.sizes" ]; then
    wget -q -c \
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes" \
        -O "${REF_DIR}/mm10.chrom.sizes" \
        && log "  [OK] mm10.chrom.sizes" || warn "  chrom.sizes download failed"
fi

JASPAR_DIR="${WORK_DIR}/references/jaspar"
mkdir -p "$JASPAR_DIR"
JASPAR_ZIP="${JASPAR_DIR}/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.zip"
if [ ! -f "$JASPAR_ZIP" ]; then
    log "  Downloading JASPAR 2022 vertebrate motifs..."
    wget -q --show-progress -c \
        "https://jaspar.elixir.no/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.zip" \
        -O "$JASPAR_ZIP" \
    && unzip -q -o "$JASPAR_ZIP" -d "${JASPAR_DIR}/" >> "$LOG" 2>&1 \
    && log "  [OK] JASPAR motifs extracted" \
    || warn "  JASPAR download failed"
else
    log "  JASPAR motifs already present"
fi

# =============================================================================
# STEP 12 — STAR INDEX JOB TEMPLATE
# =============================================================================
sep
log "Step 12 — Writing STAR index job template"

STAR_IDX="${REF_DIR}/star_index"
STAR_JOB="${WORK_DIR}/star_index_job.sh"
cat > "$STAR_JOB" << STARJOB
#!/bin/bash
#$ -N star_index
#$ -l h_vmem=42G
#$ -l h_rt=02:00:00
#$ -pe sharedmem 8
#$ -cwd
#$ -o ${WORK_DIR}/star_index.o
#$ -e ${WORK_DIR}/star_index.e

module load anaconda
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate linger_tools

mkdir -p ${STAR_IDX}
gunzip -k ${REF_DIR}/mm10.fa.gz
gunzip -k ${REF_DIR}/GRCm38.102.gtf.gz

STAR \\
  --runMode genomeGenerate \\
  --genomeDir ${STAR_IDX} \\
  --genomeFastaFiles ${REF_DIR}/mm10.fa \\
  --sjdbGTFfile ${REF_DIR}/GRCm38.102.gtf \\
  --runThreadN \${NSLOTS:-8} \\
  --genomeSAindexNbases 14

echo "STAR index complete: ${STAR_IDX}"
STARJOB
log "  qsub ${STAR_JOB}   (only if you need to align bulk RNA-seq FASTQs)"

# =============================================================================
# FINAL VERIFICATION
# =============================================================================
sep
log "Final verification..."

# Environment existence
for env in LINGER linger_preproc linger_tools; do
    if conda env list 2>/dev/null | grep -qE "^${env}[[:space:]]"; then
        ENV_PATH=$(conda env list 2>/dev/null | grep -E "^${env}[[:space:]]" | awk '{print $NF}')
        log "  [OK] ${env}  →  ${ENV_PATH}"
    else
        warn "  [MISSING] ${env}"
    fi
done

# Critical package imports
log ""
log "  Package import checks:"

conda activate LINGER
python -c "
import LingerGRN, scanpy, anndata
print('  LingerGRN: OK')
print('  scanpy   :', scanpy.__version__)
print('  anndata  :', anndata.__version__)
" 2>&1 | tee -a "$LOG" \
    && log "  [OK] LINGER env imports" \
    || warn "  LINGER env: one or more imports failed — check output above"
conda deactivate

conda activate linger_preproc
python -c "
import scanpy, scrublet, snapatac2
print('  scanpy   :', scanpy.__version__)
print('  scrublet : OK')
print('  snapatac2:', snapatac2.__version__)
" 2>&1 | tee -a "$LOG" \
    && log "  [OK] linger_preproc env imports" \
    || warn "  linger_preproc: one or more imports failed"
conda deactivate

conda activate linger_tools
STAR --version 2>&1 | head -1 | tee -a "$LOG"
samtools --version 2>&1 | head -1 | tee -a "$LOG"
conda deactivate
log "  [OK] linger_tools binaries"

# provide_data check
if [ -d "${PROVIDE_DATA_DIR}/provide_data" ]; then
    NFILES=$(find "${PROVIDE_DATA_DIR}/provide_data" -type f | wc -l)
    log "  [OK] provide_data  →  ${PROVIDE_DATA_DIR}/provide_data  (${NFILES} files)"
else
    warn "  provide_data/ directory missing — check Step 8 output above"
fi

sep
log ""
log "SUMMARY — paths you will use in your LINGER scripts:"
log "  GRNdir       = ${PROVIDE_DATA_DIR}/provide_data"
log "  LINGER_repo  = ${LINGER_REPO}"
log "  mm10 FASTA   = ${REF_DIR}/mm10.fa.gz"
log "  mm10 GTF     = ${REF_DIR}/GRCm38.102.gtf.gz"
log "  HOMER        = ${HOMER_DIR}/bin"
log ""
log "PENDING MANUAL ACTIONS:"
log "  1. Submit HOMER mm10 genome job:   qsub ${HOMER_MM10_JOB}"
log "  2. Download Cell Ranger ARC from 10x website (only if needed)"
log "  3. Submit STAR index job:          qsub ${STAR_JOB}  (only if needed)"
log ""
log "Full log: ${LOG}"
sep
