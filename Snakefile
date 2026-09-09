# =============================================================================
# LINGER Cochlear GRN Pipeline
# Port of modules 1-2 (QC, consensus peaks) from the original standalone
# scripts, plus new Module 3 (integration & cell typing). Architecture
# mirrors the tRNA-seq Three-pass Alignment Snakemake Pipeline.
# =============================================================================
#
# Pipeline stages
# ---------------
#   01  QC & preprocessing        (per-sample RNA+ATAC load, filter, Scrublet)
#   02  Consensus peak set        (per-sample peaks -> bedtools merge -> re-quantify -> sync)
#   03  Integration & cell typing (Harmony batch correction, Leiden, UMAP, marker scoring)
#   -- modules 4-10 (LINGER init, chromatin priors, GRN inference, bulk TF
#      activity, perturbation, validation, benchmarking) not yet built --
#
# Conda environments
# ------------------
#   envs/linger_preproc.yaml — scanpy, anndata, scrublet, snapatac2, harmonypy,
#                               bedtools. Covers modules 1-3 entirely.
#   (envs/linger.yaml will be added in the module 4+ build for LingerGRN itself,
#    isolated the same way envs/mimseq.yaml isolates mim-tRNAseq in the
#    reference pipeline — LINGER's own dependency set is expected to conflict
#    with the scanpy/snapatac2 stack.)
#
# Usage
# -----
#   snakemake -n --use-conda --cores 8                     # dry run
#   snakemake --use-conda --cores 8 --rerun-incomplete --keep-going
#   snakemake --profile profiles/eddie                        # Eddie/SGE, see profiles/eddie/config.yaml
#
#   # NOTE: Do NOT manually conda activate an environment before running.
#   #       Snakemake handles per-rule activation via --use-conda.
# =============================================================================

import copy
import os
import pandas as pd

configfile: "config/config_eddie.yaml"

SCRATCH = config["scratch"]

# ---------------------------------------------------------------------------
# SGE resource helper — identical pattern to the tRNA-seq pipeline's
# Snakefile. Single place that assembles the `sge_extra` SGE flag string for
# a rule, sourced from config["resources"][rule_name] rather than a literal
# hardcoded per-rule. See config_eddie.yaml's `resources:` block docstring
# for why `sge_pe` must never be set directly on a rule instead of through
# this helper (EDDIE profile runs with --cores 1, which makes
# snakemake-executor-plugin-sge silently collapse sge_pe's slot count to 1
# regardless of a rule's own thread count).
# ---------------------------------------------------------------------------
def sge_extra(rule_name):
    r     = config["resources"][rule_name]
    slots = r.get("slots", 1)
    vmem  = r["vmem_mb"]
    pe    = f" -pe sharedmem {slots}" if slots > 1 else ""
    # h_rss added 2026-07-25 — Eddie's execd enforces a SEPARATE resident-
    # memory hard limit from h_vmem, and evidently has a low site-wide
    # default when h_rss isn't explicitly requested. linger_prep_pseudobulk
    # was killed (SIGKILL, qacct: "execd enforced h_rss limit") at 24GB
    # resident despite a 62.5GB h_vmem grant. Requesting h_rss explicitly,
    # matching h_vmem, closes this gap for every rule, not just that one.
    return f"-V{pe} -l h_vmem={vmem}M -l h_rss={vmem}M"


# ---------------------------------------------------------------------------
# Path-template resolution
#
# config_eddie.yaml writes sample paths as "{multiome_root}/GSE.../foo" so
# the root can be changed in exactly one place. Resolve those templates once
# here rather than repeating .format() calls in every rule/script.
# ---------------------------------------------------------------------------
_PATH_ROOTS = {
    "multiome_root": config["multiome_root"],
    "chrom_root":     config["chrom_root"],
    "scratch":        SCRATCH,
}


def _resolve(obj):
    if isinstance(obj, str):
        return obj.format(**_PATH_ROOTS)
    if isinstance(obj, dict):
        return {k: _resolve(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_resolve(v) for v in obj]
    return obj


SAMPLES_CFG               = _resolve(copy.deepcopy(config["samples"]))
EXTRA_BULK_RNA_CFG        = _resolve(copy.deepcopy(config.get("extra_bulk_rna_inputs", [])))
EXTRA_CHROMATIN_PRIOR_CFG = _resolve(copy.deepcopy(config.get("extra_chromatin_prior_inputs", [])))
INTEGRATION_CFG           = _resolve(copy.deepcopy(config["integration"]))
LINGER_CFG                = copy.deepcopy(config["linger"])  # not path-templated — grn_dir is absolute

# Added for Modules 7-10 integration (2026-09-08) ----------------------------
# model_construction_refs — Module 7's role="baseline" bucket (see
# 07_bulk_tf_activity.smk header for why these route to Module 7, not 6).
MODEL_CONSTRUCTION_REFS_CFG = _resolve(copy.deepcopy(config.get("model_construction_refs", [])))
# chromatin_prior_extra — Module 10's benchmarking ground truth. Different
# config key from EXTRA_CHROMATIN_PRIOR_CFG above (extra_chromatin_prior_inputs,
# which feeds nothing in Modules 1-3's DAG) — do not conflate the two.
EXTRA_CHROMATIN_PRIOR_CFG_FULL = _resolve(copy.deepcopy(config.get("chromatin_prior_extra", [])))
# held_out_inputs / negative_control_input — Module 9's frozen validation
# set. Kept in their own variables (not merged) per the repo's own
# 2026-07-24 decision to score the negative control with inverted
# pass/fail logic rather than folding it into the positive checks.
HELD_OUT_CFG          = _resolve(copy.deepcopy(config.get("held_out_inputs", [])))
NEGATIVE_CONTROL_CFG  = _resolve(copy.deepcopy(config.get("negative_control_input", {})))

SAMPLES_BY_ID = {s["sample_id"]: s for s in SAMPLES_CFG}
SAMPLES       = list(SAMPLES_BY_ID.keys())


def sample_entry(sample_id):
    return SAMPLES_BY_ID[sample_id]


def raw_inputs_for_sample(wildcards):
    """
    Raw input file(s) for qc_preprocess_sample, resolved per sample `type`
    so the DAG re-triggers Module 1 if the underlying raw data changes
    (not just if the script changes). Mirrors the three loaders in the
    original 01_qc_preprocessing.py (load_h5_combined / load_mtx_combined /
    load_mtx_separate).
    """
    s = sample_entry(wildcards.sample)
    if s["type"] == "h5_combined":
        return [s["h5_path"]]
    if s["type"] == "mtx_combined":
        p = s["combined_prefix"]
        return [p + "matrix.mtx.gz", p + "barcodes.tsv.gz", p + "features.tsv.gz"]
    if s["type"] == "mtx_separate":
        rp, ap = s["rna_prefix"], s["atac_prefix"]
        return [
            rp + "matrix.mtx.gz", rp + "barcodes.tsv.gz", rp + "features.tsv.gz",
            ap + "matrix.mtx.gz", ap + "barcodes.tsv.gz", s["atac_peaks_bed"],
        ]
    raise ValueError(f"Unknown sample type: {s['type']}")


# ---------------------------------------------------------------------------
# Include rule modules
# ---------------------------------------------------------------------------
include: "workflow/rules/01_qc_preprocessing.smk"
include: "workflow/rules/02_consensus_peaks.smk"
include: "workflow/rules/03_integration_celltyping.smk"
include: "workflow/rules/04_linger_init.smk"
include: "workflow/rules/06_grn_inference.smk"
include: "workflow/rules/07_bulk_tf_activity.smk"
include: "workflow/rules/08_perturbation.smk"
include: "workflow/rules/09_validation.smk"
include: "workflow/rules/10_grn_benchmarking.smk"
# 05_chromatin_priors.smk not built — Module 5 decided as Option A, folded
# into Module 10's benchmarking (see README "Design notes"), not a
# separate pipeline stage.
#
# Modules 6-10 are NOT added to `rule all` below, matching the existing
# module6_all convention — each depends on the previous stage's real
# output existing (not just its rule being defined), and Module 6 in
# particular depends on the manual Module 3 labeling step. Run each
# module's `_all` target explicitly once its prerequisites are done:
#   module6_all  (Module 6 — GRN inference)
#   module7_all  (Module 7 — bulk TF activity)
#   module8_all  (Module 8 — perturbation; RUN linger_perturbation_sanity_check
#                 first and read its log, see 08_perturbation.smk header)
#   module9_all  (Module 9 — validation)
#   benchmark_grn_edges  (Module 10 — GRN benchmarking; no dependency on 7-9)

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        # ── Module 1 ───────────────────────────────────────────────────────
        expand(f"{SCRATCH}/{{sample}}/rna_filtered.h5ad", sample=SAMPLES),
        expand(f"{SCRATCH}/{{sample}}/atac_filtered.h5ad", sample=SAMPLES),
        # ── Module 2 ───────────────────────────────────────────────────────
        f"{SCRATCH}/consensus_peaks.bed",
        expand(f"{SCRATCH}/{{sample}}/atac_consensus.h5ad", sample=SAMPLES),
        expand(f"{SCRATCH}/{{sample}}/rna_synced.h5ad", sample=SAMPLES),
        f"{SCRATCH}/qc_sanity_checks/sanity_summary.csv",
        # ── Module 3 ───────────────────────────────────────────────────────
        # NOTE: this does NOT include labeled.h5ad. That output depends on
        # cluster_annotation.tsv, which is a MANUAL deliverable — fill it in
        # yourself after reviewing the UMAP + marker scores below, then run
        # `snakemake --use-conda --cores 4 <scratch>/module3_integration/labeled.h5ad`
        # separately. See workflow/rules/03_integration_celltyping.smk.
        f"{SCRATCH}/module3_integration/integrated.h5ad",
        f"{SCRATCH}/module3_integration/cluster_marker_scores.csv",
        f"{SCRATCH}/module3_integration/cluster_umap.png",
        f"{SCRATCH}/module3_integration/cluster_template.tsv",
