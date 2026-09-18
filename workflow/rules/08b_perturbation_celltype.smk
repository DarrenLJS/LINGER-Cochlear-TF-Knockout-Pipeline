# =============================================================================
# workflow/rules/08b_perturbation_celltype.smk
# Module 8b — cell-type-resolved in silico perturbation
#
# ADDED 2026-09-16. Module 8 (08_perturbation.smk) forward-passes ONE
# population-pooled pseudobulk through the shared, population-trained
# {chr}_net.pt files. Nothing about the net or its (call-time-recomputed)
# normalization requires that — this module swaps in a cell-type-restricted
# pseudobulk instead, per cell type, reusing the exact same nets and the
# exact same persisted TF row order Module 8's sanity check already
# validated (median Spearman rho=0.797 against the population Target). See
# linger_perturbation.py's CELL-TYPE-RESOLVED MODE docstring section and
# prep_pseudobulk_target_celltype.py's docstring for the full design
# reasoning, including why this is NOT the same thing as Module 6's
# grn_celltype_specific.py (that computes correlation-based regulatory
# scores per cell type; it never forward-passes through {chr}_net.pt and
# produces no pseudobulk expression/accessibility matrix).
#
# SCOPE, DELIBERATELY: cell-type-resolved only. The held-out-dataset-
# specific variant (same --target-path/--opn-path mechanism, substituting a
# held-out dataset's pseudobulk instead of a cell type's) uses the identical
# code path but is NOT built here — deferred until this simpler version is
# validated per cell type, rather than building an untested second
# extension on top of a first one that hasn't been checked against real
# data yet.
#
# PER-CELL-TYPE SANITY CHECK IS MANDATORY, same discipline as Module 8's
# population-level gate: the population sanity check does not establish
# that a cell type's narrower input distribution is in-distribution for a
# net trained on pooled population statistics. Read each cell type's sanity
# check log before trusting its knockout output — this module structures
# the DAG so the sanity check runs before any knockout rule for the same
# cell type, but (same limitation as Module 8) Snakemake can only gate on
# the sanity check's output file existing, not on its content passing.
# =============================================================================

import re

MODULE8B_DIR = f"{MODULE8_DIR}/celltype"
PERTURB_CT_CFG = config["perturbation"].get("celltype_knockouts", {})
CT_MIN_CELLS = PERTURB_CT_CFG.get("min_cells_for_pseudobulk", 100)
_enabled_cfg = PERTURB_CT_CFG.get("enabled_celltypes", "all")


def _celltype_knockout_types():
    """Cell types eligible for Module 8b, filtered by real cell counts from
    cluster_annotation.tsv (same file Module 6's CELL_TYPES already reads)
    rather than trusting config alone — a cell type configured 'enabled'
    that has too few cells still won't get a rule instantiated for it here.

    ASSUMPTION, not confirmed against a real post-fix cluster_annotation.tsv:
    the file keeps the 'n_cells'/'cell_type' columns integrate_harmony.py's
    cluster_template.tsv originally wrote (one row per LEIDEN cluster, several
    clusters can share one cell_type). If a real file is missing 'n_cells',
    this falls back to including every cell type from CELL_TYPES uncounted —
    the runtime --min-cells guard in prep_pseudobulk_target_celltype.py is
    the backstop either way, so this is a DAG-time convenience, not the only
    enforcement of the threshold.
    """
    ann_path = INTEGRATION_CFG["cluster_annotation_tsv"]
    if not os.path.exists(ann_path):
        return []
    df = pd.read_csv(ann_path, sep="\t")
    if "n_cells" in df.columns and "cell_type" in df.columns:
        counts = df.groupby("cell_type")["n_cells"].sum()
        types = sorted(counts[counts >= CT_MIN_CELLS].index.tolist())
    else:
        print(f"WARNING: {ann_path} has no 'n_cells' column — cannot pre-filter by cell "
              f"count at DAG-build time, falling back to all of CELL_TYPES. The runtime "
              f"--min-cells guard in prep_pseudobulk_target_celltype.py still applies.")
        types = list(CELL_TYPES)
    if _enabled_cfg not in ("all", None, []):
        types = [c for c in types if c in _enabled_cfg]
    return types


CELLTYPE_KO_TYPES = _celltype_knockout_types()

# FIX 2026-09-18 — same AmbiguousRuleException class as 08_perturbation.smk's
# ko_id fix (see that file's comment for the full mechanism). ko_id here was
# equally unconstrained; constraining it too closes off the identical
# collision the moment anything else ever nests under
# MODULE8B_DIR/{pert_celltype}/... in the future, not just today's specific
# clash with linger_perturbation_ko.
wildcard_constraints:
    pert_celltype = "|".join(re.escape(c) for c in CELLTYPE_KO_TYPES) if CELLTYPE_KO_TYPES else "(?!)",
    ko_id = "|".join(re.escape(k) for k in KNOCKOUTS.keys()) if KNOCKOUTS else "(?!)"


rule prep_pseudobulk_target_celltype:
    """Cell-type-restricted TG/RE pseudobulk — see script docstring. May
    write nothing (exit 0) if the real post-shared-barcode cell count for
    this cell type falls below --min-cells; the DAG-time filter above tries
    to avoid instantiating this rule for such a cell type in the first
    place, but the runtime guard is the real enforcement."""
    input:
        labeled = f"{SCRATCH}/module3_integration/labeled.h5ad",
        atac_consensus = expand(f"{SCRATCH}/{{sample}}/atac_consensus.h5ad", sample=SAMPLES),
    output:
        tg = f"{MODULE8B_DIR}/{{pert_celltype}}/TG_pseudobulk_{{pert_celltype}}.tsv",
        re_ = f"{MODULE8B_DIR}/{{pert_celltype}}/RE_pseudobulk_{{pert_celltype}}.tsv",
    params:
        sample_ids = SAMPLES,
        min_cells = CT_MIN_CELLS,
        outdir = lambda wc: f"{MODULE8B_DIR}/{wc.pert_celltype}",
    log:
        f"{SCRATCH}/logs/08b_prep_pseudobulk_{{pert_celltype}}.log",
    resources:
        # FIX 2026-09-18 — was reusing prep_pseudobulk_target's (16GB) resource
        # block, sized for a lightweight file-read, not for this rule's real
        # pseudo_bulk() call (same memory profile as linger_prep_pseudobulk,
        # which needs 100GB — see config_eddie.yaml's prep_pseudobulk_target_
        # celltype entry for the full incident writeup: 19/20 array tasks were
        # SIGKILLed under the old 16GB cap, confirmed via empty per-rule logs).
        runtime   = config["resources"]["prep_pseudobulk_target_celltype"]["runtime_min"],
        sge_extra = sge_extra("prep_pseudobulk_target_celltype"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        export PYTHONHASHSEED=0
        {LINGER_PYTHON} workflow/scripts/prep_pseudobulk_target_celltype.py \
            --labeled {input.labeled} \
            --atac-consensus {input.atac_consensus} \
            --sample-ids {params.sample_ids} \
            --celltype "{wildcards.pert_celltype}" \
            --min-cells {params.min_cells} \
            --output-dir "{params.outdir}"
        """


rule linger_perturbation_celltype_sanity_check:
    """RUN THIS FIRST for each cell type, same discipline as Module 8's
    population-level linger_perturbation_sanity_check. Read the log before
    trusting any knockout rule below for the same cell type."""
    input:
        exp = f"{MODULE8_DIR}/Exp.tsv",
        re_tglink = f"{MODULE8_DIR}/RE_TGlink_resolved.tsv",
        tg = f"{MODULE8B_DIR}/{{pert_celltype}}/TG_pseudobulk_{{pert_celltype}}.tsv",
        re_ = f"{MODULE8B_DIR}/{{pert_celltype}}/RE_pseudobulk_{{pert_celltype}}.tsv",
    output:
        pred = f"{MODULE8B_DIR}/{{pert_celltype}}/_sanity_check_baseline_predicted.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        module8_dir = MODULE8_DIR,
    log:
        f"{SCRATCH}/logs/08b_sanity_check_{{pert_celltype}}.log",
    resources:
        runtime   = config["resources"]["linger_perturbation_ko"]["runtime_min"],
        sge_extra = sge_extra("linger_perturbation_ko"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        export PYTHONHASHSEED=0
        {LINGER_PYTHON} workflow/scripts/linger_perturbation.py \
            --workdir {params.workdir} \
            --module8-dir {params.module8_dir} \
            --ko-id "_sanity_check_baseline_{wildcards.pert_celltype}" \
            --tf-list \
            --sanity-check \
            --target-path "{input.tg}" \
            --opn-path "{input.re_}" \
            --output-tsv "{output.pred}"
        echo ""
        echo "=== {wildcards.pert_celltype}: read the SANITY CHECK block above before trusting any knockout output for this cell type. ==="
        """


rule linger_perturbation_celltype_ko:
    input:
        exp = f"{MODULE8_DIR}/Exp.tsv",
        re_tglink = f"{MODULE8_DIR}/RE_TGlink_resolved.tsv",
        tg = f"{MODULE8B_DIR}/{{pert_celltype}}/TG_pseudobulk_{{pert_celltype}}.tsv",
        re_ = f"{MODULE8B_DIR}/{{pert_celltype}}/RE_pseudobulk_{{pert_celltype}}.tsv",
        # Same pattern as Module 8: gates on the sanity check having RUN
        # (output file exists), not on its content having passed — that
        # judgment is manual, read the log.
        sanity_check = f"{MODULE8B_DIR}/{{pert_celltype}}/_sanity_check_baseline_predicted.tsv",
    output:
        pred = f"{MODULE8B_DIR}/{{pert_celltype}}/{{ko_id}}_predicted_expression.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        module8_dir = MODULE8_DIR,
        tf_list = lambda wc: KNOCKOUTS[wc.ko_id],
    log:
        f"{SCRATCH}/logs/08b_perturb_{{pert_celltype}}_{{ko_id}}.log",
    resources:
        runtime   = config["resources"]["linger_perturbation_ko"]["runtime_min"],
        sge_extra = sge_extra("linger_perturbation_ko"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        export PYTHONHASHSEED=0
        {LINGER_PYTHON} workflow/scripts/linger_perturbation.py \
            --workdir {params.workdir} \
            --module8-dir {params.module8_dir} \
            --ko-id "{wildcards.pert_celltype}_{wildcards.ko_id}" \
            --tf-list {params.tf_list} \
            --target-path "{input.tg}" \
            --opn-path "{input.re_}" \
            --output-tsv "{output.pred}"
        """


rule module8b_all:
    input:
        expand(
            f"{MODULE8B_DIR}/{{pert_celltype}}/_sanity_check_baseline_predicted.tsv",
            pert_celltype=CELLTYPE_KO_TYPES,
        ),
        expand(
            f"{MODULE8B_DIR}/{{pert_celltype}}/{{ko_id}}_predicted_expression.tsv",
            pert_celltype=CELLTYPE_KO_TYPES,
            ko_id=list(KNOCKOUTS.keys()),
        ),
