# =============================================================================
# workflow/rules/08_perturbation.smk
# Module 8 — In silico perturbation (TF knockout), BYPASS design
#
# DECIDED 2026-09-08, not a design question anymore. Three options were on
# the table (see README's Module 8 design notes for the full writeup):
#   (1) check LINGER's GitHub/tutorials for a working mouse-mode
#       perturbation example — CHECKED, real dead end: docs/perturb.md's
#       only example is human/hg19 (H1 cell line, GRNdir='data_bulk/'),
#       and docs/scNN.md (the mouse/"other species" tutorial) stops at TF
#       activity — no perturbation section exists there at all. No GitHub
#       issue found addressing mouse-mode perturbation either.
#   (2) write a compatibility shim for perturb.py's four missing files —
#       REJECTED. Reading perturb.py's real source turned up three more
#       bugs beyond the already-known filename mismatch (hardcoded human
#       chr1-22+X chromosome list; a Target/data_merge row-alignment
#       assumption scNN's RE-linked gene subset violates; an explicit-
#       index TF/RE gather scheme that doesn't match scNN training's real
#       implicit "all TFs except self" scheme). A shim would mean
#       reimplementing scNN's real input contract anyway, then translating
#       it into a format with three more bugs layered on top.
#   (3) bypass perturb.py entirely and reimplement the forward pass
#       directly against {chr}_net.pt — CHOSEN. Modeled exactly on
#       LINGER_tr.sc_nn_NN(), the real function that trained those nets,
#       not on perturb.py's LINGER_simulation(). See
#       workflow/scripts/linger_perturbation.py's docstring for the full
#       per-gene input-reconstruction detail.
#
# FLAGGED, NOT SILENTLY ASSUMED SAFE: reconstructing which TFs went into
# each gene's trained net depends on a Python set() intersection
# (LINGER_tr.load_data_scNN()'s TFlist construction) whose iteration order
# is process-dependent unless PYTHONHASHSEED was fixed at population
# training time — that order was never persisted to disk in scNN mode. See
# prep_pseudobulk_target.py's docstring for the full explanation and
# linger_perturbation.py's --sanity-check mode, which is the only real
# check available for this. RUN linger_perturbation_sanity_check AND READ
# ITS LOG before trusting any real knockout output below.
#
# No overexpression mode is built here — per the HTML plan doc's own
# Module 9 spec ("Positive: Atoh1/Gfi1/Pou4f3 prediction vs. real
# overexpression data"), Module 9 compares this module's KNOCKOUT
# predictions directly against real overexpression ground truth (expecting
# an inverse relationship), rather than needing a separate gain-of-
# function simulation mode here.
# =============================================================================

MODULE8_DIR = f"{SCRATCH}/module8_perturbation"
KNOCKOUTS = config["perturbation"]["knockouts"]  # dict: ko_id -> [TF, TF, ...]


rule prep_pseudobulk_target:
    """One-time materialization of Exp/RE_TGlink in a fixed, self-consistent
    order, shared by every knockout run below (see script docstring)."""
    input:
        population_done = f"{SCRATCH}/module6_grn/population_training.done",
    output:
        exp = f"{MODULE8_DIR}/Exp.tsv",
        re_tglink = f"{MODULE8_DIR}/RE_TGlink_resolved.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        grn_dir = LINGER_CFG["grn_dir"],
        genome  = LINGER_CFG["genome"],
        outdir  = MODULE8_DIR,
    log:
        f"{SCRATCH}/logs/08a_prep_pseudobulk_target.log",
    resources:
        runtime   = config["resources"]["prep_pseudobulk_target"]["runtime_min"],
        sge_extra = sge_extra("prep_pseudobulk_target"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/prep_pseudobulk_target.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --output-dir {params.outdir}
        """


rule linger_perturbation_sanity_check:
    """RUN THIS FIRST. Baseline forward pass with NO knockout applied,
    scored against the real pseudobulk Target — the only available check
    on the TF-ordering risk flagged in prep_pseudobulk_target.py's
    docstring. Read the log before trusting any rule below."""
    input:
        exp = f"{MODULE8_DIR}/Exp.tsv",
        re_tglink = f"{MODULE8_DIR}/RE_TGlink_resolved.tsv",
    output:
        pred = f"{MODULE8_DIR}/_sanity_check_baseline_predicted.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        module8_dir = MODULE8_DIR,
    log:
        f"{SCRATCH}/logs/08b_sanity_check.log",
    resources:
        runtime   = config["resources"]["linger_perturbation_ko"]["runtime_min"],
        sge_extra = sge_extra("linger_perturbation_ko"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/linger_perturbation.py \
            --workdir {params.workdir} \
            --module8-dir {params.module8_dir} \
            --ko-id _sanity_check_baseline \
            --tf-list \
            --sanity-check \
            --output-tsv {output.pred}
        echo ""
        echo "=== Read the SANITY CHECK block above before trusting any real knockout output. ==="
        """


rule linger_perturbation_ko:
    input:
        exp = f"{MODULE8_DIR}/Exp.tsv",
        re_tglink = f"{MODULE8_DIR}/RE_TGlink_resolved.tsv",
        # Depending on the sanity check's OUTPUT FILE (not gating on its
        # content — Snakemake can't inspect a log for pass/fail) at least
        # guarantees it ran first and its log exists to be read. The
        # actual go/no-go judgment on trusting knockout output is manual —
        # see this rule file's module docstring.
        sanity_check = f"{MODULE8_DIR}/_sanity_check_baseline_predicted.tsv",
    output:
        pred = f"{MODULE8_DIR}/{{ko_id}}_predicted_expression.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        module8_dir = MODULE8_DIR,
        tf_list = lambda wc: KNOCKOUTS[wc.ko_id],
    log:
        f"{SCRATCH}/logs/08c_perturb_{{ko_id}}.log",
    resources:
        runtime   = config["resources"]["linger_perturbation_ko"]["runtime_min"],
        sge_extra = sge_extra("linger_perturbation_ko"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/linger_perturbation.py \
            --workdir {params.workdir} \
            --module8-dir {params.module8_dir} \
            --ko-id {wildcards.ko_id} \
            --tf-list {params.tf_list} \
            --output-tsv {output.pred}
        """


rule module8_all:
    input:
        f"{MODULE8_DIR}/_sanity_check_baseline_predicted.tsv",
        expand(f"{MODULE8_DIR}/{{ko_id}}_predicted_expression.tsv", ko_id=list(KNOCKOUTS.keys())),
