# =============================================================================
# workflow/rules/09b_validation_celltype.smk
# Module 9b — cell-type-resolved validation against held-out data
#
# ADDED 2026-09-19. NOT a replacement for Module 9 (09_validation.smk) — a
# parallel, equally-weighted second pass over the SAME held-out/negative-
# control/aging datasets, scored against Module 8b's cell-type-resolved
# knockout predictions and Module 6's cell-type-specific regulatory
# networks instead of Module 8's/Module 7's population-level ones. Per an
# explicit design decision this session: population and cell-type results
# are reported side by side with an identical schema (see
# aggregate_validation_celltype.py's `scope` column), NEITHER treated as
# the primary result — Module 9's validation_report.tsv is untouched by
# this file.
#
# REUSES validate_held_out.py UNCHANGED in logic — only two new optional
# CLI flags were added to it (--baseline-tsv, --network), both defaulting
# to the exact population behavior Module 9 already relies on. See that
# script's ADDED 2026-09-19 docstring section for the full reasoning,
# including the flagged caveat that every held-out dataset here is BULK
# RNA-seq: a cell-type-resolved run changes which MODEL the same bulk
# sample is scored against, not which cells are used — there is no real
# per-cell-type ground truth in this dataset collection. State that
# plainly wherever these results are reported, the same way Module 8's
# PYTHONHASHSEED risk and Module 9's own weak result are already flagged
# rather than glossed over.
#
# SCOPE: only the three checks that actually consume a Module 8 knockout
# prediction (atoh1_gfi1_pou4f3_overexpression, tbx2_conversion,
# negative_control_must_not_reprogram) AND the two aging checks
# (aging_vector, aging_vector_support) are replicated per cell type — all
# five of Module 9's checks, not an artificial subset. Confirmed via
# LingerGRN 1.110's real TF_activity.regulon() source (not assumed): any
# `network` argument other than "cell population"/"general" is read as
# `cell_type_specific_trans_regulatory_{network}.txt`, which Module 6
# (grn_celltype_specific.py) already writes, one per cell type.
#
# CELL-TYPE SCOPE: Module 8b's CELLTYPE_KO_TYPES (>=100-cell threshold,
# from 08b_perturbation_celltype.smk) is the base list for every check
# here, including the aging ones, so one cell-type axis is shared across
# all five checks in the combined report. Module 6's CELL_TYPES has no
# cell-count threshold and is a superset in principle, but this is
# verified per cell type at DAG-build time (_celltype_aging_types, below)
# rather than assumed — a cell type dropped by Module 6 for any other
# reason simply gets no aging rows, with a build-time warning, rather than
# a hard Snakemake failure.
# =============================================================================

import json as _json9b

MODULE9B_DIR = f"{SCRATCH}/module9b_validation_celltype"


def _celltype_aging_types():
    """CELLTYPE_KO_TYPES cell types for which Module 6 actually wrote
    cell_type_specific_trans_regulatory_{celltype}.txt into module4_linger_init/
    — checked directly against real files on disk, not assumed from Module
    6's own (unfiltered, no cell-count threshold) CELL_TYPES list."""
    m4dir = f"{SCRATCH}/module4_linger_init"
    available, missing = [], []
    for ct in CELLTYPE_KO_TYPES:
        if os.path.exists(os.path.join(m4dir, f"cell_type_specific_trans_regulatory_{ct}.txt")):
            available.append(ct)
        else:
            missing.append(ct)
    if missing:
        print(f"WARNING: 09b_validation_celltype — {len(missing)}/{len(CELLTYPE_KO_TYPES)} "
              f"Module 8b cell types have no Module 6 cell_type_specific_trans_regulatory_ "
              f"*.txt file (stale Module 6 run, or a naming mismatch against the current "
              f"cluster_annotation.tsv). These get NO aging_vector/aging_vector_support rows "
              f"in Module 9b (knockout-check rows are unaffected — those only need Module "
              f"8b, not Module 6): {missing}")
    return available


CELLTYPE_AGING_TYPES = _celltype_aging_types()

# Same three checks as Module 9's OTHER_SAMPLE_IDS (everything except the
# aging reference/support roles) — reused directly from 09_validation.smk,
# already evaluated in this same Snakemake namespace via `include:`.
CELLTYPE_KNOCKOUT_SAMPLE_IDS = OTHER_SAMPLE_IDS


def _ct_validation_inputs(wildcards):
    # NOTE: this is a function passed to `input: unpack(...)`, so it runs
    # per-job with a real `wildcards` object already bound — every path
    # below must be fully resolved here using wildcards.pert_celltype /
    # wildcards.sample_id. Unlike a static `input:`/`output:` string (which
    # Snakemake wildcard-resolves itself), a string returned from an input
    # function is used exactly as given — leaving a literal "{pert_celltype}"
    # placeholder in it would NOT get substituted and would silently point
    # at a nonexistent path.
    check = _HELD_OUT_ENTRIES[wildcards.sample_id]["check"]
    spec = CHECK_SPECS[check]
    ct = wildcards.pert_celltype
    inputs = {
        "module7_summary": f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
        "ko_pred": f"{MODULE8B_DIR}/{ct}/{spec['ko_id']}_predicted_expression.tsv",
        "sanity_check": f"{MODULE8B_DIR}/{ct}/_sanity_check_baseline_predicted.tsv",
        "tg_baseline": f"{MODULE8B_DIR}/{ct}/TG_pseudobulk_{ct}.tsv",
    }
    return inputs


rule validate_celltype_one:
    wildcard_constraints:
        pert_celltype = "|".join(re.escape(c) for c in CELLTYPE_KO_TYPES) if CELLTYPE_KO_TYPES else "(?!)",
        sample_id = _OTHER_SAMPLE_ID_PATTERN,
    input:
        unpack(_ct_validation_inputs),
    output:
        score = f"{MODULE9B_DIR}/{{pert_celltype}}/{{sample_id}}_score.tsv",
    params:
        workdir  = f"{SCRATCH}/module4_linger_init",
        grn_dir  = LINGER_CFG["grn_dir"],
        genome   = LINGER_CFG["genome"],
        module8_dir = lambda wc: f"{MODULE8B_DIR}/{wc.pert_celltype}",
        module7_summary = f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
        entry_json = _validation_entry_json,
        check_spec_json = _validation_check_spec_json,
        baseline_tsv = lambda wc: f"{MODULE8B_DIR}/{wc.pert_celltype}/TG_pseudobulk_{wc.pert_celltype}.tsv",
    log:
        f"{SCRATCH}/logs/09x_validate_celltype_{{pert_celltype}}_{{sample_id}}.log",
    resources:
        runtime   = config["resources"]["validate_celltype_one"]["runtime_min"],
        sge_extra = sge_extra("validate_celltype_one"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/validate_held_out.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --module8-dir "{params.module8_dir}" \
            --module7-summary {params.module7_summary} \
            --entry-json {params.entry_json:q} \
            --check-spec-json {params.check_spec_json:q} \
            --baseline-tsv "{params.baseline_tsv}" \
            --output-tsv "{output.score}"
        """


rule validate_celltype_aging_vector:
    wildcard_constraints:
        pert_celltype = "|".join(re.escape(c) for c in CELLTYPE_AGING_TYPES) if CELLTYPE_AGING_TYPES else "(?!)",
    input:
        module7_summary = f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
        # Real Snakemake dependency on Module 6's cell-type network, so this
        # rule can't start before it exists — mirrors validate_aging_vector's
        # population dependency on Module 7's summary.
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_type_specific_trans_regulatory_{{pert_celltype}}.txt",
    output:
        score = f"{MODULE9B_DIR}/{{pert_celltype}}/{AGING_VECTOR_SAMPLE_ID}_score.tsv",
        shift = f"{MODULE9B_DIR}/{{pert_celltype}}/{AGING_VECTOR_SAMPLE_ID}_aging_shift.tsv",
    params:
        workdir  = f"{SCRATCH}/module4_linger_init",
        grn_dir  = LINGER_CFG["grn_dir"],
        genome   = LINGER_CFG["genome"],
        module8_dir = lambda wc: f"{MODULE8B_DIR}/{wc.pert_celltype}",
        module7_summary = f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
        entry_json = _json9b.dumps(_HELD_OUT_ENTRIES[AGING_VECTOR_SAMPLE_ID]),
        check_spec_json = _json9b.dumps(CHECK_SPECS["aging_vector"]),
    log:
        f"{SCRATCH}/logs/09x_validate_celltype_{{pert_celltype}}_{AGING_VECTOR_SAMPLE_ID}.log",
    resources:
        runtime   = config["resources"]["validate_celltype_one"]["runtime_min"],
        sge_extra = sge_extra("validate_celltype_one"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/validate_held_out.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --module8-dir "{params.module8_dir}" \
            --module7-summary {params.module7_summary} \
            --entry-json {params.entry_json:q} \
            --check-spec-json {params.check_spec_json:q} \
            --network "{wildcards.pert_celltype}" \
            --output-tsv "{output.score}" \
            --output-shift-tsv "{output.shift}"
        """


rule validate_celltype_aging_support:
    wildcard_constraints:
        pert_celltype = "|".join(re.escape(c) for c in CELLTYPE_AGING_TYPES) if CELLTYPE_AGING_TYPES else "(?!)",
        sample_id = _AGING_SUPPORT_ID_PATTERN,
    input:
        module7_summary = f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_type_specific_trans_regulatory_{{pert_celltype}}.txt",
        aging_vector_shift = f"{MODULE9B_DIR}/{{pert_celltype}}/{AGING_VECTOR_SAMPLE_ID}_aging_shift.tsv",
    output:
        score = f"{MODULE9B_DIR}/{{pert_celltype}}/{{sample_id}}_score.tsv",
    params:
        workdir  = f"{SCRATCH}/module4_linger_init",
        grn_dir  = LINGER_CFG["grn_dir"],
        genome   = LINGER_CFG["genome"],
        module8_dir = lambda wc: f"{MODULE8B_DIR}/{wc.pert_celltype}",
        module7_summary = f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
        entry_json = lambda wildcards: _json9b.dumps(_HELD_OUT_ENTRIES[wildcards.sample_id]),
        check_spec_json = _json9b.dumps(CHECK_SPECS["aging_vector_support"]),
    log:
        f"{SCRATCH}/logs/09x_validate_celltype_{{pert_celltype}}_{{sample_id}}.log",
    resources:
        runtime   = config["resources"]["validate_celltype_one"]["runtime_min"],
        sge_extra = sge_extra("validate_celltype_one"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/validate_held_out.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --module8-dir "{params.module8_dir}" \
            --module7-summary {params.module7_summary} \
            --entry-json {params.entry_json:q} \
            --check-spec-json {params.check_spec_json:q} \
            --network "{wildcards.pert_celltype}" \
            --output-tsv "{output.score}" \
            --aging-vector-tsv "{input.aging_vector_shift}"
        """


def _module9b_score_targets():
    targets = []
    for ct in CELLTYPE_KO_TYPES:
        for sid in CELLTYPE_KNOCKOUT_SAMPLE_IDS:
            targets.append(f"{MODULE9B_DIR}/{ct}/{sid}_score.tsv")
    for ct in CELLTYPE_AGING_TYPES:
        targets.append(f"{MODULE9B_DIR}/{ct}/{AGING_VECTOR_SAMPLE_ID}_score.tsv")
        for sid in AGING_SUPPORT_SAMPLE_IDS:
            targets.append(f"{MODULE9B_DIR}/{ct}/{sid}_score.tsv")
    return targets


rule module9b_aggregate:
    input:
        celltype_scores = _module9b_score_targets(),
        # Real dependency on Module 9's own population report, so the
        # combined report can never be built from a stale/missing one.
        population_report = f"{MODULE9_DIR}/validation_report.tsv",
    output:
        celltype_report = f"{MODULE9B_DIR}/validation_report_celltype.tsv",
        combined_report = f"{MODULE9B_DIR}/validation_report_combined.tsv",
    params:
        module8_dir = f"{SCRATCH}/module8_perturbation",
        module8b_dir = MODULE8B_DIR,
        module4_data = f"{SCRATCH}/module4_linger_init/data",
    log:
        f"{SCRATCH}/logs/09x_aggregate.log",
    resources:
        runtime   = config["resources"]["module9b_aggregate"]["runtime_min"],
        sge_extra = sge_extra("module9b_aggregate"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        {LINGER_PYTHON} workflow/scripts/aggregate_validation_celltype.py \
            --module9b-dir {MODULE9B_DIR} \
            --population-report {input.population_report} \
            --module8-dir {params.module8_dir} \
            --module8b-dir {params.module8b_dir} \
            --module4-data {params.module4_data} \
            --output-celltype-tsv {output.celltype_report} \
            --output-combined-tsv {output.combined_report}
        """


rule module9b_all:
    input:
        f"{MODULE9B_DIR}/validation_report_combined.tsv",
