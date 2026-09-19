# =============================================================================
# workflow/rules/11_caps_score.smk
# Module 11 — CAPS (Cochlear Aging Perturbation Score)
#
# ADDED 2026-09-19. CAPS(TF, celltype) = sign_adjusted_shift x
# regulatory_centrality x confidence_weight — formula approved this session,
# after researching 3 precedent papers (CellOracle/Kamimoto et al. 2023 is
# the closest match for sign_adjusted_shift's inner-product design). Built
# on the "more rigorous" TF-activity-space option for sign_adjusted_shift:
# each knockout's predicted expression (Module 8/8b) is run through the same
# TF_activity.regulon() call Module 9/9b already trust, collapsed to a
# TF-activity vector, and diffed against the equivalent pre-KO baseline —
# NOT a simpler gene-space overlap against the aging vector's genes.
#
# EQUAL TREATMENT, same pattern as Module 9b: population and every cell type
# in CELLTYPE_AGING_TYPES (09b_validation_celltype.smk's already-verified
# subset of CELLTYPE_KO_TYPES that has both a Module 6 trans-regulatory file
# AND an aging-vector shift file) get one row per ko_id, under one identical
# schema with a `scope` column. A cell type Module 9b can't run aging checks
# for (no Module 6 file) gets no CAPS rows either — same real on-disk check,
# not re-derived here.
#
# NOTHING HERE RE-RUNS A FORWARD PASS OR RE-COMPUTES A GRN: every input is
# an existing Module 6/8/8b/9/9b output, read as-is. See
# compute_tf_activity_shift.py and compute_caps_score.py docstrings for the
# full per-term reasoning (why the sanity-check baseline, not real data, is
# diffed; why regulatory_centrality is regulon()'s own colsum term read
# directly; why confidence_weight falls back to 0.0 rather than skipping a
# scope with no Module 9b sanity_rho).
# =============================================================================

MODULE11_DIR = f"{SCRATCH}/module11_caps"

CAPS_KO_IDS = list(KNOCKOUTS.keys())

wildcard_constraints:
    caps_ko_id = "|".join(re.escape(k) for k in CAPS_KO_IDS) if CAPS_KO_IDS else "(?!)",
    caps_celltype = "|".join(re.escape(c) for c in CELLTYPE_AGING_TYPES) if CELLTYPE_AGING_TYPES else "(?!)",


# ── population scope ────────────────────────────────────────────────────────

rule caps_tf_activity_baseline_population:
    input:
        expr = f"{MODULE8_DIR}/_sanity_check_baseline_predicted.tsv",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_population_trans_regulatory.txt",
    output:
        f"{MODULE11_DIR}/population/_baseline_tf_activity.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        grn_dir = LINGER_CFG["grn_dir"],
        genome  = LINGER_CFG["genome"],
    log:
        f"{SCRATCH}/logs/11a_tf_activity_population_baseline.log",
    resources:
        runtime   = config["resources"]["caps_tf_activity"]["runtime_min"],
        sge_extra = sge_extra("caps_tf_activity"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/compute_tf_activity_shift.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --network "cell population" \
            --expression-tsv "{input.expr}" \
            --output-tsv "{output}"
        """


rule caps_tf_activity_ko_population:
    input:
        expr = f"{MODULE8_DIR}/{{caps_ko_id}}_predicted_expression.tsv",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_population_trans_regulatory.txt",
    output:
        f"{MODULE11_DIR}/population/{{caps_ko_id}}_tf_activity.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        grn_dir = LINGER_CFG["grn_dir"],
        genome  = LINGER_CFG["genome"],
    log:
        f"{SCRATCH}/logs/11a_tf_activity_population_{{caps_ko_id}}.log",
    resources:
        runtime   = config["resources"]["caps_tf_activity"]["runtime_min"],
        sge_extra = sge_extra("caps_tf_activity"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/compute_tf_activity_shift.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --network "cell population" \
            --expression-tsv "{input.expr}" \
            --output-tsv "{output}"
        """


rule caps_score_population:
    input:
        baseline_tf_activity = f"{MODULE11_DIR}/population/_baseline_tf_activity.tsv",
        ko_tf_activity = f"{MODULE11_DIR}/population/{{caps_ko_id}}_tf_activity.tsv",
        aging_shift = f"{MODULE9_DIR}/{AGING_VECTOR_SAMPLE_ID}_aging_shift.tsv",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_population_trans_regulatory.txt",
        validation_report_combined = f"{MODULE9B_DIR}/validation_report_combined.tsv",
    output:
        f"{MODULE11_DIR}/population/{{caps_ko_id}}_caps.tsv",
    params:
        tf_list = lambda wc: KNOCKOUTS[wc.caps_ko_id],
    log:
        f"{SCRATCH}/logs/11b_caps_population_{{caps_ko_id}}.log",
    resources:
        runtime   = config["resources"]["caps_score"]["runtime_min"],
        sge_extra = sge_extra("caps_score"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/compute_caps_score.py \
            --ko-id "{wildcards.caps_ko_id}" \
            --scope "population" \
            --tf-list {params.tf_list} \
            --baseline-tf-activity-tsv "{input.baseline_tf_activity}" \
            --ko-tf-activity-tsv "{input.ko_tf_activity}" \
            --aging-shift-tsv "{input.aging_shift}" \
            --trans-regulatory-tsv "{input.trans_regulatory}" \
            --validation-report-combined "{input.validation_report_combined}" \
            --output-tsv "{output}"
        """


# ── cell-type scopes ────────────────────────────────────────────────────────

rule caps_tf_activity_baseline_celltype:
    input:
        expr = f"{MODULE8B_DIR}/{{caps_celltype}}/_sanity_check_baseline_predicted.tsv",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_type_specific_trans_regulatory_{{caps_celltype}}.txt",
    output:
        f"{MODULE11_DIR}/{{caps_celltype}}/_baseline_tf_activity.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        grn_dir = LINGER_CFG["grn_dir"],
        genome  = LINGER_CFG["genome"],
        network = lambda wc: wc.caps_celltype,
    log:
        f"{SCRATCH}/logs/11a_tf_activity_celltype_{{caps_celltype}}_baseline.log",
    resources:
        runtime   = config["resources"]["caps_tf_activity"]["runtime_min"],
        sge_extra = sge_extra("caps_tf_activity"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/compute_tf_activity_shift.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --network "{params.network}" \
            --expression-tsv "{input.expr}" \
            --output-tsv "{output}"
        """


rule caps_tf_activity_ko_celltype:
    input:
        expr = f"{MODULE8B_DIR}/{{caps_celltype}}/{{caps_ko_id}}_predicted_expression.tsv",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_type_specific_trans_regulatory_{{caps_celltype}}.txt",
    output:
        f"{MODULE11_DIR}/{{caps_celltype}}/{{caps_ko_id}}_tf_activity.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        grn_dir = LINGER_CFG["grn_dir"],
        genome  = LINGER_CFG["genome"],
        network = lambda wc: wc.caps_celltype,
    log:
        f"{SCRATCH}/logs/11a_tf_activity_celltype_{{caps_celltype}}_{{caps_ko_id}}.log",
    resources:
        runtime   = config["resources"]["caps_tf_activity"]["runtime_min"],
        sge_extra = sge_extra("caps_tf_activity"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/compute_tf_activity_shift.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --network "{params.network}" \
            --expression-tsv "{input.expr}" \
            --output-tsv "{output}"
        """


rule caps_score_celltype:
    input:
        baseline_tf_activity = f"{MODULE11_DIR}/{{caps_celltype}}/_baseline_tf_activity.tsv",
        ko_tf_activity = f"{MODULE11_DIR}/{{caps_celltype}}/{{caps_ko_id}}_tf_activity.tsv",
        aging_shift = f"{MODULE9B_DIR}/{{caps_celltype}}/{AGING_VECTOR_SAMPLE_ID}_aging_shift.tsv",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_type_specific_trans_regulatory_{{caps_celltype}}.txt",
        validation_report_combined = f"{MODULE9B_DIR}/validation_report_combined.tsv",
    output:
        f"{MODULE11_DIR}/{{caps_celltype}}/{{caps_ko_id}}_caps.tsv",
    params:
        tf_list = lambda wc: KNOCKOUTS[wc.caps_ko_id],
    log:
        f"{SCRATCH}/logs/11b_caps_celltype_{{caps_celltype}}_{{caps_ko_id}}.log",
    resources:
        runtime   = config["resources"]["caps_score"]["runtime_min"],
        sge_extra = sge_extra("caps_score"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/compute_caps_score.py \
            --ko-id "{wildcards.caps_ko_id}" \
            --scope "{wildcards.caps_celltype}" \
            --tf-list {params.tf_list} \
            --baseline-tf-activity-tsv "{input.baseline_tf_activity}" \
            --ko-tf-activity-tsv "{input.ko_tf_activity}" \
            --aging-shift-tsv "{input.aging_shift}" \
            --trans-regulatory-tsv "{input.trans_regulatory}" \
            --validation-report-combined "{input.validation_report_combined}" \
            --output-tsv "{output}"
        """


# ── aggregate ────────────────────────────────────────────────────────────────

def _module11_caps_targets():
    targets = [f"{MODULE11_DIR}/population/{ko_id}_caps.tsv" for ko_id in CAPS_KO_IDS]
    for ct in CELLTYPE_AGING_TYPES:
        targets += [f"{MODULE11_DIR}/{ct}/{ko_id}_caps.tsv" for ko_id in CAPS_KO_IDS]
    return targets


rule module11_aggregate:
    input:
        caps_scores = _module11_caps_targets(),
    output:
        f"{MODULE11_DIR}/caps_scores.tsv",
    params:
        module11_dir = MODULE11_DIR,
    log:
        f"{SCRATCH}/logs/11c_aggregate.log",
    resources:
        runtime   = config["resources"]["module11_aggregate"]["runtime_min"],
        sge_extra = sge_extra("module11_aggregate"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/aggregate_caps_scores.py \
            --module11-dir "{params.module11_dir}" \
            --output-tsv "{output}"
        """


rule module11_all:
    input:
        f"{MODULE11_DIR}/caps_scores.tsv",
