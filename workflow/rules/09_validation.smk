# =============================================================================
# workflow/rules/09_validation.smk
# Module 9 — Validation against held-out data
#
# Per the HTML plan doc: "the credibility check the whole approach depends
# on... predictions are checked against real TF-overexpression experiments
# where mature cells were converted toward hair-cell identity. If LINGER
# reproduces that known result, it's reasonable to trust its extrapolation
# into the aging context, where there's no ground truth to check against
# directly." Output stated explicitly as "Correlation metrics, AUROC/AUPR"
# — see validate_held_out.py's docstring for how both are computed.
#
# CHECK_SPECS maps each config `check` value (held_out_inputs /
# negative_control_input) to how it's scored:
#   - atoh1_gfi1_pou4f3_overexpression: compares Module 8's triple_ko
#     prediction (sign-flipped: knockout simulates loss, this held-out
#     data shows gain) directly against GSE224627's real profile — per
#     the HTML's own spec ("...prediction vs. real overexpression data"),
#     NOT a separately-built overexpression simulation mode.
#   - tbx2_conversion: Module 8's tbx2_ko prediction (no sign flip — this
#     is a real loss-of-function conversion) against GSE233559.
#   - aging_vector / aging_vector_support: Module 7's role="baseline" mean
#     TF activity as the "unperturbed" reference, diffed against each
#     aging dataset's own TF activity (computed the same way Module 7
#     computes it). No Module 8 dependency for these four.
#   - negative_control_must_not_reprogram: same machinery as the
#     overexpression check, but PASS means the correlation/AUROC signal
#     should NOT look reprogramming-consistent (invert_pass=True in
#     aggregate_validation.py's pass/fail logic).
# =============================================================================

MODULE9_DIR = f"{SCRATCH}/module9_validation"

CHECK_SPECS = {
    "atoh1_gfi1_pou4f3_overexpression": {"mode": "expression_shift", "ko_id": "triple_ko", "sign": -1, "invert_pass": False, "top_k": 200},
    "tbx2_conversion":                  {"mode": "expression_shift", "ko_id": "tbx2_ko",   "sign": 1,  "invert_pass": False, "top_k": 200},
    "aging_vector":                     {"mode": "aging_tf_activity", "ko_id": None, "sign": None, "invert_pass": False, "top_k": 200},
    "aging_vector_support":             {"mode": "aging_tf_activity", "ko_id": None, "sign": None, "invert_pass": False, "top_k": 200},
    "negative_control_must_not_reprogram": {"mode": "expression_shift", "ko_id": "triple_ko", "sign": -1, "invert_pass": True, "top_k": 200},
}

_HELD_OUT_ENTRIES = {e["sample_id"]: e for e in HELD_OUT_CFG}
_HELD_OUT_ENTRIES[NEGATIVE_CONTROL_CFG["sample_id"]] = NEGATIVE_CONTROL_CFG
VALIDATION_SAMPLE_IDS = list(_HELD_OUT_ENTRIES.keys())


def _validation_entry_json(wildcards):
    import json
    return json.dumps(_HELD_OUT_ENTRIES[wildcards.sample_id])


def _validation_check_spec_json(wildcards):
    import json
    check = _HELD_OUT_ENTRIES[wildcards.sample_id]["check"]
    return json.dumps(CHECK_SPECS[check])


def _validation_inputs(wildcards):
    """Only depend on the specific Module 8 knockout this check actually
    needs (aging checks need none) — avoids forcing every knockout to
    finish before any validation can run."""
    check = _HELD_OUT_ENTRIES[wildcards.sample_id]["check"]
    spec = CHECK_SPECS[check]
    inputs = {
        "module7_summary": f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
    }
    if spec["ko_id"] is not None:
        inputs["ko_pred"] = f"{SCRATCH}/module8_perturbation/{spec['ko_id']}_predicted_expression.tsv"
        inputs["sanity_check"] = f"{SCRATCH}/module8_perturbation/_sanity_check_baseline_predicted.tsv"
    return inputs


rule validate_one:
    input:
        unpack(_validation_inputs),
    output:
        score = f"{MODULE9_DIR}/{{sample_id}}_score.tsv",
    params:
        workdir  = f"{SCRATCH}/module4_linger_init",
        grn_dir  = LINGER_CFG["grn_dir"],
        genome   = LINGER_CFG["genome"],
        module8_dir = f"{SCRATCH}/module8_perturbation",
        module7_summary = f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
        entry_json = _validation_entry_json,
        check_spec_json = _validation_check_spec_json,
    log:
        f"{SCRATCH}/logs/09a_validate_{{sample_id}}.log",
    resources:
        runtime   = config["resources"]["validate_one"]["runtime_min"],
        sge_extra = sge_extra("validate_one"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/validate_held_out.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --module8-dir {params.module8_dir} \
            --module7-summary {params.module7_summary} \
            --entry-json {params.entry_json:q} \
            --check-spec-json {params.check_spec_json:q} \
            --output-tsv {output.score}
        """


rule module9_aggregate:
    input:
        expand(f"{MODULE9_DIR}/{{sample_id}}_score.tsv", sample_id=VALIDATION_SAMPLE_IDS),
    output:
        report = f"{MODULE9_DIR}/validation_report.tsv",
    log:
        f"{SCRATCH}/logs/09b_aggregate.log",
    resources:
        runtime   = config["resources"]["module9_aggregate"]["runtime_min"],
        sge_extra = sge_extra("module9_aggregate"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        {LINGER_PYTHON} workflow/scripts/aggregate_validation.py \
            --module9-dir {MODULE9_DIR} \
            --output-tsv {output.report}
        """


rule module9_all:
    input:
        f"{MODULE9_DIR}/validation_report.tsv",
