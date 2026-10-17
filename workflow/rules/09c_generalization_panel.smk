# =============================================================================
# workflow/rules/09c_generalization_panel.smk
# Module 9c — generalisation test panel (ADDED 2026-10-15, Round B / Stage 3)
#
# Scores the generalization_test_panel entries that are eligible (Track-0 resolved in the config:
# `eligibility: eligible` or `effectors_eligible`) with the SAME rules as Module 9 (validate_held_out.py:
# gene-bootstrap excess Spearman rho over the level-only null; PASS if the CI's lower end clears
# min_excess_over_level_null, FAIL if its upper end is below it, else INCONCLUSIVE; NOT_EVALUABLE if the
# predicted shift is a no-op). POPULATION scope only. Entries whose gene is not in the trained network
# (`eligibility: not_in_network`: Casz1, Zbtb20, Prdm16) are never perturbed: they get a descriptive-only
# row in the report.
#
# COMBINATORIAL SUPPORT (GSE86204, Rbpj/Tcf7/Lef1): `combinatorial_arms` was config-only before. Each arm in the
# entry's `scoring_arms` becomes ONE ordinary entry whose target columns come from that arm's
# sample_column_regex and whose paired baseline is the control arm's columns, so no loader change is
# needed. An arm's predicted shift is the additive composite sum_i sign_i * (KO_i - baseline) over its
# `components`; one arm is also scored against a single joint Rbpj+Tcf7+Lef1 knockout pass.
#
# CHANGED 2026-10-17: panel_report.tsv gains `perturb_class` (OE / KD / mixed) and `dimensions` (the entry's test_dimensions, with dimension 3
# removed from the single-arm units) so Module 12 can organise results by the supervisor's 7 dimensions with OE and KD kept apart. Both are
# derived from keys that are already in the config (`direction`, `test_dimensions`), and they go into PANEL_META only, NOT into the loader
# entry JSON, so no panel_validate_one job is invalidated. New rule `panel_coverage` writes coverage_table.tsv / coverage_summary.tsv:
# every xlsx unit once, with its LINGER status and reason (resources/xlsx_test_units.tsv).
#
# Outputs:  module9c_panel/{panel_id}_score.tsv   module9c_panel/panel_report.tsv   module9c_panel/coverage_table.tsv
# Targets:  snakemake module9c_all
# =============================================================================

import copy as _copy

MODULE9C_DIR = f"{SCRATCH}/module9c_panel"
PANEL_CFG = _resolve(_copy.deepcopy(config.get("generalization_test_panel", [])))
_PANEL_DROP = {"notes", "scoring", "scoring_arms", "combinatorial_arms", "eligibility", "paired_baseline_entry",
               "perturbation_mode", "test_dimensions", "exploratory"}

PANEL_ENTRIES = {}      # panel_id -> loader entry (with sample_id / check / paired_baseline_entry)
PANEL_SPECS = {}        # panel_id -> check spec
PANEL_META = {}         # panel_id -> dict for the report
PANEL_DESCRIPTIVE = []  # rows for entries that are not perturbable

_ARM_CLASS = {"ctnnb1_activation_only": "OE", "notch1_ko_only": "KD"}   # single arms; the double arms are "mixed" (KD + OE)
_OE_WORDS = {"OE", "ACTIVATION", "CRISPRA", "OVEREXPRESSION"}


def _perturb_class(direction):
    ds = direction if isinstance(direction, (list, tuple)) else [direction]
    cls = {"OE" if str(x).upper() in _OE_WORDS else "KD" for x in ds}
    return cls.pop() if len(cls) == 1 else "mixed"


def _dims(entry, arm=None):
    dims = [str(x) for x in entry.get("test_dimensions", [])]
    if arm in _ARM_CLASS:            # a single-arm unit cannot address the combination dimension
        dims = [x for x in dims if x != "3"]
    return ",".join(dims)


def _panel_entry(base, sample_id, regex=None, baseline=None):
    e = {k: v for k, v in base.items() if k not in _PANEL_DROP}
    e["sample_id"] = sample_id
    e["check"] = "generalization_panel"
    if regex:
        e["sample_column_regex"] = regex
    if baseline:
        e["paired_baseline_entry"] = baseline
    return e


for _e in PANEL_CFG:
    _tid, _elig = _e["test_id"], _e.get("eligibility", "pending_check_tf_eligibility")
    if _elig == "eligible":
        _sc = _e["scoring"]
        _base = _e.get("paired_baseline_entry") or {}
        _bl = dict(_base, sample_id=_base.get("test_id", _tid + "_ref"))
        PANEL_ENTRIES[_tid] = _panel_entry(_e, _tid, baseline=_bl)
        PANEL_SPECS[_tid] = {"mode": "expression_shift", "ko_id": _sc["ko_id"], "sign": _sc["sign"],
                             "invert_pass": False, "top_k": 200}
        PANEL_META[_tid] = {"test_id": _tid, "dataset": _e["dataset"], "gene": _e["gene"], "arm": "",
                            "eligibility": _elig, "exploratory": bool(_e.get("exploratory", False)),
                            "perturb_class": _perturb_class(_e.get("direction")), "dimensions": _dims(_e)}
    elif _elig == "effectors_eligible":
        _ctrl_re = _e["combinatorial_arms"]["control"]["sample_column_regex"]
        _bl = {"sample_id": _tid + "__control", "path": _e["path"], "sample_column_regex": _ctrl_re}
        for _arm, _a in _e["scoring_arms"].items():
            _pid = f"{_tid}__{_arm}"
            PANEL_ENTRIES[_pid] = _panel_entry(_e, _pid, regex=_a["sample_column_regex"], baseline=_bl)
            _comps = [{"ko_id": c["ko_id"], "sign": c["sign"]} for c in _a["components"]]
            PANEL_SPECS[_pid] = {"mode": "expression_shift", "ko_id": _comps[0]["ko_id"], "sign": _comps[0]["sign"],
                                 "components": _comps, "invert_pass": False, "top_k": 200}
            PANEL_META[_pid] = {"test_id": _tid, "dataset": _e["dataset"], "gene": ",".join(_e["gene"]) if isinstance(_e["gene"], list) else _e["gene"],
                                "arm": _arm, "eligibility": _elig, "exploratory": bool(_e.get("exploratory", True)),
                                "perturb_class": _ARM_CLASS.get(_arm, "mixed"), "dimensions": _dims(_e, _arm)}
    else:
        PANEL_DESCRIPTIVE.append({"test_id": _tid, "dataset": _e["dataset"], "gene": _e["gene"], "arm": "",
                                  "eligibility": _elig, "exploratory": False,
                                  "perturb_class": _perturb_class(_e.get("direction")), "dimensions": _dims(_e)})

_PANEL_IDS = list(PANEL_ENTRIES.keys())
_PANEL_ID_PATTERN = "|".join(re.escape(x) for x in _PANEL_IDS) if _PANEL_IDS else "(?!)"


def _panel_kos(panel_id):
    spec = PANEL_SPECS[panel_id]
    return [c["ko_id"] for c in spec.get("components", [{"ko_id": spec["ko_id"]}])]


def _panel_inputs(wildcards):
    d = {"sanity_check": f"{SCRATCH}/module8_perturbation/_sanity_check_baseline_predicted.tsv"}
    for i, k in enumerate(_panel_kos(wildcards.panel_id)):
        d[f"ko_pred_{i}"] = f"{SCRATCH}/module8_perturbation/{k}_predicted_expression.tsv"
    return d


rule panel_validate_one:
    wildcard_constraints:
        panel_id = _PANEL_ID_PATTERN,
    input:
        unpack(_panel_inputs),
    output:
        score = f"{MODULE9C_DIR}/{{panel_id}}_score.tsv",
        shifts = f"{MODULE9C_DIR}/{{panel_id}}_shifts.tsv.gz",
    params:
        workdir  = f"{SCRATCH}/module4_linger_init",
        grn_dir  = LINGER_CFG["grn_dir"],
        genome   = LINGER_CFG["genome"],
        module8_dir = f"{SCRATCH}/module8_perturbation",
        module7_summary = f"{SCRATCH}/module7_tf_activity/tf_activity_summary.tsv",
        sanity_tsv = f"{SCRATCH}/module8_perturbation/_sanity_check_baseline_predicted.tsv",
        entry_json = lambda wc: json.dumps(PANEL_ENTRIES[wc.panel_id]),
        check_spec_json = lambda wc: json.dumps(PANEL_SPECS[wc.panel_id]),
        baseline_entry_json = lambda wc: json.dumps(PANEL_ENTRIES[wc.panel_id]["paired_baseline_entry"]),
    log:
        f"{SCRATCH}/logs/09c_panel_{{panel_id}}.log",
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
            --baseline-entry-json {params.baseline_entry_json:q} \
            --sanity-tsv {params.sanity_tsv} \
            {VALIDATION_FLAGS} \
            --output-tsv {output.score} \
            --output-shifts-tsv {output.shifts}
        """


rule module9c_aggregate:
    input:
        expand(f"{MODULE9C_DIR}/{{panel_id}}_score.tsv", panel_id=_PANEL_IDS),
    output:
        report = f"{MODULE9C_DIR}/panel_report.tsv",
    params:
        meta_json = json.dumps(PANEL_META),
        descriptive_json = json.dumps(PANEL_DESCRIPTIVE),
    log:
        f"{SCRATCH}/logs/09c_aggregate.log",
    resources:
        runtime   = config["resources"]["module9_aggregate"]["runtime_min"],
        sge_extra = sge_extra("module9_aggregate"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        {LINGER_PYTHON} workflow/scripts/aggregate_panel.py \
            --panel-dir {MODULE9C_DIR} \
            --meta-json {params.meta_json:q} \
            --descriptive-json {params.descriptive_json:q} \
            --output-tsv {output.report}
        """


rule panel_coverage:
    # ADDED 2026-10-17: every supervisor-xlsx test unit once, with LINGER's status and the reason (see build_coverage_table.py).
    input:
        report = f"{MODULE9C_DIR}/panel_report.tsv",
        units = "resources/xlsx_test_units.tsv",
    output:
        table = f"{MODULE9C_DIR}/coverage_table.tsv",
        summary = f"{MODULE9C_DIR}/coverage_summary.tsv",
    params:
        out_dir = MODULE9C_DIR,
    log:
        f"{SCRATCH}/logs/09c_coverage.log",
    resources:
        runtime   = config["resources"]["module9_aggregate"]["runtime_min"],
        sge_extra = sge_extra("module9_aggregate"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        {LINGER_PYTHON} workflow/scripts/build_coverage_table.py \
            --units-tsv {input.units} --panel-report {input.report} --out-dir {params.out_dir}
        """


rule module9c_all:
    input:
        f"{MODULE9C_DIR}/panel_report.tsv",
        f"{MODULE9C_DIR}/coverage_table.tsv",
