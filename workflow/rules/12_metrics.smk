# =============================================================================
# 12_metrics.smk  —  Module 12: metrics, baselines, gold-standard ranking (ADDED 2026-10-15)
# =============================================================================
# Reads the *_shifts.tsv.gz written by Modules 9 / 9b / 9c, the screen CAPS (Module 11b), the frozen
# directional gold table and the MSigDB GMTs. Writes three tables + a notes file. See metrics_lib.py for definitions.
# Run: snakemake module12_metrics_all ...   (needs module9_all module9b_all module9c_all module11b_screen_caps_all first)
#
# CHANGED 2026-10-17 (supervisor decks, FIG.2-FIG.5 / SF.1):
#   compute_metrics   + pcc_centroid (Systema-style), source_set / perturb_class / dimensions columns, metrics_recall_curve.tsv
#                       (cumulative Recall@K, Fig 4A) and metrics_epsilon.tsv (GSE86204 interaction residual, exploratory)
#   compute_ird       NEW. Identity / Rejuvenation / Damage indices per (ko_id, scope) and their real-vs-predicted validation (Fig 3)
#   ird_tf_activity_extra_population  NEW. TF-activity vector for the knockouts / OE ids that Module 11 does not score (OE ids and
#                       population-only knockouts), used ONLY for the population-scope rejuvenation index R. Writes under module12_metrics/ird,
#                       never into module11_caps, so caps_scores.tsv is unchanged.
# =============================================================================
MODULE12_DIR = f"{SCRATCH}/module12_metrics"
_M12 = config.get("metrics", {})
_GMT_ARGS = []
if _M12.get("hallmark_gmt"):
    _GMT_ARGS.append("hallmark=" + os.path.join(_M12["msigdb_dir"], _M12["hallmark_gmt"]))
if _M12.get("gobp_gmt"):
    _GMT_ARGS.append("gobp=" + os.path.join(_M12["msigdb_dir"], _M12["gobp_gmt"]))


rule compute_metrics:
    input:
        validation_report = f"{MODULE9B_DIR}/validation_report_combined.tsv",
        panel_report = f"{MODULE9C_DIR}/panel_report.tsv",
        screen_caps = f"{MODULE11_DIR}/screen_caps_all.tsv",
    output:
        expression = f"{MODULE12_DIR}/metrics_expression.tsv",
        pathway = f"{MODULE12_DIR}/metrics_pathway.tsv",
        ranking = f"{MODULE12_DIR}/metrics_ranking.tsv",
        recall_curve = f"{MODULE12_DIR}/metrics_recall_curve.tsv",
        epsilon = f"{MODULE12_DIR}/metrics_epsilon.tsv",
        notes = f"{MODULE12_DIR}/metrics_notes.txt",
    params:
        out_dir = MODULE12_DIR,
        gold = _M12.get("gold_tsv", ""),
        gmts = " ".join(_GMT_ARGS),
        n_boot = _M12.get("n_boot", 300),
        n_perm = _M12.get("n_perm", 200),
        n_null = _M12.get("n_null", 1000),
        gsea_perm = _M12.get("gsea_perm", 500),
        gsea_perm_large = _M12.get("gsea_perm_large", 100),
        min_tier = _M12.get("min_tier_tfs", 10),
        neg_ctrl = config.get("negative_control_input", {}).get("sample_id", "GSE281207_Tmie_KO"),
    log:
        f"{SCRATCH}/logs/12_metrics.log",
    resources:
        runtime   = config["resources"]["compute_metrics"]["runtime_min"],
        sge_extra = sge_extra("compute_metrics"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        {LINGER_PYTHON} workflow/scripts/compute_metrics.py \
            --module9-dir {MODULE9_DIR} --module9b-dir {MODULE9B_DIR} --module9c-dir {MODULE9C_DIR} \
            --validation-report-combined {input.validation_report} \
            --screen-caps-tsv {input.screen_caps} --gold-tsv {params.gold} \
            --gmt {params.gmts} --out-dir {params.out_dir} \
            --n-boot {params.n_boot} --n-perm {params.n_perm} --n-null {params.n_null} \
            --gsea-perm {params.gsea_perm} --gsea-perm-large {params.gsea_perm_large} \
            --min-tier-tfs {params.min_tier} \
            --negative-control-sample-id {params.neg_ctrl}
        """


# ── I / R / D ───────────────────────────────────────────────────────────────
_IRD = dict(config.get("ird", {}))
_IRD["ieg_genes"] = list(config.get("caps", {}).get("ieg_genes", []))
# knockouts / OE ids that Module 11 does not score (it keeps only caps.include_modes, and only non-population-only ids)
IRD_EXTRA_IDS = [k for k in KNOCKOUTS if k not in CAPS_KO_IDS]
IRD_TFACT = {k: (f"{MODULE11_DIR}/population/{k}_tf_activity.tsv" if k in CAPS_KO_IDS
                 else f"{MODULE12_DIR}/ird/population/{k}_tf_activity.tsv") for k in KNOCKOUTS}


rule ird_tf_activity_extra_population:
    wildcard_constraints:
        ird_extra_id = "|".join(re.escape(k) for k in IRD_EXTRA_IDS) if IRD_EXTRA_IDS else "(?!)",
    input:
        expr = f"{MODULE8_DIR}/{{ird_extra_id}}_predicted_expression.tsv",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_population_trans_regulatory.txt",
    output:
        f"{MODULE12_DIR}/ird/population/{{ird_extra_id}}_tf_activity.tsv",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        grn_dir = LINGER_CFG["grn_dir"],
        genome  = LINGER_CFG["genome"],
    log:
        f"{SCRATCH}/logs/12_ird_tf_activity_{{ird_extra_id}}.log",
    resources:
        runtime   = config["resources"]["caps_tf_activity"]["runtime_min"],
        sge_extra = sge_extra("caps_tf_activity"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        export PYTHONHASHSEED=0   # regulon() builds gene lists via list(set(...)) — hash-order dependent (same as Module 11)
        {LINGER_PYTHON} workflow/scripts/compute_tf_activity_shift.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --network "cell population" \
            --expression-tsv "{input.expr}" \
            --output-tsv "{output}"
        """


rule compute_ird:
    input:
        pop_expr = expand(f"{MODULE8_DIR}/{{k}}_predicted_expression.tsv", k=list(KNOCKOUTS.keys())),
        tf_act = list(IRD_TFACT.values()),
        baseline_tf_act = f"{MODULE11_DIR}/population/_baseline_tf_activity.tsv",
        aging_shift = f"{MODULE9_DIR}/{AGING_VECTOR_SAMPLE_ID}_aging_shift.tsv",
        validation_report = f"{MODULE9B_DIR}/validation_report_combined.tsv",
        panel_report = f"{MODULE9C_DIR}/panel_report.tsv",
        hallmark = os.path.join(_M12.get("msigdb_dir", ""), _M12.get("hallmark_gmt", "mh.all.v2023.2.Mm.symbols.gmt")),
    output:
        scores = f"{MODULE12_DIR}/ird_scores.tsv",
        validation = f"{MODULE12_DIR}/ird_validation.tsv",
        summary = f"{MODULE12_DIR}/ird_validation_summary.tsv",
        notes = f"{MODULE12_DIR}/ird_notes.txt",
    params:
        out_dir = MODULE12_DIR,
        ird_json = json.dumps(_IRD),
        tfact_json = json.dumps(IRD_TFACT),
        modes_json = json.dumps(KO_MODE_OF),
    log:
        f"{SCRATCH}/logs/12_ird.log",
    resources:
        runtime   = config["resources"]["compute_ird"]["runtime_min"],
        sge_extra = sge_extra("compute_ird"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        {LINGER_PYTHON} workflow/scripts/compute_ird.py \
            --module8-dir {MODULE8_DIR} --module9-dir {MODULE9_DIR} --module9b-dir {MODULE9B_DIR} --module9c-dir {MODULE9C_DIR} \
            --hallmark-gmt {input.hallmark} \
            --ird-config-json {params.ird_json:q} \
            --aging-shift-tsv {input.aging_shift} \
            --baseline-tf-activity-tsv {input.baseline_tf_act} \
            --tf-activity-json {params.tfact_json:q} \
            --ko-modes-json {params.modes_json:q} \
            --out-dir {params.out_dir}
        """


rule module12_metrics_all:
    input:
        f"{MODULE12_DIR}/metrics_ranking.tsv",
        f"{MODULE12_DIR}/metrics_recall_curve.tsv",
        f"{MODULE12_DIR}/ird_scores.tsv",
