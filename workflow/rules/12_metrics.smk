# =============================================================================
# 12_metrics.smk  —  Module 12: metrics, baselines, gold-standard ranking (ADDED 2026-10-15)
# =============================================================================
# Reads the *_shifts.tsv.gz written by Modules 9 / 9b / 9c, the screen CAPS (Module 11b), the frozen
# directional gold table and the MSigDB GMTs. Writes three tables + a notes file. See metrics_lib.py for definitions.
# Run: snakemake module12_metrics_all ...   (needs module9_all module9b_all module9c_all module11b_screen_caps_all first)
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
            --min-tier-tfs {params.min_tier}
        """


rule module12_metrics_all:
    input:
        f"{MODULE12_DIR}/metrics_ranking.tsv",
