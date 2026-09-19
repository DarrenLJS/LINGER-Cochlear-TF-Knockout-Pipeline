"""
aggregate_validation_celltype.py — Module 9b step 2, combine per-cell-type
validation scores into a report, and merge with Module 9's own population
report into one combined, EQUALLY-SCHEMA'D table.

Pure pandas — same reasoning as aggregate_validation.py/aggregate_tf_activity.py
for not needing the heavier LINGER_PYTHON env, though this script is still run
under LINGER_ENV_BIN via 09b_validation_celltype.smk to keep the pipeline's
existing per-rule env convention simple (one env prefix per module).

PASS/FAIL LOGIC: identical to aggregate_validation.py's _pass_fail — copied,
not reimplemented differently, so a population row and a cell-type row for
the same check are judged by the exact same rule. See that script's
docstring for the full reasoning (rho>0 AND auroc>0.5 for positive checks;
inverted for the negative control; aging checks are "not_scored").

sanity_rho COLUMN: per the explicit instruction this session — no automatic
confidence flag, no threshold, just the raw number for a human to judge.
Recomputed directly from files already on disk (each scope's own
"_sanity_check_baseline_predicted.tsv" vs its own real TG_pseudobulk
baseline), using the exact same Spearman-per-gene-then-median calculation
linger_perturbation.py's --sanity-check mode already does — NOT re-run via
linger_perturbation.py (that would need a fresh forced Snakemake rerun just
to persist a number that's already fully recoverable from existing output
files). One sanity_rho value per scope (population, or each cell type),
repeated across that scope's rows since it's a property of the scope, not
of any individual check.

EQUAL TREATMENT, not primary/secondary: population and every cell type are
concatenated under IDENTICAL columns with a `scope` field distinguishing
them ("population", or a cell-type name) — this is the whole point of this
script's existence rather than just writing a cell-type-only report next to
Module 9's untouched validation_report.tsv. Module 9's file itself is never
modified by this script, only read.

CAVEAT worth restating in every report built from this file's output: every
held-out/negative-control/aging dataset behind spearman_rho/auroc/aupr here
is BULK RNA-seq (see validate_held_out.py's docstring) — a cell-type scope
changes which regulatory model / pre-knockout baseline the SAME bulk sample
is scored against, not which cells are used. There is no real per-cell-type
ground truth in this dataset collection.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

p = argparse.ArgumentParser()
p.add_argument("--module9b-dir", required=True, help="MODULE9B_DIR — holds {celltype}/{sample_id}_score.tsv")
p.add_argument("--population-report", required=True, help="Module 9's own validation_report.tsv (read-only)")
p.add_argument("--module8-dir", required=True, help="module8_perturbation (population sanity-check files)")
p.add_argument("--module8b-dir", required=True, help="module8_perturbation/celltype (per-celltype sanity-check files)")
p.add_argument("--module4-data", required=True, help="module4_linger_init/data (population TG_pseudobulk.tsv)")
p.add_argument("--output-celltype-tsv", required=True)
p.add_argument("--output-combined-tsv", required=True)
args = p.parse_args()


def _sanity_rho(pred_path, target_path):
    """Median Spearman rho, predicted vs real Target, over genes with a
    trained net — identical calculation to linger_perturbation.py's
    --sanity-check block, recomputed from the two files it already wrote
    rather than re-run. Returns NaN if either file is missing/unreadable
    or no genes overlap (mirrors that script's own NaN fallback)."""
    if not (os.path.exists(pred_path) and os.path.exists(target_path)):
        return np.nan
    try:
        pred_df = pd.read_csv(pred_path, sep="\t", index_col=0)
        target = pd.read_csv(target_path, sep=",", header=0, index_col=0)
    except Exception as e:
        print(f"_sanity_rho: could not read {pred_path} / {target_path}: {e}")
        return np.nan
    common = [g for g in pred_df.index if g in target.index]
    rhos = []
    for g in common:
        rho, _ = spearmanr(pred_df.loc[g], target.loc[g])
        if not np.isnan(rho):
            rhos.append(rho)
    return float(np.median(rhos)) if rhos else np.nan


def _pass_fail(row):
    if row["mode"] == "aging_tf_activity":
        return "not_scored"
    signal = (row["spearman_rho"] > 0) and (pd.notna(row["auroc"]) and row["auroc"] > 0.5)
    if row["invert_pass"]:
        return "PASS" if not signal else "FAIL"
    return "PASS" if signal else "FAIL"


# --- cell-type rows: one sanity_rho per cell type, gathered once ---
ct_rows = []
sanity_rho_cache = {}
for path in sorted(glob.glob(os.path.join(args.module9b_dir, "*", "*_score.tsv"))):
    celltype = os.path.basename(os.path.dirname(path))
    row = pd.read_csv(path, sep="\t")
    row["scope"] = celltype
    ct_rows.append(row)
    if celltype not in sanity_rho_cache:
        sanity_rho_cache[celltype] = _sanity_rho(
            os.path.join(args.module8b_dir, celltype, "_sanity_check_baseline_predicted.tsv"),
            os.path.join(args.module8b_dir, celltype, f"TG_pseudobulk_{celltype}.tsv"),
        )

if ct_rows:
    celltype_report = pd.concat(ct_rows, ignore_index=True)
    celltype_report["sanity_rho"] = celltype_report["scope"].map(sanity_rho_cache)
    celltype_report["pass_fail"] = celltype_report.apply(_pass_fail, axis=1)
else:
    celltype_report = pd.DataFrame(columns=[
        "sample_id", "check", "mode", "ko_id", "sign", "invert_pass",
        "n_genes_or_tfs", "spearman_rho", "spearman_pval", "auroc", "aupr",
        "scope", "sanity_rho", "pass_fail",
    ])
    print("WARNING: no Module 9b per-cell-type score files found — writing empty reports")

celltype_report.to_csv(args.output_celltype_tsv, sep="\t", index=False)
print(f"Wrote {args.output_celltype_tsv}: {len(celltype_report)} rows, "
      f"{celltype_report['scope'].nunique() if len(celltype_report) else 0} cell types")

# --- population rows: read Module 9's own report, unmodified, add matching columns ---
population_report = pd.read_csv(args.population_report, sep="\t")
population_report["scope"] = "population"
population_sanity_rho = _sanity_rho(
    os.path.join(args.module8_dir, "_sanity_check_baseline_predicted.tsv"),
    os.path.join(args.module4_data, "TG_pseudobulk.tsv"),
)
population_report["sanity_rho"] = population_sanity_rho
# population_report already has its own pass_fail column from aggregate_validation.py —
# kept as-is (same _pass_fail logic, just computed there) rather than recomputed here,
# so this script never second-guesses Module 9's own published numbers.

combined_cols = [
    "scope", "sample_id", "check", "mode", "ko_id", "sign", "invert_pass",
    "n_genes_or_tfs", "spearman_rho", "spearman_pval", "auroc", "aupr",
    "sanity_rho", "pass_fail",
]
combined = pd.concat(
    [population_report[combined_cols], celltype_report[combined_cols]],
    ignore_index=True,
)
combined.to_csv(args.output_combined_tsv, sep="\t", index=False)

print(f"Wrote {args.output_combined_tsv}: {len(combined)} rows "
      f"({combined['scope'].nunique()} scopes: population + "
      f"{combined['scope'].nunique() - 1} cell types)")
print(combined[["scope", "sample_id", "check", "spearman_rho", "auroc", "sanity_rho", "pass_fail"]].to_string(index=False))
