"""
aggregate_validation_celltype.py — Module 9b step 2, combine per-cell-type
validation scores into a report, and merge with Module 9's own population
report into one combined, EQUALLY-SCHEMA'D table.

Pure pandas. Population and every cell type are concatenated under IDENTICAL
columns with a `scope` field ("population", or a cell-type name); neither is
primary. Module 9's own file is only read, never modified.

REVISED 2026-10-01
  * pass_fail comes from the per-row `verdict` written by validate_held_out.py
    (same rule at population and cell-type level; see that script and
    aggregate_validation.py). A scope that FAILS the Module 8/8b sanity gate
    (median rho < --gate-min-median-rho or frac(rho>0.3) < --gate-min-frac-gt03)
    gets pass_fail = "SCOPE_FAILS_SANITY" on its scored rows (verdict is kept
    in its own column) and, downstream, CAPS confidence_weight 0.
  * sanity columns: sanity_rho (median), sanity_frac_gt03, sanity_n_samples,
    sanity_pass — recomputed from files already on disk (each scope's
    _sanity_check_baseline_predicted.tsv vs its real TG pseudobulk).
  * aging-vector diagnostics (the aging_vector_support rows are descriptive):
      structure_rho     = median Spearman rho between this scope's aging
                          vector and every OTHER cell-type scope's vector. v1
                          cell-type vectors were near-identical (median 0.97):
                          a high value here means the vector is mostly shared
                          regulon structure, not cell-type biology.
      aging_vec_pop_rho = Spearman rho with the population aging vector.
      aging_vec_scale   = median |shift| of this scope's vector.
    These are repeated on every row of the scope.

REVISED 2026-10-03
  * gap_minus_negctl_rho (combined report): per scope, Spearman rho of the GAP
    reprogramming check minus rho of the negative control, repeated on every row
    of the scope. A descriptive discrimination measure (does the same predicted
    shift track a real reprogramming contrast more than a no-reprogramming one?);
    it has no CI because the two contrasts use different gene sets.
  * verdict_basis now also takes the value bootstrap_raw_rho (negative control
    judged on the raw rho; see validate_held_out.py).

CAVEAT worth restating in every report built from this file's output: every
held-out/negative-control/aging dataset behind these numbers is BULK or
single-cell-collapsed RNA-seq without per-cell-type ground truth — a cell-type
scope changes which regulatory model / pre-knockout baseline the SAME real
contrast is scored against, not which cells are used.
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
p.add_argument("--gate-min-median-rho", type=float, default=0.5)
p.add_argument("--gate-min-frac-gt03", type=float, default=0.7)
p.add_argument("--gap-check", default="atoh1_gfi1_pou4f3_overexpression",
               help="positive reprogramming check used for gap_minus_negctl_rho")
p.add_argument("--negctl-check", default="negative_control_must_not_reprogram",
               help="negative control used for gap_minus_negctl_rho")
p.add_argument("--output-celltype-tsv", required=True)
p.add_argument("--output-combined-tsv", required=True)
args = p.parse_args()


def _sanity_stats(pred_path, target_path):
    """(median rho, frac rho>0.3, n_samples) — same per-gene Spearman as
    linger_perturbation.py's --sanity-check, recomputed from the two files it
    already wrote. NaNs if either file is missing/unreadable."""
    nan = (np.nan, np.nan, np.nan)
    if not (os.path.exists(pred_path) and os.path.exists(target_path)):
        return nan
    try:
        pred_df = pd.read_csv(pred_path, sep="\t", index_col=0)
        target = pd.read_csv(target_path, sep=",", header=0, index_col=0)
    except Exception as e:
        print(f"_sanity_stats: could not read {pred_path} / {target_path}: {e}")
        return nan
    common = [g for g in pred_df.index if g in target.index]
    rhos = []
    for g in common:
        rho, _ = spearmanr(pred_df.loc[g], target.loc[g])
        if not np.isnan(rho):
            rhos.append(rho)
    if not rhos:
        return nan
    rhos = np.array(rhos)
    return float(np.median(rhos)), float((rhos > 0.3).mean()), int(pred_df.shape[1])


def _gate(stats):
    med, frac, _ = stats
    return bool(np.isfinite(med) and med >= args.gate_min_median_rho
                and np.isfinite(frac) and frac >= args.gate_min_frac_gt03)


def _pass_fail(row):
    # scored verdicts are overridden when the scope itself fails the sanity gate
    if row["verdict"] in ("PASS", "FAIL", "INCONCLUSIVE") and not row["sanity_pass"]:
        return "SCOPE_FAILS_SANITY"
    return row["verdict"]


def _read_vec(path):
    return pd.read_csv(path, sep="\t", index_col=0)["shift"] if os.path.exists(path) else None


# --- cell-type rows: one sanity_rho per cell type, gathered once ---
ct_rows = []
sanity_rho_cache = {}
for path in sorted(glob.glob(os.path.join(args.module9b_dir, "*", "*_score.tsv"))):
    celltype = os.path.basename(os.path.dirname(path))
    row = pd.read_csv(path, sep="\t")
    row["scope"] = celltype
    ct_rows.append(row)
    if celltype not in sanity_rho_cache:
        sanity_rho_cache[celltype] = _sanity_stats(
            os.path.join(args.module8b_dir, celltype, "_sanity_check_baseline_predicted.tsv"),
            os.path.join(args.module8b_dir, celltype, f"TG_pseudobulk_{celltype}.tsv"),
        )

# --- aging-vector diagnostics (see module docstring) -------------------------
pop_dir = os.path.dirname(os.path.abspath(args.population_report))
pop_vec_files = sorted(glob.glob(os.path.join(pop_dir, "*_aging_shift.tsv")))
pop_vec = _read_vec(pop_vec_files[0]) if pop_vec_files else None
ct_vecs = {}
for ct in sorted(sanity_rho_cache):
    fs = sorted(glob.glob(os.path.join(args.module9b_dir, ct, "*_aging_shift.tsv")))
    v = _read_vec(fs[0]) if fs else None
    if v is not None:
        ct_vecs[ct] = v
AGING_DIAG = {}
for ct, v in ct_vecs.items():
    others = [ov for oc, ov in ct_vecs.items() if oc != ct]
    rhos = []
    for ov in others:
        c = v.index.intersection(ov.index)
        c = [t for t in c if np.isfinite(v[t]) and np.isfinite(ov[t])]
        if len(c) >= 10:
            rhos.append(spearmanr(v[c], ov[c])[0])
    pr = np.nan
    if pop_vec is not None:
        c = [t for t in v.index.intersection(pop_vec.index) if np.isfinite(v[t]) and np.isfinite(pop_vec[t])]
        if len(c) >= 10:
            pr = spearmanr(v[c], pop_vec[c])[0]
    AGING_DIAG[ct] = (float(np.median(rhos)) if rhos else np.nan, pr, float(v.abs().median()))
if ct_vecs:
    print(f"aging-vector diagnostics: {len(ct_vecs)} cell-type vectors; "
          f"median structure_rho = {np.nanmedian([d[0] for d in AGING_DIAG.values()]):.2f} "
          f"(near 1 => the vectors are mostly shared regulon structure)")

if ct_rows:
    celltype_report = pd.concat(ct_rows, ignore_index=True)
    celltype_report["sanity_rho"] = celltype_report["scope"].map(lambda c: sanity_rho_cache[c][0])
    celltype_report["sanity_frac_gt03"] = celltype_report["scope"].map(lambda c: sanity_rho_cache[c][1])
    celltype_report["sanity_n_samples"] = celltype_report["scope"].map(lambda c: sanity_rho_cache[c][2])
    celltype_report["sanity_pass"] = celltype_report["scope"].map(lambda c: _gate(sanity_rho_cache[c]))
    if "verdict" not in celltype_report.columns:
        raise ValueError("cell-type score files carry no 'verdict' column — written by the pre-2026-10-01 "
                         "validate_held_out.py; delete module9b_validation_celltype/ and rerun module9b_aggregate.")
    celltype_report["pass_fail"] = celltype_report.apply(_pass_fail, axis=1)
    celltype_report["structure_rho"] = celltype_report["scope"].map(lambda c: AGING_DIAG.get(c, (np.nan,) * 3)[0])
    celltype_report["aging_vec_pop_rho"] = celltype_report["scope"].map(lambda c: AGING_DIAG.get(c, (np.nan,) * 3)[1])
    celltype_report["aging_vec_scale"] = celltype_report["scope"].map(lambda c: AGING_DIAG.get(c, (np.nan,) * 3)[2])
else:
    celltype_report = pd.DataFrame(columns=[
        "sample_id", "check", "mode", "ko_id", "sign", "invert_pass",
        "n_genes_or_tfs", "spearman_rho", "spearman_pval", "auroc", "aupr",
        "scope", "sanity_rho", "sanity_frac_gt03", "sanity_n_samples", "sanity_pass", "verdict", "pass_fail",
    ])
    print("WARNING: no Module 9b per-cell-type score files found — writing empty reports")

celltype_report.to_csv(args.output_celltype_tsv, sep="\t", index=False)
print(f"Wrote {args.output_celltype_tsv}: {len(celltype_report)} rows, "
      f"{celltype_report['scope'].nunique() if len(celltype_report) else 0} cell types")

# --- population rows: Module 9's own report + the same sanity columns / rule ---
population_report = pd.read_csv(args.population_report, sep="\t")
population_report["scope"] = "population"
pop_stats = _sanity_stats(
    os.path.join(args.module8_dir, "_sanity_check_baseline_predicted.tsv"),
    os.path.join(args.module4_data, "TG_pseudobulk.tsv"),
)
population_report["sanity_rho"], population_report["sanity_frac_gt03"], population_report["sanity_n_samples"] = pop_stats
population_report["sanity_pass"] = _gate(pop_stats)
if "verdict" not in population_report.columns:
    raise ValueError("population report has no 'verdict' column — rerun module9_all with the revised scripts.")
population_report["pass_fail"] = population_report.apply(_pass_fail, axis=1)
population_report["structure_rho"] = np.nan
population_report["aging_vec_pop_rho"] = 1.0 if pop_vec is not None else np.nan
population_report["aging_vec_scale"] = float(pop_vec.abs().median()) if pop_vec is not None else np.nan

combined_cols = [
    "scope", "sample_id", "check", "mode", "ko_id", "sign", "invert_pass",
    "n_genes_or_tfs", "spearman_rho", "spearman_pval", "auroc", "aupr",
    "level_null_rho", "excess_rho", "excess_ci_lo", "excess_ci_hi", "rho_ci_lo", "rho_ci_hi",
    "n_real_target", "n_real_reference", "pred_shift_sd",
    "sanity_rho", "sanity_frac_gt03", "sanity_n_samples", "sanity_pass",
    "structure_rho", "aging_vec_pop_rho", "aging_vec_scale",
    "null_n", "null_median_rho", "null_p95_rho", "null_pctile", "verdict_basis",
    "verdict", "pass_fail", "gap_minus_negctl_rho",
]
combined = pd.concat(
    [population_report.reindex(columns=combined_cols), celltype_report.reindex(columns=combined_cols)],
    ignore_index=True,
)
_rho = combined.pivot_table(index="scope", columns="check", values="spearman_rho", aggfunc="first")
for _c in (args.gap_check, args.negctl_check):
    if _c not in _rho.columns:
        raise ValueError(f"check {_c!r} not found in the combined report (checks present: {list(_rho.columns)}) — "
                         f"gap_minus_negctl_rho needs both --gap-check and --negctl-check")
combined["gap_minus_negctl_rho"] = combined["scope"].map(_rho[args.gap_check] - _rho[args.negctl_check])
combined.to_csv(args.output_combined_tsv, sep="\t", index=False)

print(f"Wrote {args.output_combined_tsv}: {len(combined)} rows "
      f"({combined['scope'].nunique()} scopes: population + "
      f"{combined['scope'].nunique() - 1} cell types)")
print(combined[["scope", "sample_id", "check", "spearman_rho", "excess_rho", "sanity_rho", "sanity_pass", "pass_fail"]].to_string(index=False))
