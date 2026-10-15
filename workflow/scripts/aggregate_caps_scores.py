"""
aggregate_caps_scores.py — Module 11 step 3, combine every per-(ko_id,
scope) CAPS row into one report.

EQUAL TREATMENT, same convention as aggregate_validation_celltype.py: every
row (population's own ko_id x population, and every cell type's ko_id x that
cell type) is concatenated under one identical schema with a `scope` column.

ADDED 2026-10-01: gap_verdict / validated_scope (GAP check verdict for the scope,
label only), the IEG-excluded variant (sign_adjusted_shift_noieg, CAPS_noieg,
rank_in_scope_noieg), ieg_share_aging and n_ieg_in_top3 — see
compute_caps_score.py. Status `ko_tf_no_network` rows have CAPS NaN.

ADDED 2026-10-03: `headline` (status ok AND the scope's GAP verdict is PASS) and a
second output, caps_headline_consensus.tsv, one row per ko_id over the HEADLINE
scopes only: how many validated scopes it appears in, the median CAPS / CAPS_noieg
and median within-scope rank, and how consistently the sign agrees. The
population row is carried alongside as an unvalidated reference (population_CAPS,
population_rank, population_gap_verdict), never mixed into the median. CAPS
magnitudes are not comparable across scopes, so the consensus leans on ranks and
sign agreement; median_CAPS is a descriptive summary.

REVISED 2026-10-01 — ranking is WITHIN SCOPE. Scopes differ in trans-
regulatory matrix, aging vector and sanity quality, so CAPS magnitudes are
not comparable across scopes (v1 sorted every scope together and the table was
dominated by whichever matrix had the biggest column sums). Added:
  rank_in_scope : 1 = highest CAPS among that scope's `ok` rows (NaN for
                  no_op / too_few_tfs / scope_fails_sanity rows)
  caps_popvec_rank_in_scope : same ranking by the population-vector
                  sensitivity column
Rows are sorted population first, then cell types alphabetically, then by
rank_in_scope. Rows whose status is not `ok` keep their reason in `status`.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--module11-dir", required=True, help="MODULE11_DIR — holds {scope}/{ko_id}_caps.tsv")
p.add_argument("--output-tsv", required=True)
p.add_argument("--output-consensus-tsv", required=True,
               help="per-ko_id consensus over the headline (GAP-validated) scopes")
p.add_argument("--validation-report-combined", default=None,
               help="validation_report_combined.tsv; supplies the per-scope negative-control verdict / rho")
p.add_argument("--gap-sample-id", default="GSE224627_GAP_reprogramming")
p.add_argument("--negctl-check", default="negative_control_must_not_reprogram")
args = p.parse_args()

rows = []
for path in sorted(glob.glob(os.path.join(args.module11_dir, "*", "*_caps.tsv"))):
    rows.append(pd.read_csv(path, sep="\t"))

COLS = ["scope", "ko_id", "tfs", "rank_in_scope", "status", "CAPS", "sign_adjusted_shift",
        "cosine_to_aging", "shift_l2", "regulatory_centrality", "regulatory_centrality_raw",
        "confidence_weight", "sanity_rho", "sanity_pass", "caps_popvec",
        "caps_popvec_rank_in_scope", "n_common_aging_tfs", "top_contributors",
        "gap_verdict", "validated_scope", "sign_adjusted_shift_noieg", "CAPS_noieg",
        "rank_in_scope_noieg", "ieg_share_aging", "n_ieg_in_top3", "headline",
        "rank_pct", "negctl_verdict", "gap_minus_negctl_rho",
        "ko_tf_mean_expr", "ko_tf_expressed", "noop_reason"]

CONSENSUS_COLS = ["ko_id", "tfs", "n_headline_scopes", "headline_scopes", "median_rank", "median_CAPS",
                  "median_CAPS_noieg", "majority_sign", "sign_agreement", "n_ok_scopes",
                  "population_CAPS", "population_rank", "population_gap_verdict", "median_rank_pct",
                  "n_negctl_pass", "median_gap_minus_negctl_rho"]


def _consensus(report):
    """One row per ko_id over rows with headline == True (ok status, GAP verdict PASS)."""
    out = []
    pop = report[report["scope"] == "population"].set_index("ko_id")
    for ko_id, g in report.groupby("ko_id"):
        h = g[g["headline"]]
        if h.empty:
            continue
        signs = np.sign(h["CAPS"].astype(float))
        pos, neg = int((signs > 0).sum()), int((signs < 0).sum())
        pr = pop.loc[ko_id] if ko_id in pop.index else None
        out.append({
            "ko_id": ko_id, "tfs": h["tfs"].iloc[0], "n_headline_scopes": len(h),
            "headline_scopes": ",".join(sorted(h["scope"])),
            "median_rank": float(h["rank_in_scope"].median()),
            "median_CAPS": float(h["CAPS"].median()),
            "median_CAPS_noieg": float(h["CAPS_noieg"].median()) if h["CAPS_noieg"].notna().any() else np.nan,
            "majority_sign": 1 if pos > neg else (-1 if neg > pos else 0),
            "sign_agreement": max(pos, neg) / len(h),
            "n_ok_scopes": int((g["status"] == "ok").sum()),
            "population_CAPS": float(pr["CAPS"]) if pr is not None else np.nan,
            "population_rank": float(pr["rank_in_scope"]) if pr is not None else np.nan,
            "population_gap_verdict": pr["gap_verdict"] if pr is not None else np.nan,
            "median_rank_pct": float(h["rank_pct"].median()),
            "n_negctl_pass": int((h["negctl_verdict"] == "PASS").sum()),
            "median_gap_minus_negctl_rho": float(h["gap_minus_negctl_rho"].median())
                                           if h["gap_minus_negctl_rho"].notna().any() else np.nan,
        })
    res = pd.DataFrame(out, columns=CONSENSUS_COLS)
    return res.sort_values(["median_rank", "ko_id"]).reset_index(drop=True)

if rows:
    report = pd.concat(rows, ignore_index=True)
    if "status" not in report.columns:
        raise ValueError("per-row CAPS files lack 'status' — written by the pre-2026-10-01 "
                         "compute_caps_score.py; delete module11_caps/ and rerun module11_all.")
    for _c in ("CAPS_noieg", "gap_verdict"):
        if _c not in report.columns:
            report[_c] = np.nan
    report["rank_in_scope"] = np.nan
    report["caps_popvec_rank_in_scope"] = np.nan
    report["rank_in_scope_noieg"] = np.nan
    for scope, idx in report.groupby("scope").groups.items():
        ok = report.loc[idx][report.loc[idx, "status"] == "ok"]
        report.loc[ok.index, "rank_in_scope"] = ok["CAPS"].rank(ascending=False, method="min")
        okn = ok[ok["CAPS_noieg"].notna()]
        report.loc[okn.index, "rank_in_scope_noieg"] = okn["CAPS_noieg"].rank(ascending=False, method="min")
        okp = ok[ok["caps_popvec"].notna()]
        report.loc[okp.index, "caps_popvec_rank_in_scope"] = okp["caps_popvec"].rank(ascending=False, method="min")
    # rank percentile within scope (1 = best, 0 = worst); comparable across scopes of different size
    report["rank_pct"] = np.nan
    for scope, idx in report.groupby("scope").groups.items():
        r = report.loc[idx, "rank_in_scope"]
        n = int(r.notna().sum())
        if n > 1:
            report.loc[idx, "rank_pct"] = 1.0 - (r - 1.0) / (n - 1.0)
        elif n == 1:
            report.loc[idx, "rank_pct"] = 1.0
    # negative-control verdict and GAP-minus-negctl rho per scope (label only; never gates the headline)
    report["negctl_verdict"] = pd.Series([None] * len(report), index=report.index, dtype=object)
    report["gap_minus_negctl_rho"] = np.nan
    if args.validation_report_combined and os.path.exists(args.validation_report_combined):
        vr = pd.read_csv(args.validation_report_combined, sep="\t")
        vcol = "pass_fail" if "pass_fail" in vr.columns else "verdict"
        for scope in report["scope"].unique():
            sv = vr[vr["scope"] == scope]
            nc = sv[sv["check"] == args.negctl_check] if "check" in sv.columns else sv.iloc[0:0]
            m = report["scope"] == scope
            if len(nc):
                report.loc[m, "negctl_verdict"] = str(nc[vcol].iloc[0])
            if "gap_minus_negctl_rho" in sv.columns and sv["gap_minus_negctl_rho"].notna().any():
                report.loc[m, "gap_minus_negctl_rho"] = float(sv["gap_minus_negctl_rho"].dropna().iloc[0])
    else:
        print("NOTE: no --validation-report-combined; negctl_verdict / gap_minus_negctl_rho left NaN")
    report["headline"] = (report["status"] == "ok") & (report["gap_verdict"] == "PASS")
    report["_pop"] = (report["scope"] != "population").astype(int)
    report = (report.sort_values(["_pop", "scope", "rank_in_scope", "ko_id"], na_position="last")
                    .drop(columns="_pop").reset_index(drop=True))
    report = report.reindex(columns=COLS)
else:
    report = pd.DataFrame(columns=COLS)
    print("WARNING: no Module 11 per-(ko_id, scope) CAPS files found — writing empty report")

consensus = _consensus(report) if len(report) else pd.DataFrame(columns=CONSENSUS_COLS)
consensus.to_csv(args.output_consensus_tsv, sep="\t", index=False)
print(f"Wrote {args.output_consensus_tsv}: {len(consensus)} ko_ids over "
      f"{int(report['headline'].sum()) if len(report) else 0} headline rows")
if len(report) and consensus.empty:
    print("WARNING: no scope has an `ok` CAPS row with GAP verdict PASS — the headline consensus is empty")
if len(consensus):
    print(consensus.to_string(index=False))

report.to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: {len(report)} rows, "
      f"{report['scope'].nunique() if len(report) else 0} scopes "
      f"(population + cell types), {report['ko_id'].nunique() if len(report) else 0} ko_ids")
if len(report):
    print("status counts:", report["status"].value_counts().to_dict())
    print(report[["scope", "ko_id", "rank_in_scope", "status", "CAPS", "cosine_to_aging",
                  "regulatory_centrality", "confidence_weight"]].to_string(index=False))
