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
args = p.parse_args()

rows = []
for path in sorted(glob.glob(os.path.join(args.module11_dir, "*", "*_caps.tsv"))):
    rows.append(pd.read_csv(path, sep="\t"))

COLS = ["scope", "ko_id", "tfs", "rank_in_scope", "status", "CAPS", "sign_adjusted_shift",
        "cosine_to_aging", "shift_l2", "regulatory_centrality", "regulatory_centrality_raw",
        "confidence_weight", "sanity_rho", "sanity_pass", "caps_popvec",
        "caps_popvec_rank_in_scope", "n_common_aging_tfs", "top_contributors",
        "gap_verdict", "validated_scope", "sign_adjusted_shift_noieg", "CAPS_noieg",
        "rank_in_scope_noieg", "ieg_share_aging", "n_ieg_in_top3"]

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
    report["_pop"] = (report["scope"] != "population").astype(int)
    report = (report.sort_values(["_pop", "scope", "rank_in_scope", "ko_id"], na_position="last")
                    .drop(columns="_pop").reset_index(drop=True))
    report = report.reindex(columns=COLS)
else:
    report = pd.DataFrame(columns=COLS)
    print("WARNING: no Module 11 per-(ko_id, scope) CAPS files found — writing empty report")

report.to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: {len(report)} rows, "
      f"{report['scope'].nunique() if len(report) else 0} scopes "
      f"(population + cell types), {report['ko_id'].nunique() if len(report) else 0} ko_ids")
if len(report):
    print("status counts:", report["status"].value_counts().to_dict())
    print(report[["scope", "ko_id", "rank_in_scope", "status", "CAPS", "cosine_to_aging",
                  "regulatory_centrality", "confidence_weight"]].to_string(index=False))
