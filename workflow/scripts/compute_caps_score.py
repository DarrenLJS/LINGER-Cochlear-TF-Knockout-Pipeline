"""
compute_caps_score.py — Module 11 step 2 (CAPS), one (ko_id, scope) row.

CAPS(TF, celltype) = sign_adjusted_shift x regulatory_centrality x confidence_weight
— the formula approved this session. All three terms are read from files
other Module 8/8b/9/9b/11-step-1 rules already produce; nothing here
recomputes a GRN or re-runs a forward pass.

sign_adjusted_shift = dot(tf_activity_shift[common], aging_shift[common])
  over the TF intersection of:
    - tf_activity_shift = this (ko_id, scope)'s post-KO minus pre-KO
      TF-activity vector (compute_tf_activity_shift.py's two outputs,
      diffed here)
    - aging_shift = that scope's aging-vector reference shift, already
      written by validate_aging_vector / validate_celltype_aging_vector
      (Module 9's {AGING_VECTOR_SAMPLE_ID}_aging_shift.tsv, population or
      per-celltype)
  Positive => this knockout moves the network in the SAME direction as
  real aging; negative => the opposite. CellOracle-style perturbation/
  inner-product score (Kamimoto et al. 2023), not a new formula invented
  here.

regulatory_centrality = sum over this ko_id's knocked-out TF(s) of that
  TF's column-sum in the SAME trans-regulatory file regulon() itself reads
  for this scope (cell_population_trans_regulatory.txt, or Module 6's
  cell_type_specific_trans_regulatory_{celltype}.txt) — literally the
  `colsum` term regulon() computes internally (TF_activity.py line ~109),
  read directly here rather than reimplemented. A TF missing from that
  scope's columns contributes 0 and is flagged, not silently dropped from
  the ko_id's TF list.

confidence_weight = this scope's own sanity_rho, read directly from Module
  9b's validation_report_combined.tsv (one shared value per scope, by
  construction — first non-null row is used), clipped to [0, 1]. A raw
  Spearman rho could in principle be negative or, for a scope Module 9b
  has no rows for, missing entirely — both are treated as 0.0 (zero
  confidence), not silently skipped, so a CAPS row still gets written
  rather than the whole (ko_id, scope) pair vanishing from the report.
"""
import argparse
import os

import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--ko-id", required=True)
p.add_argument("--scope", required=True, help='"population" or a cell-type name')
p.add_argument("--tf-list", required=True, nargs="+", help="TF(s) this ko_id knocks out (config['perturbation']['knockouts'][ko_id])")
p.add_argument("--baseline-tf-activity-tsv", required=True)
p.add_argument("--ko-tf-activity-tsv", required=True)
p.add_argument("--aging-shift-tsv", required=True)
p.add_argument("--trans-regulatory-tsv", required=True)
p.add_argument("--validation-report-combined", required=True)
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

print(f"--- CAPS: ko_id={args.ko_id!r} scope={args.scope!r} tfs={args.tf_list} ---")

# --- sign_adjusted_shift ---
baseline_act = pd.read_csv(args.baseline_tf_activity_tsv, sep="\t", index_col=0)["tf_activity"]
ko_act = pd.read_csv(args.ko_tf_activity_tsv, sep="\t", index_col=0)["tf_activity"]
common_tfs = baseline_act.index.intersection(ko_act.index)
tf_activity_shift = ko_act.loc[common_tfs] - baseline_act.loc[common_tfs]
print(f"tf_activity_shift: {len(tf_activity_shift)} TFs (baseline vs {args.ko_id} forward-pass output)")

aging_shift = pd.read_csv(args.aging_shift_tsv, sep="\t", index_col=0)["shift"]
common = tf_activity_shift.index.intersection(aging_shift.index)
shift_c, aging_c = tf_activity_shift.loc[common], aging_shift.loc[common]
finite = shift_c.notna() & aging_c.notna()
n_dropped = int((~finite).sum())
if n_dropped:
    print(f"sign_adjusted_shift: dropping {n_dropped}/{len(finite)} NaN TFs before dot product")
shift_c, aging_c = shift_c[finite], aging_c[finite]
if len(shift_c) < 2:
    print(f"WARNING: only {len(shift_c)} finite TFs shared between tf_activity_shift and the "
          f"aging vector — sign_adjusted_shift is unreliable, computing anyway")
sign_adjusted_shift = float(np.dot(shift_c.values, aging_c.values)) if len(shift_c) else np.nan
n_common_aging_tfs = len(shift_c)

# --- regulatory_centrality ---
trans_reg = pd.read_csv(args.trans_regulatory_tsv, sep="\t", index_col=0)
colsum = trans_reg.sum(axis=0)
missing = [tf for tf in args.tf_list if tf not in colsum.index]
if missing:
    print(f"WARNING: {len(missing)}/{len(args.tf_list)} of {args.ko_id}'s TFs not present in "
          f"{args.trans_regulatory_tsv}'s columns — contributing 0 centrality for: {missing}")
regulatory_centrality = float(sum(colsum.get(tf, 0.0) for tf in args.tf_list))

# --- confidence_weight ---
vrc = pd.read_csv(args.validation_report_combined, sep="\t")
scope_rows = vrc[vrc["scope"] == args.scope]
if len(scope_rows) == 0 or scope_rows["sanity_rho"].isna().all():
    print(f"WARNING: no sanity_rho found for scope={args.scope!r} in "
          f"{args.validation_report_combined} — using confidence_weight=0.0")
    raw_sanity_rho = np.nan
    confidence_weight = 0.0
else:
    raw_sanity_rho = float(scope_rows["sanity_rho"].dropna().iloc[0])
    confidence_weight = float(np.clip(raw_sanity_rho, 0.0, 1.0))
    if raw_sanity_rho < 0.0 or raw_sanity_rho > 1.0:
        print(f"NOTE: scope={args.scope!r} sanity_rho={raw_sanity_rho:.4f} clipped to "
              f"{confidence_weight:.4f} for confidence_weight")

CAPS = sign_adjusted_shift * regulatory_centrality * confidence_weight

result = {
    "scope": args.scope,
    "ko_id": args.ko_id,
    "tfs": ",".join(args.tf_list),
    "sign_adjusted_shift": sign_adjusted_shift,
    "n_common_aging_tfs": n_common_aging_tfs,
    "regulatory_centrality": regulatory_centrality,
    "sanity_rho": raw_sanity_rho,
    "confidence_weight": confidence_weight,
    "CAPS": CAPS,
}
pd.DataFrame([result]).to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: sign_adjusted_shift={sign_adjusted_shift}, "
      f"regulatory_centrality={regulatory_centrality}, confidence_weight={confidence_weight}, "
      f"CAPS={CAPS}")
