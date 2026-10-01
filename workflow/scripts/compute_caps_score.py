"""
compute_caps_score.py — Module 11 step 2 (CAPS), one (ko_id, scope) row.

CAPS(TF, scope) = sign_adjusted_shift x regulatory_centrality x confidence_weight
All terms are read from files other Module 8/8b/9/9b/11-step-1 rules already
produce; nothing here recomputes a GRN or re-runs a forward pass.

REVISED 2026-10-01 (every change is a response to a defect visible in the v1
caps_scores.tsv):

sign_adjusted_shift = dot(tf_activity_shift, aging_unit)  — PROJECTION
  onto the UNIT-NORMALISED aging vector (common finite TFs). v1 used the raw
  dot product, so a scope's score scaled with the arbitrary magnitude of its
  aging vector (cell-type vectors were ~15x larger than the population one).
  The projection has the units of the TF-activity shift and is comparable
  across scopes. `cosine_to_aging` (direction only) is reported next to it,
  and `shift_l2` (the size of the knockout's own effect).
  `caps_popvec` is the same projection onto the POPULATION aging vector, as a
  sensitivity column: v1 showed the 14 cell-type aging vectors were nearly
  identical to each other (median pairwise rho 0.97), i.e. mostly shared
  regulon structure rather than cell-type biology.

NO-OP GATE: if ||tf_activity_shift||_2 < --noop-shift-l2 (default 1e-5) the
  knockout did nothing in this scope (the TF is not expressed there: v1 had ~11
  rows at |shift| ~ 1e-15 whose CAPS was float noise ranked among real
  effects). CAPS is NaN, status = "no_op". Real v1 no-ops are <= 1.6e-7 and
  real effects >= ~5e-4, so the gate sits in a clean gap.

regulatory_centrality  (--centrality percentile | raw)
  percentile (default): mean, over the knocked-out TFs, of each TF's PERCENTILE
  (0-1) of its column-sum within THIS scope's trans-regulatory matrix. v1 used
  the raw summed column-sum, whose scale differs ~1000x between the population
  matrix (colsums ~13-45) and cell-type matrices (~370-46,000) and let
  Pericyte Tbx2 rank first on matrix scale alone. `regulatory_centrality_raw`
  keeps the raw sum as a column. A TF missing from the matrix contributes 0
  and is flagged.

confidence_weight = this scope's sanity median rho from Module 9b's combined
  report, clipped to [0, 1]; 0 if the scope FAILS the Module 8/8b sanity gate
  (sanity_pass False), is missing, or non-finite.

top_contributors: the three TFs with the largest |shift_i * aging_unit_i| with
  their contribution and that TF's column-sum in this scope — what actually
  drives the score.

status: ok | no_op | scope_fails_sanity | too_few_tfs | ko_tf_no_network
  ko_tf_no_network (ADDED 2026-10-01): a knocked-out TF has no column, or an
  all-zero column, in this scope's trans-regulatory matrix (e.g. Thap11 is
  absent from every cell-type matrix; 8 TFs have zero columns in the
  population matrix). Its TF activity is undefined, so CAPS = NaN instead of
  a misleading 0.

ADDED 2026-10-01 (reporting, nothing above changes):
  gap_verdict / validated_scope : this scope's verdict on the GAP reprogramming
      check (--gap-sample-id), read from the combined validation report;
      validated_scope is True only for PASS. CAPS is NOT gated on it — it is a
      label so a reader can see which scopes the model reproduces.
  sign_adjusted_shift_noieg / CAPS_noieg : the same projection with the
      stress / immediate-early TFs (--ieg-genes) removed from both the shift
      and the aging vector (vector re-normalised on the retained TFs).
  ieg_share_aging : share of the scope aging vector's squared norm carried by
      those TFs (diagnostic: 3% population, 7-12% cell types in v2).
  n_ieg_in_top3 : how many of top_contributors are in that list.
Rank within scope (not across scopes): scopes differ in matrix, vector and
sanity quality, so cross-scope CAPS magnitudes are qualitative only.
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

p = argparse.ArgumentParser()
p.add_argument("--ko-id", required=True)
p.add_argument("--scope", required=True, help='"population" or a cell-type name')
p.add_argument("--tf-list", required=True, nargs="+", help="TF(s) this ko_id knocks out (config['perturbation']['knockouts'][ko_id])")
p.add_argument("--baseline-tf-activity-tsv", required=True)
p.add_argument("--ko-tf-activity-tsv", required=True)
p.add_argument("--aging-shift-tsv", required=True)
p.add_argument("--population-aging-shift-tsv", default=None,
               help="population aging vector, for the caps_popvec sensitivity column")
p.add_argument("--trans-regulatory-tsv", required=True)
p.add_argument("--validation-report-combined", required=True)
p.add_argument("--ieg-genes", nargs="*", default=[],
               help="stress / immediate-early TFs for the IEG-excluded variant")
p.add_argument("--gap-sample-id", default="GSE224627_GAP_reprogramming",
               help="sample_id of the GAP reprogramming check whose verdict labels the scope")
p.add_argument("--noop-shift-l2", type=float, default=1e-5)
p.add_argument("--centrality", choices=["percentile", "raw"], default="percentile")
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

print(f"--- CAPS: ko_id={args.ko_id!r} scope={args.scope!r} tfs={args.tf_list} ---")

# ------------------------------------------------------------------ shift vector
baseline_act = pd.read_csv(args.baseline_tf_activity_tsv, sep="\t", index_col=0)["tf_activity"]
ko_act = pd.read_csv(args.ko_tf_activity_tsv, sep="\t", index_col=0)["tf_activity"]
common_tfs = baseline_act.index.intersection(ko_act.index)
tf_activity_shift = ko_act.loc[common_tfs] - baseline_act.loc[common_tfs]
print(f"tf_activity_shift: {len(tf_activity_shift)} TFs (baseline vs {args.ko_id} forward-pass output)")


def _project(shift, vec_path, label, exclude=()):
    """-> (projection onto unit vec, cosine, n_common, contributions Series).
    `exclude`: TFs dropped from BOTH the shift and the vector before the vector
    is re-normalised (IEG-excluded variant)."""
    vec = pd.read_csv(vec_path, sep="\t", index_col=0)["shift"]
    if len(exclude):
        vec = vec.drop(index=[g for g in exclude if g in vec.index])
        shift = shift.drop(index=[g for g in exclude if g in shift.index])
    common = shift.index.intersection(vec.index)
    s, a = shift.loc[common], vec.loc[common]
    finite = s.notna() & a.notna()
    if (~finite).sum():
        print(f"[{label}] dropping {int((~finite).sum())}/{len(finite)} NaN TFs before projecting")
    s, a = s[finite], a[finite]
    if len(s) < 2:
        print(f"WARNING: [{label}] only {len(s)} finite TFs shared with the aging vector")
        return np.nan, np.nan, len(s), pd.Series(dtype=float)
    a_norm = float(np.linalg.norm(a.values))
    s_norm = float(np.linalg.norm(s.values))
    if a_norm == 0:
        return np.nan, np.nan, len(s), pd.Series(dtype=float)
    a_unit = a / a_norm
    contrib = s * a_unit
    proj = float(contrib.sum())
    cos = float(proj / s_norm) if s_norm > 0 else np.nan
    return proj, cos, len(s), contrib


sign_adjusted_shift, cosine_to_aging, n_common_aging_tfs, contrib = _project(
    tf_activity_shift, args.aging_shift_tsv, "scope aging vector")
sign_adjusted_shift_noieg, _, _, _ = _project(
    tf_activity_shift, args.aging_shift_tsv, "scope aging vector, IEG-excluded", exclude=args.ieg_genes)
_av = pd.read_csv(args.aging_shift_tsv, sep="\t", index_col=0)["shift"].dropna()
ieg_share_aging = (float((_av.loc[_av.index.intersection(args.ieg_genes)] ** 2).sum() / (_av ** 2).sum())
                   if len(_av) and float((_av ** 2).sum()) > 0 else np.nan)
caps_popvec_proj = np.nan
if args.population_aging_shift_tsv and os.path.exists(args.population_aging_shift_tsv):
    caps_popvec_proj, _, _, _ = _project(tf_activity_shift, args.population_aging_shift_tsv, "population aging vector")

shift_finite = tf_activity_shift.dropna()
shift_l2 = float(np.linalg.norm(shift_finite.values)) if len(shift_finite) else np.nan

# ------------------------------------------------------------------ centrality
trans_reg = pd.read_csv(args.trans_regulatory_tsv, sep="\t", index_col=0)
colsum = trans_reg.sum(axis=0)
pct = pd.Series(rankdata(colsum.values) / len(colsum), index=colsum.index)
missing = [tf for tf in args.tf_list if tf not in colsum.index]
if missing:
    print(f"WARNING: {len(missing)}/{len(args.tf_list)} of {args.ko_id}'s TFs not present in "
          f"{args.trans_regulatory_tsv}'s columns — contributing 0 centrality for: {missing}")
regulatory_centrality_raw = float(sum(colsum.get(tf, 0.0) for tf in args.tf_list))
regulatory_centrality_pct = float(np.mean([pct.get(tf, 0.0) for tf in args.tf_list]))
regulatory_centrality = regulatory_centrality_pct if args.centrality == "percentile" else regulatory_centrality_raw
ko_tf_in_network = all((tf in colsum.index) and float(colsum[tf]) > 0 for tf in args.tf_list)

# ------------------------------------------------------------------ confidence
vrc = pd.read_csv(args.validation_report_combined, sep="\t")
scope_rows = vrc[vrc["scope"] == args.scope]
gap_rows = vrc[(vrc["scope"] == args.scope) & (vrc["sample_id"] == args.gap_sample_id)]
if len(gap_rows):
    _col = "pass_fail" if "pass_fail" in gap_rows.columns else "verdict"
    gap_verdict = str(gap_rows[_col].iloc[0])
else:
    gap_verdict = "missing"
validated_scope = gap_verdict == "PASS"
raw_sanity_rho, sanity_pass = np.nan, False
if len(scope_rows) and scope_rows["sanity_rho"].notna().any():
    raw_sanity_rho = float(scope_rows["sanity_rho"].dropna().iloc[0])
    if "sanity_pass" not in scope_rows.columns:
        raise ValueError(f"{args.validation_report_combined} has no 'sanity_pass' column — it was written by the "
                         f"pre-2026-10-01 aggregate_validation_celltype.py; rerun module9b_aggregate.")
    sanity_pass = bool(scope_rows["sanity_pass"].iloc[0])
else:
    print(f"WARNING: no sanity_rho found for scope={args.scope!r} — confidence_weight=0.0")
confidence_weight = float(np.clip(raw_sanity_rho, 0.0, 1.0)) if (np.isfinite(raw_sanity_rho) and sanity_pass) else 0.0

# ------------------------------------------------------------------ status + CAPS
CAPS_noieg = np.nan
if not np.isfinite(sign_adjusted_shift):
    status, CAPS, caps_popvec = "too_few_tfs", np.nan, np.nan
elif not ko_tf_in_network:
    status, CAPS, caps_popvec = "ko_tf_no_network", np.nan, np.nan
    print(f"KO TF(s) {args.tf_list} have no usable column in {args.trans_regulatory_tsv} — TF activity undefined; CAPS = NaN")
elif not np.isfinite(shift_l2) or shift_l2 < args.noop_shift_l2:
    status, CAPS, caps_popvec = "no_op", np.nan, np.nan
    print(f"NO-OP: ||shift||_2 = {shift_l2:.3g} < {args.noop_shift_l2:g} — the knockout does nothing in this scope; CAPS = NaN")
elif not sanity_pass:
    status = "scope_fails_sanity"
    CAPS = sign_adjusted_shift * regulatory_centrality * 0.0
    CAPS_noieg = sign_adjusted_shift_noieg * regulatory_centrality * 0.0 if np.isfinite(sign_adjusted_shift_noieg) else np.nan
    caps_popvec = caps_popvec_proj * regulatory_centrality * 0.0 if np.isfinite(caps_popvec_proj) else np.nan
else:
    status = "ok"
    CAPS = sign_adjusted_shift * regulatory_centrality * confidence_weight
    CAPS_noieg = (sign_adjusted_shift_noieg * regulatory_centrality * confidence_weight
                  if np.isfinite(sign_adjusted_shift_noieg) else np.nan)
    caps_popvec = (caps_popvec_proj * regulatory_centrality * confidence_weight
                   if np.isfinite(caps_popvec_proj) else np.nan)

top = ""
if len(contrib):
    t3 = contrib.reindex(contrib.abs().sort_values(ascending=False).index).head(3)
    top = ";".join(f"{tf}(contrib={v:.3g},colsum={float(colsum.get(tf, np.nan)):.3g})" for tf, v in t3.items())

result = {
    "scope": args.scope,
    "ko_id": args.ko_id,
    "tfs": ",".join(args.tf_list),
    "sign_adjusted_shift": sign_adjusted_shift,
    "cosine_to_aging": cosine_to_aging,
    "shift_l2": shift_l2,
    "n_common_aging_tfs": n_common_aging_tfs,
    "regulatory_centrality": regulatory_centrality,
    "regulatory_centrality_raw": regulatory_centrality_raw,
    "sanity_rho": raw_sanity_rho,
    "sanity_pass": sanity_pass,
    "confidence_weight": confidence_weight,
    "CAPS": CAPS,
    "caps_popvec": caps_popvec,
    "status": status,
    "top_contributors": top,
    "gap_verdict": gap_verdict,
    "validated_scope": validated_scope,
    "sign_adjusted_shift_noieg": sign_adjusted_shift_noieg,
    "CAPS_noieg": CAPS_noieg,
    "ieg_share_aging": ieg_share_aging,
    "n_ieg_in_top3": int(sum(any(t.startswith(g + "(") for g in args.ieg_genes) for t in top.split(";") if t)),
}
pd.DataFrame([result]).to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: status={status}, sign_adjusted_shift={sign_adjusted_shift}, "
      f"shift_l2={shift_l2}, centrality={regulatory_centrality} (raw {regulatory_centrality_raw}), "
      f"confidence_weight={confidence_weight}, CAPS={CAPS}")
