"""
screen_caps.py — project the Module 8c all-TF screen onto CAPS (Stage 5, ADDED 2026-10-15).

WHY: Module 11 gives CAPS for the 6 configured knockouts. The screen already holds the single-TF knockout of every one
of the 586 TFs in each screened scope, as GENE-space shifts (mean over samples of perturbed - baseline prediction;
chunk_NNN_shift.tsv.gz, columns "<mode>__<TF>" plus "__baseline_mean__"). To score all of them against the aging vector
without re-running 586 forward passes through LingerGRN's regulon(), the TF-activity step is emulated here in numpy.

THE EMULATION (and what it is not). regulon() is, per sample column r = RNA / sum(RNA):
      activity_tf = (T^T r - E^T r) / (E^T r),   E = outer(row, colsum) / sum(row),
      row = rowsum(T) + mean(rowsum(T)) * 1e-6,  colsum = colsum(T) + mean(colsum(T)) * 1e-6,  T = trans-regulatory matrix
  after a quantile normalisation ACROSS the sample columns, and Module 11 averages the activity over the ~485 columns.
  Here the baseline and each variant are single mean-expression vectors b and b + s_v, so there is no cross-sample
  quantile normalisation (a one-column quantile normalisation is the identity) and the activity of the MEAN is used
  instead of the MEAN of the activities. regulon() is non-linear (division by sum(r) and by E^T r), so this is an
  APPROXIMATION of Module 11's TF-activity shift. It is checked, not assumed: for every TF that is also a Module 11
  knockout (--ko-map-json) the script writes a consistency table comparing the screen-derived TF-activity shift and
  CAPS with Module 11's own (Spearman over TFs, sign agreement, ratio of CAPS).
      TF shift_v = activity(b + s_v) - activity(b)
  Everything after that is Module 11 unchanged: projection onto the unit scope aging vector (common finite TFs),
  no-op gate (||shift||_2 < --noop-shift-l2), regulatory centrality = percentile of the TF's column-sum in the scope
  matrix, confidence_weight = the scope's sanity rho (0 if the sanity gate fails), status ok | no_op |
  scope_fails_sanity | too_few_tfs | ko_tf_no_network, rank within scope.

OUTPUT  screen_caps_<scope>.tsv    one row per single-TF variant of --modes (default ko)
        screen_caps_<scope>_consistency.tsv   per Module 11 knockout present in the screen
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

p = argparse.ArgumentParser()
p.add_argument("--scope", required=True)
p.add_argument("--screen-dir", required=True)
p.add_argument("--trans-regulatory-tsv", required=True)
p.add_argument("--aging-shift-tsv", required=True)
p.add_argument("--validation-report-combined", required=True)
p.add_argument("--modes", nargs="+", default=["ko"])
p.add_argument("--noop-shift-l2", type=float, default=1e-5)
p.add_argument("--centrality", choices=["percentile", "raw"], default="percentile")
p.add_argument("--ieg-genes", nargs="*", default=[])
p.add_argument("--gap-sample-id", default="GSE224627_GAP_reprogramming")
p.add_argument("--ko-map-json", default="{}", help='{"Atoh1": "atoh1_ko", ...}: single-TF Module 11 knockouts to check against')
p.add_argument("--module11-scope-dir", default=None, help="module11_caps/<scope> (holds _baseline_tf_activity.tsv, <ko>_tf_activity.tsv, <ko>_caps.tsv)")
p.add_argument("--output-tsv", required=True)
p.add_argument("--output-consistency-tsv", required=True)
args = p.parse_args()

print(f"--- screen CAPS: scope={args.scope!r} modes={args.modes} ---")
T = pd.read_csv(args.trans_regulatory_tsv, sep="\t", index_col=0)
colsum_raw = T.sum(axis=0)
status = pd.read_csv(os.path.join(args.screen_dir, "screen_status.tsv"), sep="\t")
chunks = sorted(glob.glob(os.path.join(args.screen_dir, "chunk_*_shift.tsv.gz")))
if not chunks:
    raise ValueError(f"{args.screen_dir}: no chunk_*_shift.tsv.gz files")

# ---- gene-space shifts: genes x variants ------------------------------------------------------------------------------
parts, base = [], None
for f in chunks:
    d = pd.read_csv(f, sep="\t", index_col=0)
    if "__baseline_mean__" in d.columns:
        base = d["__baseline_mean__"] if base is None else base
        d = d.drop(columns="__baseline_mean__")
    parts.append(d)
S = pd.concat(parts, axis=1)
if base is None:
    raise ValueError("no __baseline_mean__ column in the screen chunks")
want = status[status["mode"].isin(args.modes)]
ids = [v for v in want["variant"] if v in S.columns]
if len(ids) < len(want):
    print(f"WARNING: {len(want) - len(ids)} variants in screen_status.tsv have no shift column")
S = S[ids]

genes = [g for g in T.index if g in S.index]
if len(genes) < 1000:
    raise ValueError(f"only {len(genes)} genes shared between the trans-regulatory matrix and the screen")
T = T.loc[genes]
S = S.loc[genes].fillna(0.0)
b = base.loc[genes].astype(float)
tf_names = list(T.columns)
Tv = T.values.astype(float)
row = Tv.sum(axis=1)
row = row + row.mean() * 1e-6
col = Tv.sum(axis=0)
col = col + col.mean() * 1e-6
row_over_total = row / row.sum()


def activity(R):
    """R: genes x n columns of expression -> TFs x n activity, regulon() maths on single columns."""
    R = np.asarray(R, dtype=float)
    r = R / R.sum(axis=0, keepdims=True)
    num = Tv.T @ r
    den = (row_over_total @ r)[None, :] * col[:, None]       # E^T r  = colsum * (row . r) / sum(row)
    return (num - den) / den


act_base = activity(b.values[:, None])[:, 0]
act_var = activity(b.values[:, None] + S.values)
tf_shift = pd.DataFrame(act_var - act_base[:, None], index=tf_names, columns=S.columns)
print(f"emulated TF-activity shifts: {tf_shift.shape[0]} TFs x {tf_shift.shape[1]} variants")

# ---- projection onto the aging vector -----------------------------------------------------------------------------------
vec_all = pd.read_csv(args.aging_shift_tsv, sep="\t", index_col=0)["shift"]


def project(shift, exclude=()):
    vec = vec_all.drop(index=[g for g in exclude if g in vec_all.index])
    sh = shift.drop(index=[g for g in exclude if g in shift.index])
    common = sh.index.intersection(vec.index)
    s_, a_ = sh.loc[common], vec.loc[common]
    ok = s_.notna() & a_.notna()
    s_, a_ = s_[ok], a_[ok]
    if len(s_) < 2 or np.linalg.norm(a_.values) == 0:
        return np.nan, np.nan, len(s_)
    au = a_ / np.linalg.norm(a_.values)
    proj = float((s_ * au).sum())
    sn = float(np.linalg.norm(s_.values))
    return proj, (proj / sn if sn > 0 else np.nan), len(s_)


# ---- confidence / centrality (same reading as compute_caps_score.py) -----------------------------------------------------
vrc = pd.read_csv(args.validation_report_combined, sep="\t")
scope_rows = vrc[vrc["scope"] == args.scope]
gap_rows = scope_rows[scope_rows["sample_id"] == args.gap_sample_id]
gcol = "pass_fail" if "pass_fail" in vrc.columns else "verdict"
gap_verdict = str(gap_rows[gcol].iloc[0]) if len(gap_rows) else "missing"
raw_rho, sanity_pass = np.nan, False
if len(scope_rows) and scope_rows["sanity_rho"].notna().any():
    raw_rho = float(scope_rows["sanity_rho"].dropna().iloc[0])
    sanity_pass = bool(scope_rows["sanity_pass"].iloc[0])
conf = float(np.clip(raw_rho, 0.0, 1.0)) if (np.isfinite(raw_rho) and sanity_pass) else 0.0
pct = pd.Series(rankdata(colsum_raw.values) / len(colsum_raw), index=colsum_raw.index)

rows = []
for v in tf_shift.columns:
    mode, tf = v.split("__", 1)
    sh = tf_shift[v]
    proj, cos, ncommon = project(sh)
    proj_ne, _, _ = project(sh, exclude=args.ieg_genes)
    l2 = float(np.linalg.norm(sh.dropna().values)) if sh.notna().any() else np.nan
    in_net = bool(tf in colsum_raw.index and float(colsum_raw[tf]) > 0)
    central = float(pct.get(tf, 0.0)) if args.centrality == "percentile" else float(colsum_raw.get(tf, 0.0))
    if not np.isfinite(proj):
        st, caps, caps_ne = "too_few_tfs", np.nan, np.nan
    elif not in_net:
        st, caps, caps_ne = "ko_tf_no_network", np.nan, np.nan
    elif not np.isfinite(l2) or l2 < args.noop_shift_l2:
        st, caps, caps_ne = "no_op", np.nan, np.nan
    elif not sanity_pass:
        st, caps, caps_ne = "scope_fails_sanity", 0.0, 0.0
    else:
        st, caps = "ok", proj * central * conf
        caps_ne = proj_ne * central * conf if np.isfinite(proj_ne) else np.nan
    rows.append({"scope": args.scope, "mode": mode, "tf": tf, "variant": v, "status": st, "CAPS": caps, "CAPS_noieg": caps_ne,
                 "sign_adjusted_shift": proj, "cosine_to_aging": cos, "shift_l2": l2, "regulatory_centrality": central,
                 "confidence_weight": conf, "sanity_rho": raw_rho, "sanity_pass": sanity_pass, "gap_verdict": gap_verdict,
                 "n_common_aging_tfs": ncommon})
out = pd.DataFrame(rows)
out["rank_in_scope"] = np.nan
for m, g in out[out["status"] == "ok"].groupby("mode"):
    out.loc[g.index, "rank_in_scope"] = g["CAPS"].rank(ascending=False, method="min")
out = out.sort_values(["mode", "rank_in_scope", "tf"], na_position="last").reset_index(drop=True)
out.to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: {len(out)} rows; status counts {out['status'].value_counts().to_dict()}")

# ---- consistency with Module 11 ------------------------------------------------------------------------------------------
ko_map = json.loads(args.ko_map_json)
cons = []
m11 = args.module11_scope_dir
if m11 and os.path.exists(os.path.join(m11, "_baseline_tf_activity.tsv")):
    base_ref = pd.read_csv(os.path.join(m11, "_baseline_tf_activity.tsv"), sep="\t", index_col=0)["tf_activity"]
    for tf, ko_id in ko_map.items():
        v = f"ko__{tf}"
        f_act = os.path.join(m11, f"{ko_id}_tf_activity.tsv")
        if v not in tf_shift.columns or not os.path.exists(f_act):
            continue
        ref_shift = pd.read_csv(f_act, sep="\t", index_col=0)["tf_activity"] - base_ref
        common = ref_shift.index.intersection(tf_shift.index)
        a_, b_ = ref_shift.loc[common], tf_shift.loc[common, v]
        ok = a_.notna() & b_.notna()
        rho = float(spearmanr(a_[ok], b_[ok])[0]) if ok.sum() > 10 else np.nan
        r_caps = np.nan
        f_caps = os.path.join(m11, f"{ko_id}_caps.tsv")
        if os.path.exists(f_caps):
            rc = pd.read_csv(f_caps, sep="\t")
            r_caps = float(rc["CAPS"].iloc[0]) if "CAPS" in rc.columns else np.nan
            r_status = str(rc["status"].iloc[0]) if "status" in rc.columns else ""
        else:
            r_status = ""
        s_caps = float(out.loc[out["variant"] == v, "CAPS"].iloc[0])
        s_status = str(out.loc[out["variant"] == v, "status"].iloc[0])
        cons.append({"scope": args.scope, "tf": tf, "ko_id": ko_id, "tf_shift_spearman": rho,
                     "tf_shift_l2_screen": float(np.linalg.norm(b_[ok].values)),
                     "tf_shift_l2_module11": float(np.linalg.norm(a_[ok].values)),
                     "CAPS_screen": s_caps, "CAPS_module11": r_caps,
                     "sign_agrees": bool(np.sign(s_caps) == np.sign(r_caps)) if np.isfinite(s_caps) and np.isfinite(r_caps) else np.nan,
                     "CAPS_ratio_screen_over_m11": (s_caps / r_caps) if np.isfinite(s_caps) and np.isfinite(r_caps) and r_caps != 0 else np.nan,
                     "status_screen": s_status, "status_module11": r_status})
else:
    print(f"NOTE: no Module 11 reference under {m11!r}; consistency table left empty")
cdf = pd.DataFrame(cons, columns=["scope", "tf", "ko_id", "tf_shift_spearman", "tf_shift_l2_screen", "tf_shift_l2_module11",
                                  "CAPS_screen", "CAPS_module11", "sign_agrees", "CAPS_ratio_screen_over_m11",
                                  "status_screen", "status_module11"])
cdf.to_csv(args.output_consistency_tsv, sep="\t", index=False)
if len(cdf):
    print("consistency with Module 11 (screen-emulated vs the real regulon() path):")
    print(cdf.to_string(index=False))
