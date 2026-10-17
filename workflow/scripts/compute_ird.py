"""
compute_ird.py — Module 12: Identity / Rejuvenation / Damage (I / R / D) indices (ADDED 2026-10-17).

WHY: the supervisor's FIG.3 places every perturbation (one point = one TF-OE/KD in one cell type) on an
identity-change x rejuvenation plane, sized by a damage score, and requires that (i) every call carries an uncertainty,
(ii) a "rejuvenation" call first passes a damage-exclusion gate. The data deck's success rule is "R up, I preserved, D not up".
The definitions below follow the three evidence lines on that slide (identity: marker / regulon score; rejuvenation: aging
module; damage: apoptosis / stress / UPR panel) and are DOCUMENTED CHOICES, not a published standard.

DEFINITIONS (all on gene-level shifts, z-scored over the TESTED genes of that shift, so indices are in SD units of that
perturbation's own shift distribution and comparable across points)
  shift      perturbation points: mean(KO or OE forward pass) - mean(model baseline) per gene, for every {ko_id} x scope found in
             module8_perturbation (population) and module8_perturbation/celltype/<scope> (cell types).
             contrasts: the real_shift and predicted_shift columns of every Module 9 / 9b / 9c *_shifts.tsv.gz.
  tested     baseline (or `level`) at or above the ird.floor_pctile percentile (default 25) of the scope's genes. Lower than the fig06 MD
             floor (50) because the identity panels are small and many marker genes are lowly expressed; the decile-matched null
             controls for expression level.
  I_signed   mean z(shift) over the hair-cell panel  -  mean z(shift) over the supporting-cell panel
             (> 0: toward hair cell; < 0: toward supporting cell). I_abs = |I_signed| is the supervisor's x-axis
             ("identity change score"). Perturbed TFs (ird.exclude_genes: Atoh1, Gfi1, Pou4f3) are removed from every panel, because a
             knockout lowers its own transcript and would otherwise score itself.
  D          mean z(shift) over the union of the damage gene sets (default Hallmark APOPTOSIS + UNFOLDED_PROTEIN_RESPONSE), minus
             the immediate-early / stress genes in caps.ieg_genes (the same exclusion as CAPS_noieg) and minus exclude_genes.
             D > 0 = damage-ward.
  R          POPULATION SCOPE ONLY, EXPLORATORY. R = - (projection of the TF-activity shift onto the unit aging vector), i.e. the
             negative of CAPS's sign_adjusted_shift: R > 0 moves TF activity AWAY from the real aging signature. Not available for
             cell-type scopes: the cell-type aging vectors are near-identical to one another (structure_rho 0.96-0.985, so they
             mostly carry shared regulon structure) and there is no aged control in the network's feature space.

UNCERTAINTY AND CALLS
  CI         95% bootstrap over the PANEL GENES (I: hair-cell and supporting-cell panels resampled separately; D: damage genes;
             R: TFs). It does NOT include replicate-to-replicate variability (the real studies have 2-5 replicates).
  null       n_null random gene sets of the same size, drawn within baseline-expression deciles of the tested genes, give a null for
             |I| and for D (for R: the aging vector shuffled across TFs). A call needs the CI to exclude 0 AND the point estimate to
             exceed the null 95th percentile. Calls are NOT corrected for multiple testing (about 1 in 20 pure-noise points will be called).
  classes    non-specific damage        D called  (this also BLOCKS any rejuvenation call: the damage-exclusion gate)
             identity reprogramming     |I| called, D not called
             identity-preserving rejuvenation   R called (R > 0, CI above 0, above the shuffled-vector null; population only), neither |I| nor D called
             no clear effect            nothing called;  no_op: the shift is identically zero
  Cell-type points can only be classified on I and D; their R column is NaN and `r_assessed` is False.

VALIDATION (ird_validation.tsv / ird_validation_summary.tsv): the SAME I and D are computed on the REAL shift and the PREDICTED shift of
every Module 9/9b/9c contrast; the summary reports sign agreement (Wilson CI) and a rank correlation across contrasts. With ~10-70
contrasts this is a coarse check, not a precise estimate.

OUTPUTS (--out-dir): ird_scores.tsv, ird_validation.tsv, ird_validation_summary.tsv, ird_notes.txt
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

p = argparse.ArgumentParser()
p.add_argument("--module8-dir", required=True)
p.add_argument("--module9-dir", required=True)
p.add_argument("--module9b-dir", required=True)
p.add_argument("--module9c-dir", default=None)
p.add_argument("--hallmark-gmt", default=None)
p.add_argument("--ird-config-json", required=True, help="JSON of config['ird'] (+ ieg_genes from config['caps'])")
p.add_argument("--aging-shift-tsv", default=None, help="population aging shift vector (index = TF, column 'shift')")
p.add_argument("--baseline-tf-activity-tsv", default=None, help="population baseline TF activity (Module 11 step 1)")
p.add_argument("--tf-activity-json", default="{}", help="JSON {ko_id: path of population {ko_id}_tf_activity.tsv}")
p.add_argument("--ko-modes-json", default="{}", help="JSON {ko_id: 'ko'|'oe'}")
p.add_argument("--out-dir", required=True)
p.add_argument("--seed", type=int, default=0)
args = p.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
CFG = json.loads(args.ird_config_json)
KO_MODES = json.loads(args.ko_modes_json)
TF_ACT = json.loads(args.tf_activity_json)
B = int(CFG.get("n_boot", 1000))
NNULL = int(CFG.get("n_null", 1000))
FLOOR = float(CFG.get("floor_pctile", 25.0))
MIN_PANEL = int(CFG.get("min_panel_genes", 3))
MIN_DMG = int(CFG.get("min_damage_genes", 10))
MIN_TESTED = int(CFG.get("min_tested_genes", 200))
EXCLUDE = set(CFG.get("exclude_genes", []))
IEG = set(CFG.get("ieg_genes", []))
HC_PANEL = list(CFG["identity_panels"]["hair_cell"])
SC_PANEL = list(CFG["identity_panels"]["supporting_cell"])
DMG_SETS = list(CFG.get("damage_sets", ["HALLMARK_APOPTOSIS", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]))
notes = []


# ------------------------------------------------------------------ helpers
def read_gmt(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 2:
                out[f[0]] = [g for g in f[2:] if g]
    return out


DMG_GENES = []
if args.hallmark_gmt and os.path.exists(args.hallmark_gmt):
    gm = read_gmt(args.hallmark_gmt)
    missing_sets = [s for s in DMG_SETS if s not in gm]
    if missing_sets:
        sys.exit(f"damage set(s) {missing_sets} not found in {args.hallmark_gmt} (available e.g. {sorted(gm)[:5]})")
    DMG_GENES = sorted({g for s in DMG_SETS for g in gm[s]} - IEG - EXCLUDE)
    notes.append(f"damage panel: union of {DMG_SETS} = {len(DMG_GENES)} genes after removing {len(IEG)} IEG/stress and {len(EXCLUDE)} perturbed-TF genes")
else:
    sys.exit(f"Hallmark GMT not found: {args.hallmark_gmt!r} — run setup_scripts/download_gold_standard.sh (it fetches mh.all.v2023.2.Mm.symbols.gmt)")
HC_PANEL = [g for g in HC_PANEL if g not in EXCLUDE]
SC_PANEL = [g for g in SC_PANEL if g not in EXCLUDE]
notes.append(f"identity panels after excluding perturbed TFs: hair_cell n={len(HC_PANEL)}, supporting_cell n={len(SC_PANEL)}")


def _zs(x):
    x = np.asarray(x, dtype=float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def _pick(genes, upper_index):
    """panel genes present in the tested gene index (case-insensitive) -> positional indices"""
    pos = []
    for g in genes:
        j = upper_index.get(g.upper())
        if j is not None:
            pos.append(j)
    return np.array(sorted(set(pos)), dtype=int)


def _null_idx(rng, deciles, panel_pos, n):
    """(n, len(panel_pos)) matrix of random gene positions drawn within the SAME expression decile as each panel gene."""
    cols = []
    by_dec = {d: np.where(deciles == d)[0] for d in np.unique(deciles)}
    for d in np.unique(deciles[panel_pos]):
        c = int((deciles[panel_pos] == d).sum())
        pool = by_dec[d]
        cols.append(pool[rng.integers(0, len(pool), size=(n, c))])
    return np.hstack(cols)


def ird_scores(shift, level, rng):
    """shift, level: pandas Series indexed by gene (aligned). -> dict of I / D scores, CIs and calls (no R)."""
    d = pd.DataFrame({"s": shift, "l": level}).replace([np.inf, -np.inf], np.nan).dropna()
    out = {"n_genes_total": int(len(d))}
    if len(d) < MIN_TESTED:
        out.update(status="too_few_genes")
        return out
    floor = np.nanpercentile(d["l"].values, FLOOR)
    d = d[d["l"] >= floor]
    out["n_genes_tested"] = int(len(d))
    if len(d) < MIN_TESTED or d["s"].std() == 0:
        out.update(status="no_op" if len(d) >= MIN_TESTED else "too_few_genes")
        return out
    z = _zs(d["s"].values)
    idx = {g.upper(): i for i, g in enumerate(d.index)}
    deciles = pd.qcut(pd.Series(d["l"].values).rank(method="first"), 10, labels=False).values
    hc, sc, dm = _pick(HC_PANEL, idx), _pick(SC_PANEL, idx), _pick(DMG_GENES, idx)
    out.update(n_hc=int(len(hc)), n_sc=int(len(sc)), n_dmg=int(len(dm)), status="ok")
    raw = d["s"].values
    if len(hc) >= MIN_PANEL and len(sc) >= MIN_PANEL:
        I = float(z[hc].mean() - z[sc].mean())
        bh = z[hc][rng.integers(0, len(hc), (B, len(hc)))].mean(1)
        bs = z[sc][rng.integers(0, len(sc), (B, len(sc)))].mean(1)
        ib = bh - bs
        nh = z[_null_idx(rng, deciles, hc, NNULL)].mean(1)
        ns = z[_null_idx(rng, deciles, sc, NNULL)].mean(1)
        null_abs = np.abs(nh - ns)
        lo, hi = np.percentile(ib, [2.5, 97.5])
        thr = float(np.percentile(null_abs, 95))
        out.update(I_signed=I, I_abs=abs(I), I_ci_lo=float(lo), I_ci_hi=float(hi), I_null_p95_abs=thr,
                   I_call=bool((lo > 0 or hi < 0) and abs(I) > thr),
                   HC_mean_z=float(z[hc].mean()), SC_mean_z=float(z[sc].mean()),
                   HC_mean_raw=float(raw[hc].mean()), SC_mean_raw=float(raw[sc].mean()),
                   direction="toward hair cell" if I > 0 else "toward supporting cell")
    else:
        out.update(I_call=False, status="identity_panel_too_small")
    if len(dm) >= MIN_DMG:
        D = float(z[dm].mean())
        db = z[dm][rng.integers(0, len(dm), (B, len(dm)))].mean(1)
        dn = z[_null_idx(rng, deciles, dm, NNULL)].mean(1)
        lo, hi = np.percentile(db, [2.5, 97.5])
        thr = float(np.percentile(dn, 95))
        out.update(D=D, D_ci_lo=float(lo), D_ci_hi=float(hi), D_null_p95=thr, D_raw=float(raw[dm].mean()),
                   D_call=bool(lo > 0 and D > thr))
    else:
        out.update(D_call=False)
        out["status"] = out["status"] if out["status"] != "ok" else "damage_panel_too_small"
    return out


def rejuvenation(ko_id, rng):
    """Population-only R = - projection of the TF-activity shift on the unit aging vector. -> (R, lo, hi, n_tfs) or NaNs."""
    nan = (np.nan, np.nan, np.nan, 0, np.nan, np.nan)
    path = TF_ACT.get(ko_id)
    if not (path and os.path.exists(path) and args.baseline_tf_activity_tsv and os.path.exists(args.baseline_tf_activity_tsv)
            and args.aging_shift_tsv and os.path.exists(args.aging_shift_tsv)):
        return nan
    base = pd.read_csv(args.baseline_tf_activity_tsv, sep="\t", index_col=0)["tf_activity"]
    ko = pd.read_csv(path, sep="\t", index_col=0)["tf_activity"]
    vec = pd.read_csv(args.aging_shift_tsv, sep="\t", index_col=0)["shift"]
    common = base.index.intersection(ko.index).intersection(vec.index)
    s = (ko.loc[common] - base.loc[common]).values.astype(float)
    a = vec.loc[common].values.astype(float)
    fin = np.isfinite(s) & np.isfinite(a)
    s, a = s[fin], a[fin]
    if len(s) < 10 or np.linalg.norm(a) == 0:
        return nan
    R = float(-(s * a).sum() / np.linalg.norm(a))
    i = rng.integers(0, len(s), (B, len(s)))
    sb, ab = s[i], a[i]
    nb = np.linalg.norm(ab, axis=1)
    nb[nb == 0] = np.nan
    Rb = -(sb * ab).sum(1) / nb
    lo, hi = np.nanpercentile(Rb, [2.5, 97.5])
    # null: the aging vector's entries shuffled across TFs (same size, same marginal distribution, no TF-specific alignment)
    perm = np.argsort(rng.random((NNULL, len(a))), axis=1)
    Rn = -(s[None, :] * a[perm]).sum(1) / np.linalg.norm(a)
    return R, float(lo), float(hi), int(len(s)), float(np.percentile(Rn, 95)), float(np.percentile(Rn, 5))


def classify(row):
    if row.get("status") == "no_op":
        return "no-op"
    if row.get("status") in (None, "too_few_genes") or (isinstance(row.get("status"), float)):
        return "not assessed"
    if row.get("D_call") is True:
        return "non-specific damage"
    if row.get("I_call") is True:
        return "identity reprogramming"
    if row.get("r_assessed") and row.get("R_call") is True:
        return "identity-preserving rejuvenation"
    return "no clear effect"


# ------------------------------------------------------------------ A. perturbation points (module 8 / 8b)
rng = np.random.default_rng(args.seed)


def _shift_from(baseline_path, ko_path):
    b = pd.read_csv(baseline_path, sep="\t", index_col=0)
    k = pd.read_csv(ko_path, sep="\t", index_col=0)
    genes = b.index.intersection(k.index)
    cols = b.columns.intersection(k.columns)
    bm, km = b.loc[genes, cols].mean(axis=1), k.loc[genes, cols].mean(axis=1)
    return km - bm, bm


points = []
scopes = [("population", args.module8_dir)]
for d in sorted(glob.glob(os.path.join(args.module8_dir, "celltype", "*"))):
    if os.path.isdir(d):
        scopes.append((os.path.basename(d), d))
for scope, d in scopes:
    base_path = os.path.join(d, "_sanity_check_baseline_predicted.tsv")
    if not os.path.exists(base_path):
        notes.append(f"scope {scope}: no _sanity_check_baseline_predicted.tsv in {d} — skipped")
        continue
    for ko_path in sorted(glob.glob(os.path.join(d, "*_predicted_expression.tsv"))):
        ko_id = re.sub(r"_predicted_expression\.tsv$", "", os.path.basename(ko_path))
        sh, lv = _shift_from(base_path, ko_path)
        row = {"scope": scope, "ko_id": ko_id, "mode": KO_MODES.get(ko_id, "ko" if ko_id.endswith("_ko") else ("oe" if ko_id.endswith("_oe") else "")),
               **ird_scores(sh, lv, rng)}
        if scope == "population":
            R, lo, hi, nt, rp95, rp05 = rejuvenation(ko_id, rng)
            row.update(R=R, R_ci_lo=lo, R_ci_hi=hi, R_null_p95=rp95, R_null_p05=rp05, n_tfs_R=nt, r_assessed=bool(np.isfinite(R)),
                       R_call=bool(np.isfinite(R) and lo > 0 and R > rp95))
        else:
            row.update(R=np.nan, R_ci_lo=np.nan, R_ci_hi=np.nan, R_null_p95=np.nan, R_null_p05=np.nan, n_tfs_R=0, r_assessed=False, R_call=False)
        row["ird_class"] = classify(row)
        points.append(row)
pts = pd.DataFrame(points)
if len(pts):
    first = ["scope", "ko_id", "mode", "status", "ird_class", "I_signed", "I_abs", "I_ci_lo", "I_ci_hi", "I_null_p95_abs", "I_call", "direction",
             "D", "D_ci_lo", "D_ci_hi", "D_null_p95", "D_call", "R", "R_ci_lo", "R_ci_hi", "R_null_p95", "R_call", "r_assessed"]
    pts = pts[[c for c in first if c in pts.columns] + [c for c in pts.columns if c not in first]]
pts.to_csv(os.path.join(args.out_dir, "ird_scores.tsv"), sep="\t", index=False)
print(f"ird_scores.tsv: {len(pts)} points ({pts['scope'].nunique() if len(pts) else 0} scopes)")
if len(pts):
    print(pts["ird_class"].value_counts().to_string())

# ------------------------------------------------------------------ B. validation on real vs predicted contrast shifts
files = []
for f in sorted(glob.glob(os.path.join(args.module9_dir, "*_shifts.tsv.gz"))):
    files.append(("population", os.path.basename(f)[: -len("_shifts.tsv.gz")], "validation", f))
for f in sorted(glob.glob(os.path.join(args.module9b_dir, "*", "*_shifts.tsv.gz"))):
    files.append((os.path.basename(os.path.dirname(f)), os.path.basename(f)[: -len("_shifts.tsv.gz")], "validation", f))
if args.module9c_dir:
    for f in sorted(glob.glob(os.path.join(args.module9c_dir, "*_shifts.tsv.gz"))):
        files.append(("population", os.path.basename(f)[: -len("_shifts.tsv.gz")], "panel", f))
vrows = []
for scope, sid, grp, f in files:
    d = pd.read_csv(f, sep="\t").dropna()
    if len(d) < MIN_TESTED:
        continue
    g = d["gene"].astype(str).values
    lv = pd.Series(d["level"].values, index=g)
    rs = ird_scores(pd.Series(d["real_shift"].values, index=g), lv, rng)
    ps = ird_scores(pd.Series(d["predicted_shift"].values, index=g), lv, rng)
    row = {"scope": scope, "contrast": sid, "group": grp}
    for pre, sc_ in (("real", rs), ("pred", ps)):
        for k in ("status", "n_genes_tested", "n_hc", "n_sc", "n_dmg", "I_signed", "I_ci_lo", "I_ci_hi", "I_call", "D", "D_ci_lo", "D_ci_hi", "D_call"):
            row[f"{pre}_{k}"] = sc_.get(k, np.nan)
    vrows.append(row)
val = pd.DataFrame(vrows)


def _wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (float(c - h), float(c + h))


def _clear(lo, hi):
    return (lo > 0) | (hi < 0)


summ = []
if len(val):
    for label, sub in (("all_contrasts", val), ("population_scope", val[val["scope"] == "population"]),
                       ("cell_type_scopes", val[val["scope"] != "population"])):
        for axis, col_ in (("I", "I_signed"), ("D", "D")):
            r_, p_ = f"real_{col_}", f"pred_{col_}"
            ok = sub.dropna(subset=[r_, p_])
            ax_ = "I" if axis == "I" else "D"
            clear = ok[_clear(ok[f"real_{ax_}_ci_lo"], ok[f"real_{ax_}_ci_hi"])]
            agree = int((np.sign(clear[r_]) == np.sign(clear[p_])).sum())
            lo, hi = _wilson(agree, len(clear))
            rho = spearmanr(ok[r_], ok[p_])[0] if len(ok) >= 5 else np.nan
            summ.append({"subset": label, "axis": axis, "n_contrasts": int(len(ok)), "n_real_clear": int(len(clear)),
                         "n_sign_agree": agree, "sign_agreement": (agree / len(clear)) if len(clear) else np.nan,
                         "sign_agree_ci_lo": lo, "sign_agree_ci_hi": hi, "spearman_rho_real_vs_pred": rho,
                         "note": "sign agreement counted only where the REAL index's CI excludes 0; chance = 0.5; gene-bootstrap CIs"})
val.to_csv(os.path.join(args.out_dir, "ird_validation.tsv"), sep="\t", index=False)
pd.DataFrame(summ).to_csv(os.path.join(args.out_dir, "ird_validation_summary.tsv"), sep="\t", index=False)
print(f"ird_validation.tsv: {len(val)} contrasts; ird_validation_summary.tsv: {len(summ)} rows")
if summ:
    print(pd.DataFrame(summ)[["subset", "axis", "n_contrasts", "n_real_clear", "sign_agreement", "spearman_rho_real_vs_pred"]].to_string(index=False))

notes.append("R: population scope only, exploratory (see module docstring); cell-type scopes have R = NaN and r_assessed = False.")
notes.append("CIs resample panel genes (or TFs for R) only; replicate-level variability is not captured.")
notes.append("I/R/D definitions are documented choices following the supervisor's FIG.3 evidence lines; panels are in config['ird'].")
open(os.path.join(args.out_dir, "ird_notes.txt"), "w").write("\n".join(notes) + "\n")
print("\n".join(notes))
