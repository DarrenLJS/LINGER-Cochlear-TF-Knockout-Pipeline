"""
compute_metrics.py — Module 12: metrics, baselines and gold-standard ranking (ADDED 2026-10-15). See metrics_lib.py for the
exact definitions of every number.

INPUTS (all produced by the pipeline; missing optional inputs are skipped and said so, never silently zeroed):
  module9_validation/*_shifts.tsv.gz, module9b_validation_celltype/<scope>/*_shifts.tsv.gz, module9c_panel/*_shifts.tsv.gz
        gene-aligned real / predicted shift + level vectors written by validate_held_out.py --output-shifts-tsv
  module9b_validation_celltype/validation_report_combined.tsv     (verdicts / scored flag)
  module11_caps/screen_caps_all.tsv                                (screen CAPS, Stage 5)
  resources/directional_gold.tsv                                   (Stage 4, built once by build_directional_gold.py)
  <msigdb_dir>/*.gmt                                               (Hallmark, GO:BP)

OUTPUTS (--out-dir)
  metrics_expression.tsv   one row per (scope, contrast): model pcc / rel_rmse / deg_conc with bootstrap CIs, the four baselines
                           (no_change, permuted, level_only, loco_mean) and model-minus-best-baseline with CI
  metrics_pathway.tsv      one row per (scope, contrast, gene-set collection): NES Spearman and strong-set sign agreement
  metrics_ranking.tsv      one row per (scope, tier, direction, subset): sign accuracy, AP / nDCG@K / Recall@K with
                           matched-random nulls, bootstrap CIs, leave-source-out rows
  metrics_recall_curve.tsv cumulative Recall@K (K = 1..200) per (scope, tier, direction, subset) with the matched-random null band
                           and a 95% bootstrap CI (ADDED 2026-10-17, supervisor Fig 4A)
  metrics_epsilon.tsv      GSE86204 interaction residual eps = delta_double - (delta_notch + delta_ctnnb1), real vs predicted
                           (ADDED 2026-10-17; EXPLORATORY, see the note written into the file's `note` column)
  metrics_notes.txt        what was skipped and why (a ridge baseline is NOT implemented: no training set exists for a
                           single-model pipeline; I/R/D now live in compute_ird.py — see metrics_lib.py)

CHANGED 2026-10-17: metrics_expression.tsv gained `pcc_centroid` (Systema-style centroid-subtracted PCC-delta, with CI) and the
columns `source_set` (xlsx_unit | supporting_held_out | negative_control), `perturb_class` (OE | KD), `dimensions` and `gene`
for panel (xlsx) contrasts, so results can be organised by the supervisor's 7 test dimensions with OE and KD kept apart.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metrics_lib as M

p = argparse.ArgumentParser()
p.add_argument("--module9-dir", required=True)
p.add_argument("--module9b-dir", required=True)
p.add_argument("--module9c-dir", default=None)
p.add_argument("--validation-report-combined", required=True)
p.add_argument("--screen-caps-tsv", default=None)
p.add_argument("--gold-tsv", default=None)
p.add_argument("--gmt", nargs="*", default=[], help="name=path pairs, e.g. hallmark=/path/mh.all...gmt")
p.add_argument("--out-dir", required=True)
p.add_argument("--n-boot", type=int, default=300)
p.add_argument("--n-perm", type=int, default=200)
p.add_argument("--n-null", type=int, default=1000, help="matched-random permutations for the ranking metrics")
p.add_argument("--gsea-perm", type=int, default=500)
p.add_argument("--gsea-perm-large", type=int, default=100, help="permutations for collections with > 500 sets (GO:BP)")
p.add_argument("--ks", nargs="*", type=int, default=[10, 20, 50])
p.add_argument("--min-tier-tfs", type=int, default=10)
p.add_argument("--seed", type=int, default=0)
p.add_argument("--negative-control-sample-id", default="GSE281207_Tmie_KO",
               help="sample_id of the negative-control contrast (labelled source_set=negative_control)")
p.add_argument("--recall-kmax", type=int, default=200)
args = p.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
notes = ["I/R/D (identity / rejuvenation / damage): computed by compute_ird.py (module12_metrics/ird_scores.tsv), not here.",
         "Ridge baseline: NOT computed — a ridge fit needs a perturbation training set, which a single-model (LINGER-only) pipeline does not have. "
         "Baselines present: no_change, permuted, level_only, loco_mean (+ pcc_centroid, the centroid-subtracted PCC-delta).",
         "pcc_centroid: Systema-style definition approximated with the leave-one-contrast-out mean as centroid; confirm the exact definition with the supervisor."]

# ------------------------------------------------------------------ load shift files
files = []   # (scope, sample_id, group, path)
for f in sorted(glob.glob(os.path.join(args.module9_dir, "*_shifts.tsv.gz"))):
    files.append(("population", os.path.basename(f)[:-len("_shifts.tsv.gz")], "validation", f))
for f in sorted(glob.glob(os.path.join(args.module9b_dir, "*", "*_shifts.tsv.gz"))):
    files.append((os.path.basename(os.path.dirname(f)), os.path.basename(f)[:-len("_shifts.tsv.gz")], "validation", f))
if args.module9c_dir:
    for f in sorted(glob.glob(os.path.join(args.module9c_dir, "*_shifts.tsv.gz"))):
        files.append(("population", os.path.basename(f)[:-len("_shifts.tsv.gz")], "panel", f))
if not files:
    sys.exit("no *_shifts.tsv.gz files found — rerun validate_one / validate_celltype_one / panel_validate_one with the revised rules")
print(f"{len(files)} shift files")
data = {}
real_by_sample = {}
for scope, sid, grp, f in files:
    d = pd.read_csv(f, sep="\t").dropna()
    data[(scope, sid)] = (d, grp)
    real_by_sample.setdefault(sid, pd.Series(d["real_shift"].values, index=d["gene"].values))
print(f"{len(real_by_sample)} distinct contrasts: {sorted(real_by_sample)}")


def dataset_of(sid):
    m = re.match(r"(GSE\d+)", sid)
    return m.group(1) if m else sid


def loco_vector(sid, genes):
    """mean of the z-scored real shifts of every OTHER dataset's contrasts, on `genes` (NaN-aware)."""
    cols = []
    for o, ser in real_by_sample.items():
        if dataset_of(o) == dataset_of(sid):
            continue
        s = ser.reindex(genes)
        if s.notna().sum() > 100:
            v = s.values.astype(float)
            m, sd = np.nanmean(v), np.nanstd(v)
            cols.append((v - m) / sd if sd > 0 else v * 0)
    if not cols:
        return None
    arr = np.vstack(cols)
    return np.nanmean(arr, axis=0)


rep = pd.read_csv(args.validation_report_combined, sep="\t")
verdict_of = {(r["scope"], r["sample_id"]): (r.get("pass_fail"), r.get("verdict_basis")) for _, r in rep.iterrows()}
panel_rep = None
if args.module9c_dir and os.path.exists(os.path.join(args.module9c_dir, "panel_report.tsv")):
    panel_rep = pd.read_csv(os.path.join(args.module9c_dir, "panel_report.tsv"), sep="\t", dtype={"dimensions": str}).set_index("panel_id")

# ------------------------------------------------------------------ expression metrics
def _meta_for(sid, grp):
    """source_set / perturb_class / dimensions / gene for one contrast (ADDED 2026-10-17)."""
    if grp == "panel":
        src = "xlsx_unit"
    elif sid == args.negative_control_sample_id:
        src = "negative_control"
    else:
        src = "supporting_held_out"
    out = {"source_set": src, "perturb_class": "", "dimensions": "", "gene": ""}
    if panel_rep is not None and sid in panel_rep.index:
        r_ = panel_rep.loc[sid]
        for k_, c_ in (("perturb_class", "perturb_class"), ("dimensions", "dimensions"), ("gene", "gene")):
            if c_ in panel_rep.columns and pd.notna(r_.get(c_)):
                out[k_] = r_[c_]
    return out


rows = []
for (scope, sid), (d, grp) in data.items():
    r, pr, lv = d["real_shift"].values, d["predicted_shift"].values, d["level"].values
    if len(r) < 200 or np.std(pr) == 0:
        rows.append({"scope": scope, "contrast": sid, "group": grp, "n_genes": len(r), "note": "too few genes or a no-op prediction"})
        continue
    loco = loco_vector(sid, d["gene"])
    m = M.contrast_metrics(r, pr, lv, loco, ks=(100, 300), n_perm=args.n_perm, seed=args.seed)
    m.update(M.bootstrap_contrast(r, pr, lv, loco, ks=(100, 300), B=args.n_boot, seed=args.seed))
    if grp == "validation":
        m["verdict"] = verdict_of.get((scope, sid), (None, None))[0]
    elif panel_rep is not None and sid in panel_rep.index:
        m["verdict"] = panel_rep.loc[sid, "verdict"]
    rows.append({"scope": scope, "contrast": sid, "group": grp, **_meta_for(sid, grp), **m})
expr = pd.DataFrame(rows)
expr.to_csv(os.path.join(args.out_dir, "metrics_expression.tsv"), sep="\t", index=False)
print(f"metrics_expression.tsv: {len(expr)} rows")

# ------------------------------------------------------------------ pathways
gmts = {}
for spec in args.gmt:
    name, path = spec.split("=", 1)
    if os.path.exists(path):
        gmts[name] = path
    else:
        notes.append(f"gene-set file for {name!r} not found at {path} — pathway metrics for it skipped (download_gold_standard.sh fetches it)")
prow = []
if gmts:
    universe = set()
    for (_, _), (d, _) in data.items():
        universe |= set(d["gene"])
    sets_by = {n: M.read_gmt(pth, universe=universe) for n, pth in gmts.items()}
    for n, s in sets_by.items():
        print(f"gene sets {n}: {len(s)} (size 15-500 after intersecting the gene universe)")
    for (scope, sid), (d, grp) in data.items():
        if len(d) < 200 or d["predicted_shift"].std() == 0:
            continue
        for n, sets in sets_by.items():
            perm = args.gsea_perm if len(sets) <= 500 else args.gsea_perm_large
            out, _ = M.pathway_agreement(d["real_shift"].values, d["predicted_shift"].values, d["gene"].values, sets,
                                         n_perm=perm, seed=args.seed)
            prow.append({"scope": scope, "contrast": sid, "group": grp, "collection": n, "gsea_perm": perm, **out})
pd.DataFrame(prow).to_csv(os.path.join(args.out_dir, "metrics_pathway.tsv"), sep="\t", index=False)
print(f"metrics_pathway.tsv: {len(prow)} rows")

# ------------------------------------------------------------------ ranking vs gold
rrows = []
if args.screen_caps_tsv and args.gold_tsv and os.path.exists(args.screen_caps_tsv) and os.path.exists(args.gold_tsv):
    caps = pd.read_csv(args.screen_caps_tsv, sep="\t")
    caps = caps[caps["mode"] == "ko"]
    gold = pd.read_csv(args.gold_tsv, sep="\t")
    tiers = {t: g.set_index("tf") for t, g in gold.groupby("tier")}
    scopes = list(caps["scope"].unique())
    # consensus over the screened scopes whose GAP verdict is PASS: median CAPS percentile rank per TF
    head = caps[(caps["status"] == "ok") & (caps["gap_verdict"] == "PASS")]
    if head["scope"].nunique() >= 2:
        pr_ = head.assign(pct=head.groupby("scope")["CAPS"].rank(pct=True))
        cons = pr_.groupby("tf").agg(CAPS=("pct", "median"), n=("pct", "size"), regulatory_centrality=("regulatory_centrality", "median")).reset_index()
        cons["CAPS"] = cons["CAPS"] - 0.5          # centred so the sign means above / below the within-scope median
        cons["status"] = "ok"
        cons["scope"] = "consensus_headline"
        notes.append(f"consensus_headline = median within-scope CAPS percentile rank (centred at 0.5) over {head['scope'].nunique()} GAP-PASS screened scopes")
        caps = pd.concat([caps, cons], ignore_index=True)
        scopes.append("consensus_headline")
    for scope in scopes:
        sc = caps[(caps["scope"] == scope) & (caps["status"] == "ok")][["tf", "CAPS", "regulatory_centrality"]].copy()
        if len(sc) < 50:
            continue
        for tier in sorted(tiers):
            t = tiers[tier]
            others = [tiers[o].index for o in ("B", "C") if o != tier and o in tiers] if tier in ("B", "C") else []
            for subset in (["all"] + (["leave_source_out"] if tier in ("B", "C") and others else [])):
                d = sc.copy()
                if subset == "leave_source_out":
                    excl = set().union(*[set(ix) for ix in others])
                    d = d[~d["tf"].isin(excl)]
                if tier == "D":
                    d["is_member"] = d["tf"].isin(t.index)
                    n_gold = int(d["is_member"].sum())
                    r = {"scope": scope, "tier": tier, "direction": "membership", "subset": subset, "n_gold": n_gold,
                         "underpowered": bool(n_gold < args.min_tier_tfs)}
                    r.update(M.ranking_metrics(d, ks=args.ks, B_null=args.n_null, B_boot=args.n_boot, seed=args.seed, use_abs=True))
                    rrows.append(r)
                    continue
                d["y"] = d["tf"].map(t["expected_ko_sign"])
                n_gold = int(d["y"].notna().sum())
                sa = M.sign_accuracy(d)
                rrows.append({"scope": scope, "tier": tier, "direction": "sign", "subset": subset, "n_gold": n_gold,
                              "underpowered": bool(n_gold < args.min_tier_tfs), **sa})
                for pos, label in ((1, "ko_raises_caps(+1)"), (-1, "ko_lowers_caps(-1)")):
                    n_pos = int((d["y"] == pos).sum())
                    r = {"scope": scope, "tier": tier, "direction": label, "subset": subset, "n_gold": n_gold,
                         "underpowered": bool(n_pos < args.min_tier_tfs)}
                    r.update(M.ranking_metrics(d, ks=args.ks, B_null=args.n_null, B_boot=args.n_boot, seed=args.seed, positive=pos))
                    rrows.append(r)
else:
    notes.append("ranking metrics skipped: --screen-caps-tsv and/or --gold-tsv missing")
pd.DataFrame(rrows).to_csv(os.path.join(args.out_dir, "metrics_ranking.tsv"), sep="\t", index=False)
print(f"metrics_ranking.tsv: {len(rrows)} rows")

# ------------------------------------------------------------------ cumulative recall curves (supervisor Fig 4A)
crow = []
if args.screen_caps_tsv and args.gold_tsv and os.path.exists(args.screen_caps_tsv) and os.path.exists(args.gold_tsv):
    for scope in scopes:
        sc = caps[(caps["scope"] == scope) & (caps["status"] == "ok")][["tf", "CAPS", "regulatory_centrality"]].copy()
        if len(sc) < 50:
            continue
        for tier in sorted(tiers):
            t = tiers[tier]
            d = sc.copy()
            if tier == "D":
                d["is_member"] = d["tf"].isin(t.index)
                cur = M.recall_curve(d, kmax=args.recall_kmax, seed=args.seed, use_abs=True)
                if len(cur):
                    crow.append(cur.assign(scope=scope, tier=tier, direction="membership", subset="all"))
                continue
            d["y"] = d["tf"].map(t["expected_ko_sign"])
            for pos, label in ((1, "ko_raises_caps(+1)"), (-1, "ko_lowers_caps(-1)")):
                cur = M.recall_curve(d, kmax=args.recall_kmax, seed=args.seed + (1 if pos > 0 else 2), positive=pos)
                if len(cur):
                    crow.append(cur.assign(scope=scope, tier=tier, direction=label, subset="all"))
curve_df = pd.concat(crow, ignore_index=True) if crow else pd.DataFrame(
    columns=["K", "recall", "null_mean", "null_p05", "null_p95", "ci_lo", "ci_hi", "n_pos", "n_ranked", "scope", "tier", "direction", "subset"])
curve_df.to_csv(os.path.join(args.out_dir, "metrics_recall_curve.tsv"), sep="\t", index=False)
print(f"metrics_recall_curve.tsv: {len(curve_df)} rows")

# ------------------------------------------------------------------ GSE86204 interaction residual (EXPLORATORY)
EPS_PREFIX = "GSE86204_Notch1_Ctnnb1_combo_P7__"
EPS_NOTE = ("EXPLORATORY. real eps = real(double) - [real(notch1) + real(ctnnb1_activation)], in each arm's own units. "
            "pred_additive uses the additive 'double_perturbation' arm, so its predicted eps is zero by construction (no interaction "
            "is modelled). pred_joint uses the joint Rbpj+Tcf7+Lef1 pass, the only non-additive prediction, but that pass encodes all "
            "three effectors as KNOCKOUTS whereas Ctnnb1 stabilisation is a Wnt GAIN, so its sign is not a clean counterpart; "
            "treat any correlation as descriptive. 2-3 replicates per arm; CIs resample genes only.")
erow = []
need = ["notch1_ko_only", "ctnnb1_activation_only", "double_perturbation", "double_joint_effector_ko"]
arms = {a: data.get(("population", EPS_PREFIX + a)) for a in need}
if all(v is not None for v in arms.values()):
    ser = {a: v[0].set_index("gene") for a, v in arms.items()}
    genes = ser[need[0]].index
    for a in need[1:]:
        genes = genes.intersection(ser[a].index)
    if len(genes) >= 200:
        rr = lambda a: ser[a].loc[genes, "real_shift"].values
        pp = lambda a: ser[a].loc[genes, "predicted_shift"].values
        eps_real = rr("double_perturbation") - (rr("notch1_ko_only") + rr("ctnnb1_activation_only"))
        eps_joint = pp("double_joint_effector_ko") - (pp("notch1_ko_only") + pp("ctnnb1_activation_only"))
        eps_add = pp("double_perturbation") - (pp("notch1_ko_only") + pp("ctnnb1_activation_only"))
        zr_ = M.zs(eps_real)
        rng_ = np.random.default_rng(args.seed)

        def _boot(vec):
            vals = []
            for _ in range(args.n_boot):
                i = rng_.integers(0, len(genes), len(genes))
                vals.append(M.pcc(M.zs(vec[i]), M.zs(eps_real[i])))
            return [float(x) for x in np.nanpercentile(vals, [2.5, 97.5])] if np.isfinite(vals).any() else [np.nan, np.nan]

        for name, vec in (("pred_joint", eps_joint), ("pred_additive", eps_add)):
            const = np.std(vec) == 0
            lo, hi = ([np.nan, np.nan] if const else _boot(vec))
            erow.append({"n_genes": int(len(genes)), "real_eps_sd": float(np.std(eps_real)), "prediction": name,
                         "pred_eps_sd": float(np.std(vec)), "pcc_eps": np.nan if const else M.pcc(M.zs(vec), zr_),
                         "pcc_eps_ci_lo": lo, "pcc_eps_ci_hi": hi,
                         "real_eps_to_sum_sd_ratio": float(np.std(eps_real) / max(np.std(rr("notch1_ko_only") + rr("ctnnb1_activation_only")), 1e-12)),
                         "note": EPS_NOTE})
    else:
        notes.append(f"metrics_epsilon: only {len(genes)} genes shared by the four GSE86204 arms (need >= 200) — skipped")
else:
    notes.append("metrics_epsilon: not all four GSE86204 arm shift files found in module9c_panel — skipped")
pd.DataFrame(erow, columns=["n_genes", "real_eps_sd", "prediction", "pred_eps_sd", "pcc_eps", "pcc_eps_ci_lo", "pcc_eps_ci_hi",
                            "real_eps_to_sum_sd_ratio", "note"]).to_csv(os.path.join(args.out_dir, "metrics_epsilon.tsv"), sep="\t", index=False)
print(f"metrics_epsilon.tsv: {len(erow)} rows")
open(os.path.join(args.out_dir, "metrics_notes.txt"), "w").write("\n".join(notes) + "\n")
print("\n".join(notes))
