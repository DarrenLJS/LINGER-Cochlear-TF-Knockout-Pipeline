"""
linger_perturbation.py — Module 8 step 2 (bypass design, not perturb.py).

Reimplements the forward pass directly against Module 6's trained
{chr}_net.pt files, modeled EXACTLY on LINGER_tr.sc_nn_NN() — the real
function that TRAINED those nets (confirmed from source) — rather than on
perturb.py's LINGER_simulation()/get_simulation(), which were built for a
different (human/hg19/hg38 'LINGER' method) input-index scheme and are
confirmed structurally incompatible with scNN's real training contract in
three ways (see prep_pseudobulk_target.py's docstring and the pipeline
README's Module 8 design notes for the full comparison).

PER-GENE FORWARD PASS, replicated from LINGER_tr.sc_nn_NN():
  TFtemp = Exp.drop([gene]).values   if gene is itself a TF, else Exp.values
  REtemp = Opn.loc[RE_TGlink[gene]].values
  inputs = vstack(TFtemp, REtemp)
  inputs = per-row z-score across samples  (see NORMALIZATION below)
  y_pred = net(inputs.T)             where net = netall[gene_row_index],
                                      netall = torch.load({chr}_net.pt)

NORMALIZATION (--norm-ref)  [REVISED 2026-10-01]
  population (DEFAULT): every input row is z-scored with the mean/std of the
      POPULATION pseudobulk row — exactly the statistics the nets were
      trained against: z = (x - mu_pop) / (sd_pop + eps). A cell-type input
      therefore lands where the net actually saw it during training, instead
      of being re-centred on its own (narrower) per-subset distribution.
      Confirmed by diag_a3: live per-subset z-scoring depressed the
      cell-type sanity rho (Hair cell 0.32); population reference restores it.
      For the population run itself the two choices are identical.
  live (legacy): mean/std recomputed from THIS call's inputs, across the
      samples passed in (v1 behaviour).
  Rows whose population sd is <= eps are degenerate (constant in training,
  so z was always 0): they are set to z = 0 and counted in the log.

KNOCKOUT (--ko-mode)  [REVISED 2026-10-01]
  training_zero (DEFAULT): the knocked-out TF's raw expression is set to 0 and
      its z-score is (0 - mu_pop) / (sd_pop + eps), i.e. "this TF has zero
      expression" as the net understands it. This is level-aware: a TF that
      is highly expressed gets a large negative z, a TF that is not expressed
      in the cell type (raw already ~0) is a true no-op.
  legacy (v1): raw TF set to 0 and the row is re-normalised live; an
      all-zero row has mean 0 and sd 0 so z = 0 (mean-clamp), which is
      level-blind and amplifies near-zero rows to unit variance. Kept only
      so v1 can be reproduced (--ko-mode legacy --norm-ref live).

OVEREXPRESSION (--mode oe)  [ADDED 2026-10-01]
  The perturbed TF's raw expression is raised, per sample, to at least a
  cap: raw_new = max(raw_scope, cap), cap = the --oe-quantile (default 0.99)
  quantile of that TF's row in the POPULATION Exp. It is level-aware: a
  sample that already expresses the TF above the cap is NOT lowered, and a
  TF that is unexpressed everywhere in the population (cap ~ 0) is a no-op.
  The cap keeps the input inside the range the nets were trained on, so OE
  here is "model-capped" and will usually be weaker than a real transgene.
  OE needs --norm-ref population and --ko-mode training_zero (the z-score of
  the raised value must be taken against the training statistics).

BATCH MODE (--variants-json)  [ADDED 2026-10-01]
  Runs many perturbations in ONE pass: each chromosome net is loaded once and,
  per gene, the baseline and every variant go through the net as one stacked
  batch. Variants are a JSON list of {"id", "tfs", "mode": "ko"|"oe"}.
  Outputs --shift-tsv (genes x variants, mean over samples of
  [perturbed - baseline] prediction; a "__baseline_mean__" column holds the
  baseline mean) and --status-tsv (one row per variant). Needs population
  normalisation + training_zero (the same maths as the single run; checked
  numerically against it). No per-sample prediction tables are written.

--sanity-check MODE (run this FIRST, before trusting any real knockout):
runs the SAME forward pass with an UNMODIFIED Exp (no knockout) and
reports Spearman correlation between predicted and real (Target) values
for every gene that got a trained net, plus a numeric gate
(--gate-min-median-rho / --gate-min-frac-gt03). The gate is reported, not
enforced here: Module 9b/11 read the same predicted file and give a scope
that fails the gate confidence_weight 0.

CELL-TYPE-RESOLVED MODE (--target-path/--opn-path):
Exp for the substituted Target is built by REINDEXING to the exact TF row
order in the persisted population Exp.tsv, NOT by recomputing a fresh TFlist
via set() intersection (see prep_pseudobulk_target_celltype.py). A TF absent
from a cell type's pseudobulk is "not detected", i.e. 0. REs absent from a
cell type's Opn are likewise 0.

Fail-loud input guards: a requested knockout TF that is not in the
population TF list raises; NaN/empty predictions raise.
"""
import argparse
import ast
import json
import os
import time

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

# Same PyTorch 2.6+ weights_only fix as grn_population_training.py /
# grn_celltype_specific.py — these are our own freshly-trained checkpoints,
# not untrusted downloads, matching the trust judgment LingerGRN's own
# later releases already made for these exact call sites.
_orig_torch_load = torch.load
def _torch_load_default_unsafe(*a, **kw):
    kw.setdefault("weights_only", False)
    return _orig_torch_load(*a, **kw)
torch.load = _torch_load_default_unsafe

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True, help="Module 4/6 workdir, where {chr}_net.pt files live")
p.add_argument("--module8-dir", required=True, help="dir with Exp.tsv/RE_TGlink_resolved.tsv from prep_pseudobulk_target.py")
p.add_argument("--ko-id", required=True)
p.add_argument("--tf-list", nargs="*", default=[],
                help="TFs to perturb (zeroed for --mode ko, raised to the cap for --mode oe). "
                     "Empty list == sanity-check baseline (no perturbation).")
p.add_argument("--mode", choices=["ko", "oe"], default="ko",
                help="perturbation applied to --tf-list in a single run (see OVEREXPRESSION in the docstring).")
p.add_argument("--oe-quantile", type=float, default=0.99,
                help="OE cap = this quantile of the TF's population Exp row (default 0.99).")
p.add_argument("--variants-json", default=None,
                help="BATCH MODE: JSON list of {id, tfs, mode}; replaces --tf-list/--mode. See docstring.")
p.add_argument("--shift-tsv", default=None, help="batch mode: genes x variants mean-shift output.")
p.add_argument("--status-tsv", default=None, help="batch mode: per-variant status table.")
p.add_argument("--sanity-check", action="store_true",
                help="Also compute predicted-vs-real correlation against Target.")
p.add_argument("--target-path", default=None,
                help="Override Target pseudobulk (default: WORKDIR/data/TG_pseudobulk.tsv, population-level). "
                     "Pass a cell-type pseudobulk (see prep_pseudobulk_target_celltype.py) for cell-type-resolved mode.")
p.add_argument("--opn-path", default=None,
                help="Override Opn (RE accessibility) pseudobulk (default: WORKDIR/data/RE_pseudobulk.tsv, population-level).")
p.add_argument("--pop-opn-path", default=None,
                help="Population Opn used ONLY for the population normalisation statistics "
                     "(default: WORKDIR/data/RE_pseudobulk.tsv).")
p.add_argument("--norm-ref", choices=["population", "live"], default="population",
                help="z-score reference: 'population' = training statistics (default); 'live' = v1 per-call statistics.")
p.add_argument("--ko-mode", choices=["training_zero", "legacy"], default="training_zero",
                help="knockout semantics: 'training_zero' (default) or v1 'legacy' mean-clamp (see docstring).")
p.add_argument("--gate-min-median-rho", type=float, default=0.5)
p.add_argument("--gate-min-frac-gt03", type=float, default=0.7)
p.add_argument("--output-tsv", default=None, help="single-run predicted expression (required unless --variants-json)")
args = p.parse_args()

BATCH = args.variants_json is not None
if BATCH:
    if not (args.shift_tsv and args.status_tsv):
        raise ValueError("--variants-json requires --shift-tsv and --status-tsv")
    if args.tf_list:
        raise ValueError("--variants-json replaces --tf-list; do not pass both")
    if args.sanity_check:
        raise ValueError("--sanity-check is a single-run option, not available in batch mode")
    if args.norm_ref != "population" or args.ko_mode != "training_zero":
        raise ValueError("batch mode needs --norm-ref population and --ko-mode training_zero")
else:
    if not args.output_tsv:
        raise ValueError("--output-tsv is required unless --variants-json is given")
if args.mode == "oe" and (args.norm_ref != "population" or args.ko_mode != "training_zero"):
    raise ValueError("--mode oe needs --norm-ref population and --ko-mode training_zero")
if not (0.5 <= args.oe_quantile <= 1.0):
    raise ValueError(f"--oe-quantile must be in [0.5, 1], got {args.oe_quantile}")

WORKDIR = args.workdir
eps = 1e-6
_T0 = time.perf_counter()


def _tick(msg):
    print(f"[timing] {msg}: {time.perf_counter() - _T0:.1f}s since start", flush=True)

# Persisted population Exp — always loaded, since its row (TF) order is what
# every {chr}_net.pt was trained against and its row statistics are the
# population normalisation reference.
Exp_population = pd.read_csv(os.path.join(args.module8_dir, "Exp.tsv"), sep="\t", index_col=0)

target_path = args.target_path or os.path.join(WORKDIR, "data", "TG_pseudobulk.tsv")
opn_path = args.opn_path or os.path.join(WORKDIR, "data", "RE_pseudobulk.tsv")
Opn = pd.read_csv(opn_path, sep=",", header=0, index_col=0)
Target = pd.read_csv(target_path, sep=",", header=0, index_col=0)
CELLTYPE_MODE = args.target_path is not None

if not CELLTYPE_MODE:
    Exp = Exp_population
else:
    # Reindex to the population Exp's exact TF row order; missing TFs are 0
    # (not detected). Leaving NaN would silently poison every prediction
    # (a NaN in any input feature propagates through the dense MLP).
    missing_tfs = [tf for tf in Exp_population.index if tf not in Target.index]
    if missing_tfs:
        print(f"WARNING: {len(missing_tfs)} of {len(Exp_population.index)} population TFs "
              f"not present in the substituted Target ({target_path}) — filling as 0 "
              f"(not detected/expressed in this cell type's pseudobulk). Missing: {missing_tfs}")
    Exp = Target.reindex(Exp_population.index).fillna(0)
    print(f"Cell-type-resolved mode: Exp reindexed to population's {len(Exp_population.index)}-TF "
          f"order from {target_path} ({Exp.shape[1]} pseudo-samples).")

RE_TGlink = pd.read_csv(os.path.join(args.module8_dir, "RE_TGlink_resolved.tsv"), sep="\t")
re_col = [c for c in RE_TGlink.columns if c not in ("gene", "chr")][0]
RE_TGlink[re_col] = RE_TGlink[re_col].apply(ast.literal_eval)

print(f"--- {args.ko_id} ({'BATCH of ' + args.variants_json if BATCH else args.mode + ' TFs: ' + str(args.tf_list or '[none — baseline]')}) "
      f"[ko-mode={args.ko_mode}, norm-ref={args.norm_ref}] ---")
_tick("inputs loaded")

# ---- fail-loud knockout guard -------------------------------------------------
missing = [tf for tf in args.tf_list if tf not in Exp.index]
if missing:
    raise ValueError(f"knockout TF(s) {missing} are not in the population TF list "
                     f"({Exp.shape[0]} TFs) — refusing to run a silent no-op knockout.")
KO_SET = set(args.tf_list)

# ---- population normalisation statistics -------------------------------------
# TF rows: from the persisted population Exp (pre-knockout). RE rows: from the
# population Opn, restricted to the REs actually used (chunked read keeps the
# cell-type jobs from holding the whole population matrix in memory).
needed_res = sorted({r for lst in RE_TGlink[re_col] for r in lst})

def _row_stats(arr):
    a = np.asarray(arr, dtype=np.float64)
    return a.mean(axis=1), a.std(axis=1, ddof=1)

tf_mu_all, tf_sd_all = _row_stats(Exp_population.values)

re_mu = pd.Series(dtype=np.float64)
re_sd = pd.Series(dtype=np.float64)
if args.norm_ref == "population" or (args.ko_mode == "training_zero" and KO_SET):
    if not CELLTYPE_MODE or (args.opn_path is None):
        sub = Opn.reindex(needed_res).dropna(how="all")
        m, s = _row_stats(sub.values)
        re_mu, re_sd = pd.Series(m, index=sub.index), pd.Series(s, index=sub.index)
    else:
        pop_opn = args.pop_opn_path or os.path.join(WORKDIR, "data", "RE_pseudobulk.tsv")
        need = set(needed_res)
        mus, sds = [], []
        for chunk in pd.read_csv(pop_opn, sep=",", header=0, index_col=0, chunksize=40000):
            sub = chunk.loc[chunk.index.isin(need)]
            if len(sub):
                m, s = _row_stats(sub.values)
                mus.append(pd.Series(m, index=sub.index)); sds.append(pd.Series(s, index=sub.index))
        re_mu, re_sd = pd.concat(mus), pd.concat(sds)
    re_mu = re_mu[~re_mu.index.duplicated()]
    re_sd = re_sd[~re_sd.index.duplicated()]
    print(f"population normalisation statistics: {len(tf_mu_all)} TF rows, {len(re_mu)} RE rows")

_tick("population statistics ready")

# ---- OE cap: quantile of the population Exp row (see OVEREXPRESSION) ----------
oe_cap = pd.Series(np.quantile(Exp_population.values.astype(np.float64), args.oe_quantile, axis=1),
                   index=Exp_population.index)

# =============================================================================
# BATCH MODE — every variant + the baseline through each net in one stacked batch
# =============================================================================
if BATCH:
    with open(args.variants_json) as fh:
        variants = json.load(fh)
    if not isinstance(variants, list) or not variants:
        raise ValueError(f"{args.variants_json}: expected a non-empty JSON list of {{id, tfs, mode}}")
    ids = [v["id"] for v in variants]
    if len(set(ids)) != len(ids):
        raise ValueError(f"duplicate variant ids in {args.variants_json}")
    if "__baseline_mean__" in ids:
        raise ValueError("variant id '__baseline_mean__' is reserved")
    tf_index = {t: i for i, t in enumerate(Exp.index)}
    for v in variants:
        if v.get("mode") not in ("ko", "oe"):
            raise ValueError(f"variant {v.get('id')!r}: mode must be 'ko' or 'oe', got {v.get('mode')!r}")
        if not v.get("tfs"):
            raise ValueError(f"variant {v['id']!r}: empty tfs list")
        miss = [t for t in v["tfs"] if t not in tf_index]
        if miss:
            raise ValueError(f"variant {v['id']!r}: TF(s) {miss} not in the population TF list "
                             f"({Exp.shape[0]} TFs) — refusing to run a silent no-op.")

    S = Exp.shape[1]
    raw = Exp.values.astype(np.float64)                       # (T, S)
    deg_tf = tf_sd_all <= eps
    z_base = (raw - tf_mu_all[:, None]) / (tf_sd_all[:, None] + eps)
    z_base[deg_tf] = 0.0

    def _z_row(i, raw_row):
        return np.zeros(S) if deg_tf[i] else (raw_row - tf_mu_all[i]) / (tf_sd_all[i] + eps)

    Zs = [z_base]                                             # index 0 = baseline
    n_changed = []
    for v in variants:
        Zv = z_base.copy()
        changed = 0
        for t in v["tfs"]:
            i = tf_index[t]
            new_raw = np.zeros(S) if v["mode"] == "ko" else np.maximum(raw[i], oe_cap.iloc[i])
            changed += int((new_raw != raw[i]).sum())
            Zv[i] = _z_row(i, new_raw)
        Zs.append(Zv)
        n_changed.append(changed)
    Zs_t = [torch.tensor(z, dtype=torch.float32) for z in Zs]
    NV = len(Zs)                                              # baseline + variants

    shift_cols = {v: {} for v in ids}
    base_mean = {}
    missing_nets_b = 0
    for chrom in sorted(RE_TGlink["chr"].unique()):
        net_path = os.path.join(WORKDIR, f"{chrom}_net.pt")
        if not os.path.exists(net_path):
            print(f"WARNING: {net_path} not found — skipping chromosome {chrom} entirely "
                  f"({(RE_TGlink['chr'] == chrom).sum()} genes affected)")
            continue
        netall = torch.load(net_path)
        chr_rows = RE_TGlink[RE_TGlink["chr"] == chrom].reset_index(drop=True)
        for idx in range(chr_rows.shape[0]):
            if idx not in netall:
                missing_nets_b += 1
                continue
            gene = chr_rows.loc[idx, "gene"]
            re_list = chr_rows.loc[idx, re_col]
            keep = (Exp.index != gene) if gene in Exp.index else np.ones(Exp.shape[0], dtype=bool)
            re_series = Opn.reindex(re_list)
            REtemp = re_series.fillna(0).values
            mu_re = re_mu.reindex(re_list).fillna(0.0).values
            sd_re = re_sd.reindex(re_list).fillna(0.0).values
            z_re = (REtemp - mu_re[:, None]) / (sd_re[:, None] + eps)
            deg = sd_re <= eps
            z_re[deg] = 0.0
            z_re_t = torch.tensor(z_re, dtype=torch.float32)
            keep_t = torch.tensor(keep)
            batch = torch.cat([torch.cat((Zt[keep_t], z_re_t), dim=0).T for Zt in Zs_t], dim=0)  # (NV*S, F)
            net = netall[idx]
            net.eval()
            with torch.no_grad():
                y = net(batch).detach().numpy().reshape(NV, S)
            m = y.mean(axis=1)
            base_mean[gene] = m[0]
            for k, vid in enumerate(ids):
                shift_cols[vid][gene] = m[k + 1] - m[0]
        _tick(f"{chrom} done ({chr_rows.shape[0]} genes, {NV} forward variants incl. baseline)")

    if missing_nets_b:
        print(f"NOTE: {missing_nets_b} genes had no trained net — excluded, not an error here.")
    if not base_mean:
        raise RuntimeError("no gene produced a prediction — check net files / RE_TGlink_resolved.tsv")
    shift_df = pd.DataFrame(shift_cols)
    shift_df["__baseline_mean__"] = pd.Series(base_mean)
    if shift_df.isna().any().any():
        raise RuntimeError(f"{int(shift_df.isna().sum().sum())} NaN values in batch output — "
                           f"refusing to write a poisoned file.")
    rows = []
    for k, v in enumerate(variants):
        # A variant whose input rows are identical to the baseline's (TF unchanged in
        # this scope, or a degenerate z=0 row) is an EXACT no-op: set its shift to 0
        # so float32 batch-offset noise (~1e-7) cannot make it look like a tiny effect.
        exact_noop = n_changed[k] == 0 or all(deg_tf[tf_index[t]] for t in v["tfs"])
        if exact_noop:
            shift_df[v["id"]] = 0.0
        col = shift_df[v["id"]].values
        rows.append({"variant": v["id"], "mode": v["mode"], "tfs": ",".join(v["tfs"]),
                     "n_samples_changed": n_changed[k],
                     "shift_l2": float(np.sqrt((col ** 2).sum())),
                     "shift_max_abs": float(np.abs(col).max()),
                     "noop": bool(exact_noop or np.abs(col).max() < 1e-6)})
    shift_df.to_csv(args.shift_tsv, sep="\t")
    pd.DataFrame(rows).to_csv(args.status_tsv, sep="\t", index=False)
    n_noop = sum(r["noop"] for r in rows)
    print(f"Wrote {args.shift_tsv}: {shift_df.shape[0]} genes x {len(ids)} variants (+ baseline); "
          f"{n_noop} variant(s) are exact no-ops (perturbed TF unchanged in this scope). Status: {args.status_tsv}")
    _tick("batch finished")
    raise SystemExit(0)

# ---- apply the perturbation to a COPY ------------------------------------------
Exp_ko = Exp.copy()
for tf in args.tf_list:
    if args.mode == "ko":
        Exp_ko.loc[tf] = 0
    else:
        new_raw = np.maximum(Exp.loc[tf].values, oe_cap[tf])
        print(f"OE {tf}: cap (q{args.oe_quantile}) = {oe_cap[tf]:.4g}; "
              f"{int((new_raw != Exp.loc[tf].values).sum())}/{Exp.shape[1]} samples raised")
        Exp_ko.loc[tf] = new_raw

predicted = {}  # gene -> np.array of predicted values across samples
missing_nets = 0
n_degenerate_rows = 0

for chrom in sorted(RE_TGlink["chr"].unique()):
    net_path = os.path.join(WORKDIR, f"{chrom}_net.pt")
    if not os.path.exists(net_path):
        print(f"WARNING: {net_path} not found — skipping chromosome {chrom} entirely "
              f"({(RE_TGlink['chr'] == chrom).sum()} genes affected)")
        continue
    netall = torch.load(net_path)
    chr_rows = RE_TGlink[RE_TGlink["chr"] == chrom].reset_index(drop=True)

    for idx in range(chr_rows.shape[0]):
        if idx not in netall:
            missing_nets += 1
            continue  # gene had no net trained (sc_nn_NN's own good==0 path, if it ever applies here)
        gene = chr_rows.loc[idx, "gene"]
        re_list = chr_rows.loc[idx, re_col]

        if gene in Exp_ko.index:
            keep = Exp_ko.index != gene
        else:
            keep = np.ones(Exp_ko.shape[0], dtype=bool)
        TFtemp = Exp_ko.values[keep]
        tf_names = Exp_ko.index[keep]
        # A peak absent from a cell-type-restricted Opn means "no accessibility
        # signal observed here" (closed chromatin) — 0, not unknown.
        re_series = Opn.reindex(re_list)
        missing_res = re_series.index[re_series.isna().any(axis=1)].tolist()
        if missing_res:
            print(f"WARNING: gene {gene}: {len(missing_res)} of {len(re_list)} REs not "
                  f"present in Opn ({opn_path}) — filling as 0 (no accessibility signal "
                  f"in this pseudobulk).")
        REtemp = re_series.fillna(0).values
        inputs = torch.tensor(np.vstack((TFtemp, REtemp)), dtype=torch.float32)

        # population statistics for this gene's input rows (same row order)
        if args.norm_ref == "population" or args.ko_mode == "training_zero":
            mu_pop = np.concatenate([tf_mu_all[keep], re_mu.reindex(re_list).fillna(0.0).values])
            sd_pop = np.concatenate([tf_sd_all[keep], re_sd.reindex(re_list).fillna(0.0).values])
            mu_t = torch.tensor(mu_pop, dtype=torch.float32)
            sd_t = torch.tensor(sd_pop, dtype=torch.float32)

        if args.norm_ref == "population":
            degenerate = sd_t <= eps
            z = (inputs - mu_t[:, None]) / (sd_t[:, None] + eps)
            z[degenerate] = 0.0
            n_degenerate_rows += int(degenerate.sum())
            inputs = z
        else:
            mean = inputs.mean(dim=1)
            std = inputs.std(dim=1)
            inputs = ((inputs.T - mean) / (std + eps)).T
            if args.ko_mode == "training_zero" and KO_SET:
                ko_rows = [i for i, n in enumerate(tf_names) if n in KO_SET]
                if ko_rows:
                    inputs[ko_rows] = ((0.0 - mu_t[ko_rows]) / (sd_t[ko_rows] + eps))[:, None]

        net = netall[idx]
        net.eval()
        with torch.no_grad():
            y_pred = net(inputs.T)
        predicted[gene] = y_pred.detach().numpy().reshape(-1)
    _tick(f"{chrom} done ({chr_rows.shape[0]} genes)")

if missing_nets:
    print(f"NOTE: {missing_nets} genes had no trained net for their chromosome "
          f"(sc_nn's own per-gene training-failure path) — excluded from output, not an error here.")
if args.norm_ref == "population":
    print(f"NOTE: {n_degenerate_rows} degenerate input rows (population sd <= eps) set to z=0 across all genes.")

if not predicted:
    raise RuntimeError("no gene produced a prediction — check net files / RE_TGlink_resolved.tsv")
pred_df = pd.DataFrame(predicted).T
pred_df.columns = Target.columns
if pred_df.isna().any().any():
    raise RuntimeError(f"{int(pred_df.isna().sum().sum())} NaN predictions — an input feature was NaN; "
                       f"refusing to write a poisoned output.")
pred_df.to_csv(args.output_tsv, sep="\t")
print(f"Wrote {args.output_tsv}: {pred_df.shape[0]} genes x {pred_df.shape[1]} samples")
_tick("single run finished")

if args.sanity_check:
    common = [g for g in pred_df.index if g in Target.index]
    rhos = []
    for g in common:
        rho, _ = spearmanr(pred_df.loc[g], Target.loc[g])
        if not np.isnan(rho):
            rhos.append(rho)
    if rhos:
        rhos = np.array(rhos)
        med, frac = float(np.median(rhos)), float((rhos > 0.3).mean())
        print(f"\nSANITY CHECK — predicted vs real Target, {len(rhos)} genes with a trained net, "
              f"{pred_df.shape[1]} pseudo-samples:")
        print(f"  median Spearman rho = {med:.3f}   mean = {rhos.mean():.3f}   frac rho>0.3 = {frac:.2f}")
        ok = (med >= args.gate_min_median_rho) and (frac >= args.gate_min_frac_gt03)
        print(f"  SANITY GATE (median rho >= {args.gate_min_median_rho}, frac rho>0.3 >= {args.gate_min_frac_gt03}): "
              f"{'PASS' if ok else 'FAIL'}")
        if CELLTYPE_MODE:
            print("  INTERPRETATION: cell-type mode. The TF ordering is inherited from the population run, so a "
                  "low value here points at the input distribution (this cell type's pseudobulk being far from "
                  "what the nets were trained on) rather than at TF ordering. A scope that FAILS the gate is "
                  "given confidence_weight 0 downstream; its numbers are not to be trusted.")
        else:
            print("  INTERPRETATION: these nets were fit to reproduce Target from this exact Exp/Opn input, so a "
                  "low correlation here most likely means the reconstructed Exp/TF ordering (see "
                  "prep_pseudobulk_target.py's docstring, and check PYTHONHASHSEED=0) does NOT match what "
                  "{chr}_net.pt was trained against. Do not trust real knockout output until this passes.")
    else:
        print("\nSANITY CHECK: no overlapping genes scored — check Target/predicted gene name alignment.")
        print("  SANITY GATE: FAIL")
