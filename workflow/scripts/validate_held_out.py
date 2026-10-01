"""
validate_held_out.py — Module 9, one held-out/negative-control dataset's
validation score.

=============================================================================
REVISED 2026-10-01 — read this block first; it supersedes the older text
below wherever they disagree (the older text is kept as design history).
=============================================================================
For mode == "expression_shift" (reprogramming checks + negative control):
  * SCALE. Both sides are now ln(1 + CP10k) — the scale Module 4's TG
    pseudobulk (the nets' training Target) and therefore Module 8's
    predictions live on. The real side is depth-normalised PER SAMPLE/CELL
    (counts / total * 1e4, then ln1p) BEFORE averaging. v1 averaged raw
    counts and then took log2(1 + mean), i.e. no depth normalisation and a
    different log base/transform from the predictions.
  * PREDICTED SHIFT = (mean KO prediction) - (mean model baseline) x sign,
    where the model baseline is Module 8/8b's own
    _sanity_check_baseline_predicted.tsv (--sanity-tsv), i.e. the SAME net
    run without the knockout. v1 subtracted log2(1 + mean of the real
    pseudobulk), a different quantity on a different scale, so ~99.9% of
    v1's "predicted shift" was a knockout-independent gene-level term and
    only `sign` distinguished one check from another.
  * REAL SHIFT is always a within-study contrast:
      - `real_contrast` on the entry (GSE224627): cell groups picked by plate
        + cluster name from the dataset's own cluster file (see
        _load_cell_contrast);
      - else `--baseline-entry-json` (paired within-study control, e.g.
        GSE281207 HET);
      - else a cross-study fallback (population/cell-type TG pseudobulk) that
        is allowed ONLY for checks with spec["scored"] == False (descriptive).
  * VERDICT (replaces "rho > 0 and AUROC > 0.5"): gene-bootstrap CI (B =
    --bootstrap-n, seed --bootstrap-seed) of the EXCESS Spearman rho over a
    LEVEL-ONLY NULL (rho between the real shift and sign x the model's own
    baseline level — what a predictor that knows nothing but gene level
    achieves), judged against --min-excess:
        positive check : PASS if excess CI lower bound  >  min_excess,
                         FAIL if CI upper bound         <  min_excess,
                         else INCONCLUSIVE
        negative control (invert_pass): PASS if CI upper bound < min_excess,
                         FAIL if CI lower bound > min_excess, else INCONCLUSIVE
        NOT_EVALUABLE  if the predicted shift is (numerically) constant
                         (SD < --noop-pred-sd: e.g. the knocked-out TF is not
                         expressed in this scope) or < --min-genes genes
                         overlap; not_scored if spec["scored"] is False.
    The bootstrap resamples genes, which are not independent, so the CI is
    optimistic; --min-excess is the practical-significance guard. AUROC/AUPR
    are still reported as columns but no longer decide the verdict.
For mode == "aging_tf_activity": unchanged logic, but the per-cell-type
module (9b) now forwards the paired baseline too (see 09b smk) so the shift is
a same-network old-vs-young contrast.
=============================================================================

(Original docstring follows.)

Methodology per the HTML plan doc's Module 9 spec, taken at face value
rather than invented from scratch:
  - Output: "Correlation metrics, AUROC/AUPR" (stated explicitly).
  - Positive checks compare a Module 8 KNOCKOUT prediction directly
    against real overexpression/conversion data ("Atoh1/Gfi1/Pou4f3
    prediction vs. real overexpression data") — NOT a separately-built
    overexpression simulation. This script therefore reuses whichever
    Module 8 `*_predicted_expression.tsv` the check config names, applying
    a sign flip for checks where the held-out ground truth is the
    OPPOSITE direction of what Module 8 simulated (TF gain vs TF loss).
  - Aging checks use Module 7's role="baseline" regulon scores as the
    "unperturbed" reference point, per the plan doc and last session's
    revised-plan table, diffed against each aging dataset's own
    TF-activity profile (computed here the same way Module 7 does, via
    TF_activity.regulon(network="cell population")).
    REVISED 2026-09-13 per linger_cochlear_pipeline_progress02.html: the
    plan doc names GSE274279 the "ground truth for aging vector" and the
    other four aging datasets "aging-vector-support" — i.e. the four
    support datasets are meant to be checked AGAINST GSE274279's shift
    vector, not merely diffed against baseline in isolation. Two roles,
    set via check_spec["aging_role"]:
      - "reference" (GSE274279 / check=aging_vector): computes its own
        real_shift vs Module 7 baseline and writes it to
        --output-shift-tsv as "the aging vector". No external ground
        truth exists for this row, so its own rho/auroc/aupr are NaN —
        it's the reference, not a result to score.
      - "support" (check=aging_vector_support): computes its own
        real_shift the same way, then Spearman-correlates it against the
        reference vector (read from --aging-vector-tsv) and computes
        AUROC/AUPR treating the reference's top-K |shift| TFs as the
        "hit" ground truth and this dataset's |shift| as the score. This
        is the actual "does this independent aging dataset agree with
        the aging vector's direction" check the plan doc describes.
    FIXED 2026-09-26 (Bug Class A) — every aging_vector/aging_vector_support
    entry that spans multiple real ages/timepoints (GSE274279, GSE154833,
    GSE153882, GSE196870) used to load ALL of them into one entry and
    average across ages before this shift was ever computed, diluting the
    real young-vs-old signal these checks exist to detect. Each such
    dataset is now config-split into a target (old) entry plus a
    `paired_baseline_entry` (its own young/baseline timepoint, same
    study), and real_shift is computed within that study via
    --baseline-entry-json (see its help text) instead of always diffing
    against Module 7's cross-study role="baseline" mean.
  - Negative control: same shift-comparison machinery as the positive
    reprogramming checks, but the PASS condition is inverted — this
    dataset's real profile should NOT correlate with the reprogramming
    direction the way GSE224627/GSE233559's real profiles do.

  FIXED 2026-09-26 — paired within-study baseline for expression_shift.
  Previously every expression_shift check (including the GSE281207 Tmie
  negative control) computed real_shift as this dataset's OWN overall mean
  minus Module 4's population pseudobulk baseline — a cross-study
  reference built from a completely different set of experiments/batches/
  platforms. For a dataset like GSE281207 that ships its own matched
  control (P21-HET) alongside the condition of interest (P21-KO) in the
  SAME study, that is strictly worse than using the dataset's own control:
  it introduces a cross-study batch confound into exactly the comparison
  meant to test whether the model's predictions generalize, and — because
  the negative-control entry previously loaded ALL HET+KO files as one
  blended mean (Bug Class A) — the "real profile" being scored was never
  a clean KO profile in the first place.
  `--baseline-entry-json` (optional): a second config entry (same schema as
  `--entry-json`) pointing at this dataset's own within-study control
  condition (e.g. GSE281207's P21-HET entry). When given, expression_shift
  loads and log2-means BOTH entries via load_rna_only and computes
  real_shift as target_mean - control_mean, entirely within one study,
  instead of diffing against the cross-study Module 4 baseline. When
  omitted, behavior is unchanged (falls back to --baseline-tsv / Module 4's
  TG_pseudobulk.tsv), so every check that has no natural within-study
  control of its own (the reprogramming positive checks, the standalone
  held-out datasets) is unaffected by this addition.

TWO METRICS COMPUTED FOR EVERY DATASET, matching the HTML's stated output
exactly (not picking one over the other):
  1. Spearman correlation between the predicted shift vector (from Module
     8's knockout, sign-adjusted per check) and the real shift vector
     (held-out expression minus Module 4's population pseudobulk
     baseline), over the shared gene set.
  2. AUROC/AUPR, reusing Module 10's existing metric choice rather than
     inventing a new one: the top-K real-shifted genes (|real_shift|,
     K configurable) are treated as the binary ground-truth "hit" set,
     ranked against |predicted_shift| as the score.

FLAGGED, NOT GUESSED: exactly which cells/samples in each held-out
dataset represent the "converted"/"reprogrammed" state vs a baseline
state isn't specified anywhere upstream (these are whole downloaded
datasets, not pre-split into before/after). This script uses the
dataset's OVERALL mean expression profile against Module 4's population
pseudobulk baseline as the real shift — a reasonable default, but WORTH
CONFIRMING against each dataset's actual structure (e.g. if GSE224627
has explicit condition labels, a proper per-condition diff would be more
rigorous than an overall-mean shift). Uses bulk_rna_loader's same
best-effort format auto-detector as Module 7, with the same caveat about
unverified per-dataset layouts.

ADDED 2026-09-19 — --baseline-tsv / --network, for Module 9b (cell-type-
resolved validation, 09b_validation_celltype.smk). Both default to the
exact population behavior above (Module 9 itself calls this script with
neither flag set, so it is completely unaffected), so this file has ONE
implementation shared by both modules rather than a forked duplicate:
  - --baseline-tsv overrides the "before" reference used by
    mode=="expression_shift" (default: Module 4's population
    TG_pseudobulk.tsv). Module 9b points this at a cell type's own
    TG_pseudobulk_{celltype}.tsv (already written by Module 8b's
    prep_pseudobulk_target_celltype) instead — diffing a cell type's
    knockout prediction against the POPULATION baseline would conflate
    genuine cross-cell-type expression differences with the knockout's
    actual effect; diffing against that cell type's own pre-knockout
    profile isolates the knockout effect the same way Module 9 already
    does at the population level.
  - --network overrides the LingerGRN `network` argument to
    TF_activity.regulon() for mode=="aging_tf_activity" (default: "cell
    population"). Passing a cell-type name here instead is a real,
    supported LingerGRN 1.110 code path — confirmed by reading
    TF_activity.py directly: any value other than "cell population"/
    "general" is read as `cell_type_specific_trans_regulatory_{network}.
    txt`, which Module 6 (grn_celltype_specific.py) already writes, one
    per cell type, into the same --workdir this script already reads.

CAVEAT that applies to EVERY Module 9b check, not just the aging ones,
worth restating here since it's easy to lose sight of once the plumbing
works: every held-out/negative-control/aging dataset loaded by this
script (via bulk_rna_loader.load_rna_only) is BULK RNA-seq, not single-
cell/multiome — there is no per-cell-type label to subset by on the real
side of any comparison. Pointing --network or --baseline-tsv at a
specific cell type changes which MODEL (regulatory network / pre-
knockout reference) the same bulk sample is scored against, not which
CELLS from that sample are used — it answers "does this bulk sample look
consistent with cell type X's regulatory activity / knockout response",
not "what did cell type X's cells in this sample actually do" (no such
ground truth exists in this dataset collection). Report results from
this mode with that distinction stated, not implied.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from bulk_rna_loader import load_rna_only  # noqa: E402

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True, help="Module 4 workdir — data/TG_pseudobulk.tsv baseline lives here")
p.add_argument("--grn-dir", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--module8-dir", required=True)
p.add_argument("--module7-summary", required=True, help="module7_tf_activity/tf_activity_summary.tsv")
p.add_argument("--entry-json", required=True, help="held-out/negative-control config entry")
p.add_argument("--baseline-entry-json", default=None,
               help="Optional second config entry (same schema as --entry-json) pointing at "
                    "this dataset's own within-study control condition (e.g. GSE281207's "
                    "P21-HET entry vs its P21-KO --entry-json, or GSE274279's 3M-Cochlea entry "
                    "vs its 24M-Cochlea --entry-json). Used by BOTH modes: for "
                    "mode=='expression_shift', real_shift is this dataset's mean minus THIS "
                    "entry's mean instead of minus the cross-study Module 4 TG_pseudobulk.tsv "
                    "baseline; for mode=='aging_tf_activity', this entry's TF activity (via the "
                    "same regulon() call) replaces Module 7's cross-study role=\"baseline\" mean. "
                    "Mutually exclusive with --baseline-tsv (expression_shift only).")
p.add_argument("--check-spec-json", required=True,
                help='{"mode": "expression_shift"|"aging_tf_activity", '
                     '"ko_id": "<module8 ko_id>"|null, "sign": 1|-1, '
                     '"invert_pass": true|false, "top_k": 200}')
p.add_argument("--output-tsv", required=True)
p.add_argument("--output-shift-tsv", default=None,
               help="Only used when check_spec['aging_role']=='reference': writes this "
                    "dataset's real TF-activity shift vector so aging_vector_support "
                    "checks can be correlated against it.")
p.add_argument("--aging-vector-tsv", default=None,
               help="Only used when check_spec['aging_role']=='support': path to the "
                    "reference ('aging_vector') shift vector written via --output-shift-tsv.")
p.add_argument("--baseline-tsv", default=None,
               help="REVISED 2026-10-01: only the cross-study FALLBACK real reference, and only "
                    "for descriptive (spec['scored'] == False) expression_shift checks "
                    "(default: WORKDIR/data/TG_pseudobulk.tsv; Module 9b passes the cell type's "
                    "TG_pseudobulk_{celltype}.tsv). It is no longer used for the predicted shift.")
p.add_argument("--sanity-tsv", default=None,
               help="Module 8/8b _sanity_check_baseline_predicted.tsv — the model's own no-knockout "
                    "baseline. REQUIRED for mode=='expression_shift' (predicted shift = KO - this).")
p.add_argument("--bootstrap-n", type=int, default=1000)
p.add_argument("--bootstrap-seed", type=int, default=0)
p.add_argument("--min-excess", type=float, default=0.05,
               help="minimum excess Spearman rho over the level-only null that counts as signal")
p.add_argument("--noop-pred-sd", type=float, default=1e-6,
               help="predicted-shift SD below this => NOT_EVALUABLE")
p.add_argument("--min-genes", type=int, default=100)
p.add_argument("--network", default="cell population",
               help="ADDED for Module 9b. Overrides the LingerGRN `network` arg to "
                    "TF_activity.regulon() for mode=='aging_tf_activity' (default: 'cell "
                    "population'). Pass a cell type name to use Module 6's "
                    "cell_type_specific_trans_regulatory_{network}.txt instead — a real "
                    "LingerGRN 1.110 code path, not a population-only special case.")
args = p.parse_args()

entry = json.loads(args.entry_json)
spec = json.loads(args.check_spec_json)
sample_id = entry.get("sample_id")
check = entry.get("check")
top_k = spec.get("top_k", 200)

print(f"--- {sample_id} (check={check}, mode={spec['mode']}) ---")


def _auroc_aupr(real_shift, predicted_shift, top_k):
    common = real_shift.index.intersection(predicted_shift.index)
    real = real_shift.loc[common]
    pred = predicted_shift.loc[common]
    finite = real.notna() & pred.notna()
    n_dropped = int((~finite).sum())
    if n_dropped:
        print(f"_auroc_aupr: dropping {n_dropped}/{len(finite)} entries with NaN "
              f"(undefined TF-activity/shift, e.g. zero-variance regulon in a small "
              f"dataset) before scoring: {list(real.index[~finite])}")
    real, pred = real[finite], pred[finite]
    if len(real) < 10:
        return np.nan, np.nan, len(real)
    k = min(top_k, len(real) // 2)
    if k < 2:
        return np.nan, np.nan, len(real)
    threshold = real.abs().sort_values(ascending=False).iloc[k - 1]
    labels = (real.abs() >= threshold).astype(int).values
    scores = pred.abs().values
    if len(set(labels)) < 2:
        return np.nan, np.nan, len(real)
    return roc_auc_score(labels, scores), average_precision_score(labels, scores), len(real)


from scipy.stats import rankdata  # noqa: E402


def _ln_cp10k_mean(adata, label):
    """Per-sample depth normalisation, then mean over samples/cells.
    counts/total*1e4 -> ln1p -> mean. Returns (Series over var_names, n_used)."""
    X = adata.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    X = np.asarray(X, dtype=np.float64)
    X[X < 0] = 0
    tot = X.sum(axis=1)
    bad = tot <= 0
    if bad.any():
        print(f"[{label}] dropping {int(bad.sum())} sample(s)/cell(s) with zero total signal")
        X, tot = X[~bad], tot[~bad]
    if X.shape[0] == 0:
        raise ValueError(f"[{label}] no usable samples after depth normalisation")
    ln = np.log1p(X / tot[:, None] * 1e4)
    print(f"[{label}] {X.shape[0]} sample(s)/cell(s), {X.shape[1]} genes; total-signal range "
          f"{tot.min():.0f}-{tot.max():.0f}")
    return pd.Series(ln.mean(axis=0), index=[str(g) for g in adata.var_names]), int(X.shape[0])


def _load_cell_contrast(entry):
    """real_contrast on the entry (GSE224627): two cell groups chosen by
    plate + cluster name from the dataset's own cluster file. Fail-loud on any
    name that does not exist, and on a plate/cluster assignment that does not
    reproduce the biology it is supposed to carry (self-check genes must go UP
    in the target group)."""
    rc = entry["real_contrast"]
    d = entry["path"]
    cl = pd.read_csv(os.path.join(d, rc["cluster_file"]))
    cell_col, clus_col = rc.get("cell_column", "Cell_ID"), rc.get("cluster_column", "Cluster_names")
    X = pd.read_csv(os.path.join(d, rc["matrix_file"]), sep="\t", index_col=0)
    sym_col = rc["symbol_column"]
    if sym_col not in X.columns:
        raise ValueError(f"symbol_column {sym_col!r} not in matrix columns {list(X.columns)[:5]}...")
    sym = X.pop(sym_col).astype(str)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0)
    n_dup = int(sym.duplicated().sum())
    X = X.groupby(sym.values).sum()  # duplicated symbols are SUMMED, not dropped
    print(f"[contrast] matrix {X.shape[0]} unique symbols x {X.shape[1]} cells ({n_dup} duplicate symbol rows summed)")
    sep = rc.get("plate_sep", "_")
    cl["plate"] = cl[cell_col].astype(str).str.split(sep).str[0]

    def _pick(group, name):
        sub = cl[cl["plate"] == group["plate"]]
        if sub.empty:
            raise ValueError(f"[contrast:{name}] plate {group['plate']!r} not found; plates present: "
                             f"{sorted(cl['plate'].unique())}")
        have = set(sub[clus_col])
        missing = [c for c in group["clusters"] if c not in have]
        if missing:
            raise ValueError(f"[contrast:{name}] cluster name(s) {missing} not found on plate "
                             f"{group['plate']!r}. Present there: {sorted(have)}")
        cells = sub.loc[sub[clus_col].isin(group["clusters"]), cell_col].tolist()
        absent = [c for c in cells if c not in X.columns]
        if absent:
            raise ValueError(f"[contrast:{name}] {len(absent)} labelled cells missing from the matrix, e.g. {absent[:3]}")
        print(f"[contrast:{name}] plate={group['plate']} clusters={group['clusters']} -> {len(cells)} cells")
        return cells

    t_cells, r_cells = _pick(rc["target"], "target"), _pick(rc["reference"], "reference")
    if set(t_cells) & set(r_cells):
        raise ValueError("[contrast] target and reference share cells")
    mats = {}
    for nm, cells in (("target", t_cells), ("reference", r_cells)):
        a = ad_mod.AnnData(X=sp_mod.csr_matrix(X[cells].values.T))
        a.var_names = list(X.index)
        mats[nm] = _ln_cp10k_mean(a, f"contrast:{nm}")
    (t_mean, n_t), (r_mean, n_r) = mats["target"], mats["reference"]
    # biology self-check: the plate/cluster assignment must reproduce the
    # transgene/hair-cell signature, otherwise the labels are not what the
    # config says and every number downstream would be meaningless.
    min_lfc = rc.get("min_self_check_lfc", 0.5)
    for g in rc.get("self_check_up_genes", []):
        if g not in t_mean.index:
            raise ValueError(f"[contrast] self-check gene {g!r} absent from the matrix")
        lfc = float(t_mean[g] - r_mean[g])
        print(f"[contrast] self-check {g}: ln-CP10k target {t_mean[g]:.2f} vs reference {r_mean[g]:.2f} (lfc {lfc:+.2f})")
        if lfc < min_lfc:
            raise ValueError(f"[contrast] self-check FAILED for {g}: lfc {lfc:+.2f} < {min_lfc}. The target "
                             f"group is not the reprogrammed population the config describes (plate labels "
                             f"swapped?). Refusing to score.")
    for g in rc.get("report_genes", []):
        if g in t_mean.index:
            print(f"[contrast] report {g}: lfc {float(t_mean[g] - r_mean[g]):+.2f}")
    return t_mean, r_mean, n_t, n_r


def _spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den > 0 else np.nan


def _bootstrap_excess(real, pred, level, B, seed):
    rng = np.random.default_rng(seed)
    n = len(real)
    rho_b, exc_b = np.empty(B), np.empty(B)
    for i in range(B):
        idx = rng.integers(0, n, n)
        r = _spearman(real[idx], pred[idx])
        nl = _spearman(real[idx], level[idx])
        rho_b[i] = r
        exc_b[i] = r - max(nl if np.isfinite(nl) else 0.0, 0.0)
    ok = np.isfinite(exc_b)
    return (np.percentile(rho_b[ok], [2.5, 97.5]), np.percentile(exc_b[ok], [2.5, 97.5]))


def _verdict(scored, invert, evaluable, exc_lo, exc_hi, min_excess):
    if not scored:
        return "not_scored"
    if not evaluable:
        return "NOT_EVALUABLE"
    if invert:
        if exc_hi < min_excess:
            return "PASS"
        if exc_lo > min_excess:
            return "FAIL"
        return "INCONCLUSIVE"
    if exc_lo > min_excess:
        return "PASS"
    if exc_hi < min_excess:
        return "FAIL"
    return "INCONCLUSIVE"


extra = {}  # additional output columns for expression_shift rows

if spec["mode"] == "expression_shift":
    import anndata as ad_mod
    import scipy.sparse as sp_mod

    scored = spec.get("scored", True)
    if not args.sanity_tsv:
        raise ValueError("mode=='expression_shift' requires --sanity-tsv (the model's own no-knockout baseline)")

    # ---- real shift: ALWAYS a within-study contrast (or descriptive fallback) ----
    if entry.get("real_contrast"):
        t_mean, r_mean, n_t, n_r = _load_cell_contrast(entry)
        real_reference = "within-study cell-group contrast (real_contrast)"
    else:
        t_adata = load_rna_only(entry)
        t_mean, n_t = _ln_cp10k_mean(t_adata, f"{sample_id}:target")
        if args.baseline_entry_json:
            b_entry = json.loads(args.baseline_entry_json)
            b_adata = load_rna_only(b_entry)
            r_mean, n_r = _ln_cp10k_mean(b_adata, f"{b_entry.get('sample_id') or b_entry.get('ref_id')}:paired-control")
            real_reference = "paired within-study control"
        else:
            if scored:
                raise ValueError(
                    f"{sample_id}: a SCORED expression_shift check needs a within-study reference "
                    f"(real_contrast or paired_baseline_entry). Refusing to fall back to a cross-study "
                    f"pseudobulk — set check_spec['scored']=false to run it descriptively.")
            baseline_path = args.baseline_tsv or os.path.join(args.workdir, "data", "TG_pseudobulk.tsv")
            ref = pd.read_csv(baseline_path, sep=",", header=0, index_col=0)
            r_mean, n_r = ref.mean(axis=1), int(ref.shape[1])
            real_reference = f"CROSS-STUDY fallback ({os.path.basename(baseline_path)}) — descriptive only"
    print(f"real reference: {real_reference}; n_target={n_t}, n_reference={n_r}")
    common_genes = t_mean.index.intersection(r_mean.index)
    real_shift = (t_mean.loc[common_genes] - r_mean.loc[common_genes])
    print(f"real_shift: {len(real_shift)} genes (ln-CP10k units)")

    # ---- predicted shift: KO minus the model's OWN baseline, same scale ----
    ko_id = spec["ko_id"]
    ko_path = os.path.join(args.module8_dir, f"{ko_id}_predicted_expression.tsv")
    ko_pred = pd.read_csv(ko_path, sep="\t", index_col=0)
    base_pred = pd.read_csv(args.sanity_tsv, sep="\t", index_col=0)
    genes_p = ko_pred.index.intersection(base_pred.index)
    base_mean = base_pred.loc[genes_p].mean(axis=1)
    predicted_shift = (ko_pred.loc[genes_p].mean(axis=1) - base_mean) * spec["sign"]
    print(f"predicted_shift: {len(predicted_shift)} genes (KO {ko_id} - model baseline {os.path.basename(args.sanity_tsv)}, "
          f"sign={spec['sign']}); SD={float(predicted_shift.std()):.3g}, max|.|={float(predicted_shift.abs().max()):.3g}")

    common = real_shift.index.intersection(predicted_shift.index)
    real_c, pred_c = real_shift.loc[common], predicted_shift.loc[common]
    lvl_c = (base_mean.loc[common] * spec["sign"])
    finite = real_c.notna() & pred_c.notna() & lvl_c.notna()
    if (~finite).sum():
        print(f"dropping {int((~finite).sum())}/{len(finite)} NaN entries before scoring")
    real_c, pred_c, lvl_c = real_c[finite], pred_c[finite], lvl_c[finite]

    pred_sd = float(pred_c.std()) if len(pred_c) > 1 else 0.0
    evaluable = (len(real_c) >= args.min_genes) and (pred_sd >= args.noop_pred_sd)
    if len(real_c) >= 10 and pred_sd > 0:
        rho, pval = spearmanr(real_c, pred_c)
    else:
        rho, pval = np.nan, np.nan
    auroc, aupr, n_common = _auroc_aupr(real_shift, predicted_shift, top_k)
    n_common = len(real_c)

    level_null, exc, exc_lo, exc_hi, rho_lo, rho_hi = (np.nan,) * 6
    if evaluable:
        level_null = _spearman(real_c.values, lvl_c.values)
        exc = rho - max(level_null if np.isfinite(level_null) else 0.0, 0.0)
        (rho_lo, rho_hi), (exc_lo, exc_hi) = _bootstrap_excess(
            real_c.values, pred_c.values, lvl_c.values, args.bootstrap_n, args.bootstrap_seed)
        print(f"rho={rho:.3f} [{rho_lo:.3f}, {rho_hi:.3f}]  level-only null rho={level_null:.3f}  "
              f"excess={exc:.3f} [{exc_lo:.3f}, {exc_hi:.3f}] (min_excess={args.min_excess})")
    else:
        print(f"NOT EVALUABLE: n_genes={len(real_c)} (min {args.min_genes}), predicted-shift SD={pred_sd:.3g} "
              f"(min {args.noop_pred_sd:g}) — the knockout has no effect in this scope")
    verdict = _verdict(scored, spec.get("invert_pass", False), evaluable, exc_lo, exc_hi, args.min_excess)
    print(f"VERDICT: {verdict}")
    extra = {
        "verdict": verdict, "n_real_target": n_t, "n_real_reference": n_r,
        "real_reference": real_reference, "pred_shift_sd": pred_sd,
        "rho_ci_lo": rho_lo, "rho_ci_hi": rho_hi, "level_null_rho": level_null,
        "excess_rho": exc, "excess_ci_lo": exc_lo, "excess_ci_hi": exc_hi,
        "min_excess": args.min_excess,
    }

elif spec["mode"] == "aging_tf_activity":
    # --- real: held-out dataset's own TF activity, via the same LingerGRN call Module 7 uses ---
    import LingerGRN.TF_activity as TF_activity
    adata = load_rna_only(entry)
    os.chdir(args.workdir)
    real_tf_activity = TF_activity.regulon(args.workdir + "/", adata, args.grn_dir, args.network, args.genome)
    print(f"regulon network: {args.network!r}")
    real_mean = real_tf_activity.mean(axis=1)

    # FIXED 2026-09-26 (Bug Class A) — baseline used to always be Module 7's
    # role="baseline" mean TF activity (a cross-study reference built from
    # model_construction_refs, different tissue/platform/age from most
    # aging datasets). Every aging_vector/aging_vector_support entry that
    # actually contains its own young/baseline timepoint (GSE274279,
    # GSE154833, GSE153882, GSE196870) was previously loading ALL its ages
    # at once and averaging them into real_mean above BEFORE this diff —
    # blending young+old into one "aging shift" the same way GSE281207
    # blended HET+KO. Now that each such entry is config-split into a
    # target (old) + paired_baseline_entry (young, same study), reuse
    # --baseline-entry-json here too: when given, the baseline's own TF
    # activity is computed via the identical regulon() call instead of
    # Module 7's cross-study mean, so real_shift is a genuine within-study
    # old-vs-young contrast. Entries with no natural within-study young
    # timepoint of their own keep the original cross-study baseline.
    if args.baseline_entry_json:
        baseline_entry = json.loads(args.baseline_entry_json)
        baseline_adata = load_rna_only(baseline_entry)
        baseline_tf_activity = TF_activity.regulon(
            args.workdir + "/", baseline_adata, args.grn_dir, args.network, args.genome
        )
        baseline_mean = baseline_tf_activity.mean(axis=1)
        print(
            f"aging baseline: paired within-study control "
            f"sample_id={baseline_entry.get('sample_id') or baseline_entry.get('ref_id')!r} "
            f"({baseline_adata.shape[0]} sample(s)) — NOT the cross-study Module 7 baseline"
        )
    else:
        m7 = pd.read_csv(args.module7_summary, sep="\t", index_col=0)
        baseline_rows = m7[m7["role"] == "baseline"].drop(columns=["role"])
        baseline_mean = baseline_rows.mean(axis=0)  # mean over every model_construction_refs entry, per TF
        baseline_mean.index.name = None

    common_tfs = real_mean.index.intersection(baseline_mean.index)
    real_shift = real_mean.loc[common_tfs] - baseline_mean.loc[common_tfs]
    baseline_desc = "paired within-study control" if args.baseline_entry_json else "Module 7 cross-study baseline"
    print(f"aging real_shift: {len(real_shift)} TFs (held-out {sample_id} vs {baseline_desc})")

    aging_role = spec.get("aging_role")
    if aging_role == "reference":
        # This dataset IS "the aging vector" (GSE274279 per the plan doc) —
        # there's no external ground truth to score it against, so it's
        # reported descriptively (NaN rho/auroc/aupr) and written out for
        # the "support" datasets to correlate against.
        if args.output_shift_tsv:
            real_shift.rename("shift").to_frame().to_csv(args.output_shift_tsv, sep="\t")
            print(f"aging vector: wrote {len(real_shift)}-TF reference shift to {args.output_shift_tsv}")
        else:
            print("WARNING: aging_role=='reference' but --output-shift-tsv not given — "
                  "aging_vector_support checks will have nothing to correlate against")
        rho, pval = np.nan, np.nan
        auroc, aupr, n_common = np.nan, np.nan, len(real_shift)

    elif aging_role == "support":
        if not args.aging_vector_tsv:
            raise ValueError("check_spec['aging_role']=='support' requires --aging-vector-tsv "
                              "(the reference aging_vector dataset's shift file)")
        aging_vector = pd.read_csv(args.aging_vector_tsv, sep="\t", index_col=0)["shift"]
        common = real_shift.index.intersection(aging_vector.index)
        print(f"aging TF-activity shift computed for {len(real_shift)} TFs "
              f"({len(common)} shared with the reference aging vector)")
        av_c, rs_c = aging_vector.loc[common], real_shift.loc[common]
        finite = av_c.notna() & rs_c.notna()
        if (~finite).sum():
            print(f"spearman: dropping {int((~finite).sum())}/{len(finite)} NaN entries "
                  f"before correlating (undefined TF-activity, likely zero-variance "
                  f"regulon genes in this dataset): {list(av_c.index[~finite])}")
        av_c, rs_c = av_c[finite], rs_c[finite]
        if len(av_c) < 10:
            rho, pval = np.nan, np.nan
        else:
            rho, pval = spearmanr(av_c, rs_c)
        # Reference vector's top-K |shift| TFs = ground-truth "hit" set;
        # this dataset's |shift| = score — same convention as _auroc_aupr
        # uses elsewhere (real_shift=ground truth, predicted_shift=score).
        auroc, aupr, n_common = _auroc_aupr(aging_vector, real_shift, top_k)

    else:
        raise ValueError("check_spec['mode']=='aging_tf_activity' requires "
                          "check_spec['aging_role'] to be 'reference' or 'support', "
                          f"got: {aging_role!r}")

else:
    raise ValueError(f"Unknown check spec mode: {spec['mode']}")

result = {
    "sample_id": sample_id,
    "check": check,
    "mode": spec["mode"],
    "ko_id": spec.get("ko_id"),
    "sign": spec.get("sign"),
    "invert_pass": spec.get("invert_pass", False),
    "n_genes_or_tfs": n_common,
    "spearman_rho": rho,
    "spearman_pval": pval,
    "auroc": auroc,
    "aupr": aupr,
}
result.update(extra)
result.setdefault("verdict", "not_scored")
pd.DataFrame([result]).to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: rho={rho}, auroc={auroc}, aupr={aupr}")
