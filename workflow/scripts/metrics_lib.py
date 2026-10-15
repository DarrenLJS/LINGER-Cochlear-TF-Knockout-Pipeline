"""
metrics_lib.py — pure functions behind compute_metrics.py (ADDED 2026-10-15). numpy / pandas / scipy only.

DEFINITIONS (written down because the metric names come from the project decks, whose exact formulas were not available
when this was built; where a choice was made it is stated here and in compute_metrics.py's header):

Expression-shift contrasts (a real within-study contrast r vs the model's predicted shift p, aligned on genes):
  * Both vectors are z-scored (mean 0, SD 1) before any error metric. Real shifts are in ln-CP10k units and predicted shifts in
    the model's expression units, so raw RMSE across the two is not meaningful. After z-scoring the NO-CHANGE baseline
    (p = 0) has RMSE exactly 1, so `rel_rmse = RMSE(z(p), z(r))` IS the RMSE relative to no-change; < 1 beats it.
  * pcc_delta = Pearson correlation of p and r (the Spearman counterpart is already in the validation report).
  * deg_concordance@K = among the K genes with the largest |r|, the fraction whose predicted shift has the same sign.
  * Baselines scored on the SAME transformed vectors:
      no_change   p = 0                                             (rel_rmse 1, pcc 0 by construction)
      permuted    p shuffled over genes (mean and 95th percentile over `n_perm` shuffles)
      level_only  p = z(baseline expression level), sign chosen to agree with r  (the best a level-only predictor can do;
                  the same |rho| convention as the verdict's level null)
      loco_mean   leave-one-contrast-out matched mean: p = mean of the z-scored REAL shifts of every OTHER contrast
                  (a generic 'typical perturbation' shift)
    A ridge-regression baseline is NOT implemented: its definition is not known.
  * CIs: gene bootstrap (resample genes, recompute) for the model's value and for model-minus-best-baseline.

Pathway level: prerank GSEA (weighted KS, p = 1) on the signed z-score of r and of p separately, same gene sets; NES =
ES / mean(null ES of the same sign) with a size-matched random-gene-set null. Reported: Spearman correlation of the two NES
vectors across gene sets, and sign agreement among sets with |NES_real| >= nes_min.

Ranking against the directional gold table (expected KO sign +1 / -1 / membership 0):
  * score = screen CAPS of the single-TF knockout; only status == ok rows are ranked.
  * sign_accuracy: share of gold TFs (tier B/C, expected sign +/-1) whose CAPS has the expected sign; Wilson CI; chance 0.5 and
    the majority-class rate are reported next to it.
  * AP / nDCG@K / Recall@K: ranking all ok TFs by CAPS, the TFs expected +1 are the positives (unknown is treated as
    negative, which understates precision — every number is therefore reported next to its matched-random null).
    For tier D (membership only) the ranking score is |CAPS| and positives are the members.
  * matched-random null: the positive labels are permuted WITHIN deciles of regulatory centrality (so a hub-TF bias cannot
    manufacture enrichment); enrichment = observed / null mean, p = share of null >= observed.
  * leave-source-out: each tier's metric is also computed on the TFs that are NOT in the other tier.
  * CIs: bootstrap over TFs.
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, norm


# ------------------------------------------------------------------------------------------------------ expression metrics
def zs(x):
    x = np.asarray(x, dtype=float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def pcc(a, b):
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def rel_rmse(zp, zr):
    return float(np.sqrt(np.mean((zp - zr) ** 2)))


def deg_concordance(p, r, k):
    k = min(k, len(r))
    idx = np.argsort(-np.abs(r))[:k]
    pr = np.sign(p[idx])
    keep = pr != 0
    return float((pr[keep] == np.sign(r[idx][keep])).mean()) if keep.any() else np.nan


def contrast_metrics(r, p, level, loco, ks=(100, 300), n_perm=200, seed=0):
    """r, p, level: aligned 1-D arrays (genes). loco: z-scored leave-one-contrast-out mean vector or None.
    -> dict of model and baseline values (point estimates)."""
    rng = np.random.default_rng(seed)
    zr, zp, zl = zs(r), zs(p), zs(level)
    out = {"n_genes": int(len(r)), "pcc": pcc(zp, zr), "rel_rmse": rel_rmse(zp, zr)}
    for k in ks:
        out[f"deg_conc@{k}"] = deg_concordance(p, r, k)
    # no-change
    out["nochange_rel_rmse"] = rel_rmse(np.zeros_like(zr), zr)
    # level only (sign agrees with r)
    pl = pcc(zl, zr)
    s = 1.0 if (np.isnan(pl) or pl >= 0) else -1.0
    out["level_pcc"], out["level_rel_rmse"] = (abs(pl) if not np.isnan(pl) else np.nan), rel_rmse(s * zl, zr)
    for k in ks:
        out[f"level_deg_conc@{k}"] = deg_concordance(s * zl, r, k)
    # permuted
    pp, pr_, pd_ = [], [], {k: [] for k in ks}
    for _ in range(n_perm):
        sh = rng.permutation(zp)
        pp.append(pcc(sh, zr)); pr_.append(rel_rmse(sh, zr))
        for k in ks:
            pd_[k].append(deg_concordance(sh, r, k))
    out["perm_pcc_mean"], out["perm_pcc_p95"] = float(np.nanmean(pp)), float(np.nanpercentile(pp, 95))
    out["perm_rel_rmse_mean"], out["perm_rel_rmse_p05"] = float(np.mean(pr_)), float(np.percentile(pr_, 5))
    for k in ks:
        out[f"perm_deg_conc@{k}_p95"] = float(np.nanpercentile(pd_[k], 95))
    # leave-one-contrast-out matched mean
    if loco is not None and np.std(loco) > 0:
        zloco = zs(loco)
        out["loco_pcc"], out["loco_rel_rmse"] = pcc(zloco, zr), rel_rmse(zloco, zr)
        for k in ks:
            out[f"loco_deg_conc@{k}"] = deg_concordance(zloco, r, k)
    else:
        out["loco_pcc"] = out["loco_rel_rmse"] = np.nan
        for k in ks:
            out[f"loco_deg_conc@{k}"] = np.nan
    return out


def bootstrap_contrast(r, p, level, loco, ks=(100, 300), B=300, seed=0):
    """Gene bootstrap -> 95% CIs for the model's pcc / rel_rmse / deg_conc and for (model - best baseline)."""
    rng = np.random.default_rng(seed)
    n = len(r)
    zr0 = zs(r)
    rows = []
    for _ in range(B):
        i = rng.integers(0, n, n)
        rb, pb, lb = r[i], p[i], level[i]
        zr, zp, zl = zs(rb), zs(pb), zs(lb)
        row = {"pcc": pcc(zp, zr), "rel_rmse": rel_rmse(zp, zr)}
        pl = pcc(zl, zr)
        s = 1.0 if (np.isnan(pl) or pl >= 0) else -1.0
        base_pcc = [abs(pl) if not np.isnan(pl) else np.nan]
        base_rm = [rel_rmse(s * zl, zr), 1.0]
        if loco is not None and np.std(loco) > 0:
            zlo = zs(loco[i])
            base_pcc.append(pcc(zlo, zr))
            base_rm.append(rel_rmse(zlo, zr))
        row["best_baseline_pcc"] = float(np.nanmax(base_pcc)) if not np.all(np.isnan(base_pcc)) else np.nan
        row["best_baseline_rel_rmse"] = float(np.nanmin(base_rm))
        row["d_pcc"] = row["pcc"] - row["best_baseline_pcc"]
        row["d_rel_rmse"] = row["rel_rmse"] - row["best_baseline_rel_rmse"]
        for k in ks:
            row[f"deg_conc@{k}"] = deg_concordance(pb, rb, k)
        rows.append(row)
    df = pd.DataFrame(rows)
    out = {}
    for c in df.columns:
        lo, hi = np.nanpercentile(df[c], [2.5, 97.5])
        out[f"{c}_ci_lo"], out[f"{c}_ci_hi"] = float(lo), float(hi)
    return out


# ------------------------------------------------------------------------------------------------------ prerank GSEA
def read_gmt(path, min_size=15, max_size=500, universe=None):
    sets = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            genes = [g for g in f[2:] if g]
            if universe is not None:
                genes = [g for g in genes if g in universe]
            if min_size <= len(genes) <= max_size:
                sets[f[0]] = genes
    return sets


def prerank_es(score_sorted_abs, hit_idx, n):
    """weighted (p=1) enrichment score. score_sorted_abs: |score| in ranked order; hit_idx: positions of hits."""
    k = len(hit_idx)
    if k == 0 or k >= n:
        return 0.0
    hit_idx = np.sort(hit_idx)
    w = score_sorted_abs[hit_idx]
    nr = w.sum()
    if nr <= 0:
        return 0.0
    cum_hit = np.cumsum(w) / nr                        # P_hit after each hit
    misses_before = hit_idx - np.arange(k)             # number of non-hits before hit i
    p_miss_at_hit = misses_before / (n - k)            # P_miss just before / at hit i
    up = cum_hit - p_miss_at_hit                       # deviation right after hit i
    prev_cum = np.concatenate(([0.0], cum_hit[:-1]))
    down = prev_cum - p_miss_at_hit                    # deviation just before hit i (after the misses preceding it)
    # deviation at the end of the list is 0
    mx, mn = max(up.max(), 0.0), min(down.min(), 0.0)
    return float(mx if abs(mx) >= abs(mn) else mn)


def prerank_nes(scores, gene_sets, n_perm=1000, seed=0):
    """scores: pd.Series (gene -> signed score). Returns DataFrame(set, size, ES, NES)."""
    rng = np.random.default_rng(seed)
    s = scores.dropna()
    order = np.argsort(-s.values, kind="stable")
    genes = s.index.values[order]
    absw = np.abs(s.values[order])
    n = len(genes)
    gpos = {g: i for i, g in enumerate(genes)}
    null_cache = {}
    rows = []
    for name, gl in gene_sets.items():
        idx = np.array([gpos[g] for g in gl if g in gpos])
        k = len(idx)
        if k < 5 or k >= n - 5:
            continue
        es = prerank_es(absw, idx, n)
        if k not in null_cache:
            ns = np.empty(n_perm)
            for b in range(n_perm):
                ns[b] = prerank_es(absw, rng.choice(n, size=k, replace=False), n)
            null_cache[k] = ns
        ns = null_cache[k]
        same = ns[ns >= 0] if es >= 0 else ns[ns < 0]
        denom = np.abs(same).mean() if len(same) else np.nan
        rows.append({"set": name, "size": k, "ES": es, "NES": es / denom if denom and denom > 0 else np.nan})
    return pd.DataFrame(rows)


def pathway_agreement(r, p, genes, gene_sets, n_perm=500, nes_min=1.5, seed=0):
    """NES of the real and predicted shift on the same sets -> Spearman of NES and sign agreement among strong sets."""
    sr = pd.Series(zs(r), index=genes)
    sp_ = pd.Series(zs(p), index=genes)
    a = prerank_nes(sr, gene_sets, n_perm, seed).set_index("set")
    b = prerank_nes(sp_, gene_sets, n_perm, seed + 1).set_index("set")
    j = a[["NES"]].join(b[["NES"]], lsuffix="_real", rsuffix="_pred", how="inner").dropna()
    out = {"n_sets": int(len(j))}
    if len(j) > 3:
        out["nes_spearman"] = float(spearmanr(j["NES_real"], j["NES_pred"])[0])
        strong = j[j["NES_real"].abs() >= nes_min]
        out["n_strong_real"] = int(len(strong))
        out["nes_sign_agree_strong"] = float((np.sign(strong["NES_real"]) == np.sign(strong["NES_pred"])).mean()) if len(strong) else np.nan
    else:
        out["nes_spearman"] = out["nes_sign_agree_strong"] = np.nan
        out["n_strong_real"] = 0
    return out, j


# ------------------------------------------------------------------------------------------------------ ranking metrics
def average_precision(score, y):
    y = np.asarray(y, dtype=bool)
    if y.sum() == 0:
        return np.nan
    order = np.argsort(-np.asarray(score), kind="stable")
    yy = y[order]
    hits = np.cumsum(yy)
    prec = hits / (np.arange(len(yy)) + 1)
    return float((prec * yy).sum() / yy.sum())


def ndcg_at_k(score, y, k):
    y = np.asarray(y, dtype=float)
    if y.sum() == 0:
        return np.nan
    order = np.argsort(-np.asarray(score), kind="stable")[:k]
    disc = 1.0 / np.log2(np.arange(2, len(order) + 2))
    dcg = (y[order] * disc).sum()
    ideal = np.sort(y)[::-1][:k]
    idcg = (ideal * (1.0 / np.log2(np.arange(2, len(ideal) + 2)))).sum()
    return float(dcg / idcg) if idcg > 0 else np.nan


def recall_at_k(score, y, k):
    y = np.asarray(y, dtype=bool)
    if y.sum() == 0:
        return np.nan
    order = np.argsort(-np.asarray(score), kind="stable")[:k]
    return float(y[order].sum() / y.sum())


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (float(c - h), float(c + h))


def matched_permutation(score, y, strata, metric, B=1000, seed=0):
    """Permute labels within strata (centrality deciles). -> (observed, null mean, p = P(null >= obs))."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=bool)
    obs = metric(score, y)
    groups = [np.where(strata == s)[0] for s in np.unique(strata)]
    nulls = np.empty(B)
    for b in range(B):
        yp = y.copy()
        for g in groups:
            yp[g] = rng.permutation(y[g])
        nulls[b] = metric(score, yp)
    return obs, float(np.nanmean(nulls)), float((np.sum(nulls >= obs) + 1) / (B + 1))


def ranking_metrics(df, ks=(10, 20, 50), B_null=1000, B_boot=500, seed=0, positive=1, use_abs=False, strata_col="regulatory_centrality"):
    """df: one row per ok TF with columns CAPS, y (expected sign: +1/-1/0 or NaN for non-gold), regulatory_centrality.
    positive: the expected sign treated as the positive class (use +1 or -1; use_abs for membership tier, y==0 -> member)."""
    d = df.dropna(subset=["CAPS"]).copy()
    score = d["CAPS"].values * (-1.0 if positive == -1 else 1.0)
    if use_abs:
        score = np.abs(d["CAPS"].values)
        y = d["is_member"].values.astype(bool)
    else:
        y = (d["y"].values == positive)
    q = pd.qcut(d[strata_col].rank(method="first"), 10, labels=False, duplicates="drop").values if strata_col in d else np.zeros(len(d))
    out = {"n_ranked": int(len(d)), "n_pos": int(y.sum())}
    if y.sum() == 0 or y.sum() == len(y):
        return out
    ap, ap_null, ap_p = matched_permutation(score, y, q, average_precision, B_null, seed)
    out.update({"AP": ap, "AP_null_mean": ap_null, "AP_enrichment": ap / ap_null if ap_null and ap_null > 0 else np.nan, "AP_p": ap_p,
                "AP_prevalence": float(y.mean())})
    for k in ks:
        out[f"nDCG@{k}"] = ndcg_at_k(score, y, k)
        out[f"Recall@{k}"] = recall_at_k(score, y, k)
        _, rn, rp = matched_permutation(score, y, q, lambda s_, y_, k=k: recall_at_k(s_, y_, k), max(200, B_null // 5), seed + k)
        out[f"Recall@{k}_null_mean"], out[f"Recall@{k}_p"] = rn, rp
    rng = np.random.default_rng(seed + 7)
    bs = {"AP": [], **{f"nDCG@{k}": [] for k in ks}, **{f"Recall@{k}": [] for k in ks}}
    n = len(score)
    for _ in range(B_boot):
        i = rng.integers(0, n, n)
        sb, yb = score[i], y[i]
        if yb.sum() == 0:
            continue
        bs["AP"].append(average_precision(sb, yb))
        for k in ks:
            bs[f"nDCG@{k}"].append(ndcg_at_k(sb, yb, k))
            bs[f"Recall@{k}"].append(recall_at_k(sb, yb, k))
    for m, v in bs.items():
        if len(v) > 20:
            out[f"{m}_ci_lo"], out[f"{m}_ci_hi"] = (float(x) for x in np.nanpercentile(v, [2.5, 97.5]))
    return out


def sign_accuracy(df):
    d = df.dropna(subset=["CAPS", "y"])
    d = d[d["y"] != 0]
    n = len(d)
    if n == 0:
        return {"n_signed": 0}
    k = int((np.sign(d["CAPS"].values) == d["y"].values).sum())
    lo, hi = wilson(k, n)
    maj = max((d["y"] == 1).mean(), (d["y"] == -1).mean())
    return {"n_signed": n, "n_correct": k, "sign_accuracy": k / n, "sign_acc_ci_lo": lo, "sign_acc_ci_hi": hi,
            "chance": 0.5, "majority_class_rate": float(maj)}
