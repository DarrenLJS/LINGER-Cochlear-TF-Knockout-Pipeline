"""
sanity_checks.py — post-Module-2 sanity checks across all samples, before
trusting them as Module 3 input. Ported from the original sanity_checks.py.

Two changes vs. the original:
  1. Reads rna_synced.h5ad (Module 2d's new output) instead of overwritten
     rna_filtered.h5ad — see sync_rna_atac_barcodes.py docstring.
  2. EXPECTED_N_PEAKS is computed from consensus_peaks.bed's line count at
     run time instead of a hardcoded constant. The original constant
     (224719) had already gone stale once (was 235228 before GSE182202_P1
     was excluded from SAMPLES) and required a manual code edit + comment
     every time SAMPLES membership changed. Deriving it from the actual
     consensus file removes that maintenance burden and that whole class of
     bug entirely.

Checks performed, per sample:
  1. rna_synced.h5ad and atac_consensus.h5ad have IDENTICAL barcodes.
  2. atac_consensus.h5ad has exactly EXPECTED_N_PEAKS n_vars (computed from
     consensus_peaks.bed, not hardcoded).
  3. Cell retention vs Module 1's atac_filtered.h5ad (pre-consensus).
  4. Doublet score distribution (from Module 1's Scrublet .uns, carried
     through) — flags any sample above 2x expected_doublet_rate.
  5. Per-cell QC metric distributions — one grid figure per sample.
  6. Consensus peak width distribution — one pooled figure.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QC_PARAMS = snakemake.params.qc_params
QC_DIR = snakemake.params.qc_dir
DOUBLET_FLAG_MULTIPLIER = 2.0

with open(snakemake.input.consensus_bed) as f:
    EXPECTED_N_PEAKS = sum(1 for _ in f)
print(f"EXPECTED_N_PEAKS (from consensus_peaks.bed): {EXPECTED_N_PEAKS}")


def barcode_check(rna, atac):
    if rna is None or atac is None:
        return {"barcodes_match": "MISSING_FILE", "n_obs_rna": None, "n_obs_atac": None}
    same_set = set(rna.obs_names) == set(atac.obs_names)
    same_order = list(rna.obs_names) == list(atac.obs_names)
    status = "OK" if same_order else ("SAME_SET_DIFF_ORDER" if same_set else "MISMATCH")
    return {"barcodes_match": status, "n_obs_rna": rna.n_obs, "n_obs_atac": atac.n_obs}


def doublet_stats(rna):
    if rna is None or "doublet_rate_observed" not in rna.uns:
        return {"doublet_rate": None, "doublet_flag": None}
    rate = float(rna.uns["doublet_rate_observed"])
    expected = QC_PARAMS.get("scrublet_expected_doublet_rate", 0.06)
    flag = rate > DOUBLET_FLAG_MULTIPLIER * expected
    return {"doublet_rate": round(rate, 4), "doublet_flag": flag}


def qc_grid_plot(sample_id, rna, atac):
    if rna is None:
        return
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle(f"{sample_id} — per-cell QC (n={rna.n_obs})")

    metric_map = [
        ("n_genes_by_counts", "n genes / cell", axes[0, 0], rna),
        ("total_counts", "RNA counts / cell", axes[0, 1], rna),
        ("pct_counts_mt", "% mito", axes[1, 0], rna),
    ]
    for col, label, ax, adata in metric_map:
        if col in adata.obs.columns:
            ax.hist(adata.obs[col], bins=50, color="#4C72B0")
            ax.set_xlabel(label)
            ax.set_ylabel("cells")
        else:
            ax.text(0.5, 0.5, f"'{col}' not in .obs", ha="center", va="center")
            ax.axis("off")

    ax = axes[1, 1]
    if atac is not None and "n_peaks" in atac.obs.columns:
        ax.hist(atac.obs["n_peaks"], bins=50, color="#55A868")
        ax.set_xlabel("n peaks / cell (consensus)")
        ax.set_ylabel("cells")
    elif atac is not None:
        peaks_per_cell = np.asarray((atac.X > 0).sum(axis=1)).ravel()
        ax.hist(peaks_per_cell, bins=50, color="#55A868")
        ax.set_xlabel("n peaks / cell (consensus, computed from X)")
        ax.set_ylabel("cells")
    else:
        ax.text(0.5, 0.5, "atac_consensus.h5ad missing", ha="center", va="center")
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(f"{QC_DIR}/{sample_id}_qc_grid.png", dpi=150)
    plt.close(fig)


rows = []
all_peak_widths = []
doublet_hist_data = {}

for sample_id, rna_path, atac_path, atac_pre_path in zip(
    snakemake.params.sample_ids,
    snakemake.input.rna_synced,
    snakemake.input.atac_consensus,
    snakemake.input.atac_filtered,
):
    print(f"--- {sample_id} ---")
    rna = sc.read_h5ad(rna_path)
    atac = sc.read_h5ad(atac_path)
    atac_pre = sc.read_h5ad(atac_pre_path)

    row = {"sample_id": sample_id}
    row.update(barcode_check(rna, atac))
    row.update(doublet_stats(rna))

    row["n_vars_atac"] = atac.n_vars
    row["n_vars_matches_expected"] = (atac.n_vars == EXPECTED_N_PEAKS)
    if "chrom" in atac.var.columns and "start" in atac.var.columns and "end" in atac.var.columns:
        widths = (atac.var["end"] - atac.var["start"]).values
    else:
        try:
            widths = []
            for v in atac.var_names:
                chrom, coords = v.split(":")
                s, e = coords.split("-")
                widths.append(int(e) - int(s))
            widths = np.array(widths)
        except Exception:
            widths = np.array([])
    if widths.size:
        row["peak_width_median"] = int(np.median(widths))
        row["peak_width_max"] = int(np.max(widths))
        all_peak_widths.append(widths)

    row["n_obs_pre_consensus"] = atac_pre.n_obs
    row["n_obs_post_consensus"] = atac.n_obs
    row["pct_retained"] = round(100 * atac.n_obs / atac_pre.n_obs, 2) if atac_pre.n_obs else None

    rows.append(row)
    qc_grid_plot(sample_id, rna, atac)

    if "doublet_score" in rna.obs.columns:
        doublet_hist_data[sample_id] = rna.obs["doublet_score"].values

summary = pd.DataFrame(rows)
summary.to_csv(snakemake.output.summary, index=False)
print(f"\nWrote {snakemake.output.summary}")

problems = summary[
    (summary["barcodes_match"] != "OK")
    | (summary["n_vars_matches_expected"] != True)
    | (summary["doublet_flag"] == True)
]
if len(problems):
    print("\n*** FLAGGED SAMPLES — inspect before Module 3 ***")
    print(problems.to_string(index=False))
else:
    print("\nNo samples flagged: barcodes match, peak counts match, doublet rates within 2x expected.")

if all_peak_widths:
    pooled = np.concatenate(all_peak_widths)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(pooled, bins=80, color="#4C72B0")
    ax.set_xlabel("consensus peak width (bp)")
    ax.set_ylabel("peaks")
    ax.set_title(f"Consensus peak width distribution (n={len(pooled)})")
    plt.tight_layout()
    plt.savefig(f"{QC_DIR}/peak_width_distribution.png", dpi=150)
    plt.close(fig)

if doublet_hist_data:
    fig, ax = plt.subplots(figsize=(8, 5))
    for sid, scores in doublet_hist_data.items():
        ax.hist(scores, bins=40, histtype="step", label=sid, linewidth=1.5)
    expected = QC_PARAMS.get("scrublet_expected_doublet_rate", 0.06)
    ax.set_xlabel("Scrublet doublet score")
    ax.set_ylabel("cells")
    ax.set_title(f"Doublet score distributions (expected rate={expected})")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{QC_DIR}/doublet_scores.png", dpi=150)
    plt.close(fig)

print(f"\nAll plots written to {QC_DIR}")
