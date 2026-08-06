"""
score_markers.py — Module 3: per-Leiden-cluster mean expression of cochlear
marker genes (Myo7a/Pou4f3/Gfi1 hair cell, Slc26a5 "Prestin" OHC, Sox2/Hes1
supporting cell — from config_eddie.yaml integration.marker_genes), to give
you the evidence needed to fill in cluster_template.tsv by hand. Clustering
is automated; the actual IHC vs OHC vs supporting-cell call is not, by design.
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

marker_groups = snakemake.params.marker_genes  # dict: category -> [genes]
adata = sc.read_h5ad(snakemake.input.integrated)

all_markers = [g for genes in marker_groups.values() for g in genes]
present = [g for g in all_markers if g in adata.raw.var_names]
missing = sorted(set(all_markers) - set(present))
if missing:
    print(f"WARNING: markers not found in this dataset's genes (skipped): {missing}")

raw_adata = adata.raw.to_adata()
raw_adata.obs["leiden"] = adata.obs["leiden"].values

mean_expr = sc.get.obs_df(raw_adata, keys=present + ["leiden"]).groupby("leiden").mean()
mean_expr.to_csv(snakemake.output.scores_csv)
print(f"Wrote per-cluster mean marker expression -> {snakemake.output.scores_csv}")

if present:
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(present)), max(4, 0.4 * mean_expr.shape[0])))
    im = ax.imshow(mean_expr[present].values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, rotation=90)
    ax.set_yticks(range(mean_expr.shape[0]))
    ax.set_yticklabels(mean_expr.index)
    ax.set_xlabel("marker gene")
    ax.set_ylabel("leiden cluster")
    plt.colorbar(im, label="mean log1p expression")
    plt.tight_layout()
    plt.savefig(snakemake.output.heatmap_png, dpi=150)
    plt.close(fig)
