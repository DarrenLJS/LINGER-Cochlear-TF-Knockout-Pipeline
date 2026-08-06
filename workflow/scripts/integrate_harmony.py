"""
integrate_harmony.py — Module 3: pool all samples' QC'd RNA, integrate across
batch (sample_id) with Harmony, cluster (Leiden), and project to UMAP.

Cell-type annotation itself stays manual by design (per the pipeline plan) —
this script produces the clustering + a cluster_template.tsv for you to fill
in, not final labels. See apply_cell_type_labels.py for the next step.
"""
import anndata as ad
import scanpy as sc
import pandas as pd

cfg = snakemake.params.integration_cfg
sample_ids = snakemake.params.sample_ids

# ---------------------------------------------------------------------------
# Load + concatenate. rna_synced.h5ad already carries QC'd, doublet-removed,
# consensus-barcode-matched cells (Module 1+2 output) — nothing further to
# filter here, just pool and batch-correct.
# ---------------------------------------------------------------------------
adatas = {}
for sid, path in zip(sample_ids, snakemake.input.rna_synced):
    a = sc.read_h5ad(path)
    a.obs["sample_id"] = sid
    a.obs["dataset"] = sid.split("_")[0]  # e.g. GSE182202_P8 -> GSE182202
    adatas[sid] = a

adata = ad.concat(
    adatas, join="outer", label="sample_id_concat_key", index_unique="-", fill_value=0
)
print(f"Pooled: {adata.shape[0]} cells x {adata.shape[1]} genes across {len(adatas)} samples "
      f"(genes UNIONED across samples via join='outer', fill_value=0 — see 2026-07-24 fix note "
      f"below; do not change back to join='inner')")
# ---------------------------------------------------------------------------
# FIX 2026-07-24 — was join="inner" (gene INTERSECTION across samples).
# Real-run evidence: cluster_marker_scores.csv came back with only 3 of the
# 6 configured marker_genes (Myo7a/Sox2/Hes1 present; Pou4f3/Gfi1/Slc26a5
# silently absent). Root cause: hair cells are ~1-5% of cells per sample, so
# a hair-cell marker can fail min_cells_per_gene=3 filtering in even ONE
# sample — inner join then drops it from the pooled dataset for ALL
# samples, disproportionately destroying exactly the rare-cell-type marker
# genes this pipeline most needs to identify hair cells by. outer join
# keeps the union instead; a gene absent from one sample is zero-filled
# there rather than deleted everywhere.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Standard normalize -> HVG -> PCA
# ---------------------------------------------------------------------------
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=cfg["n_hvgs"], batch_key=cfg["batch_key"])
adata.raw = adata
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=cfg["n_pcs"])
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

# ---------------------------------------------------------------------------
# Harmony batch correction (chosen over Scanorama, 2026-07-21)
# ---------------------------------------------------------------------------
sc.external.pp.harmony_integrate(adata, key=cfg["batch_key"])

# ---------------------------------------------------------------------------
# Neighbors on the Harmony-corrected embedding, Leiden clustering, UMAP
# ---------------------------------------------------------------------------
sc.pp.neighbors(adata, use_rep="X_pca_harmony")
sc.tl.leiden(adata, resolution=cfg["leiden_resolution"])
sc.tl.umap(adata)

n_clusters = adata.obs["leiden"].nunique()
print(f"Leiden clustering (resolution={cfg['leiden_resolution']}): {n_clusters} clusters")

adata.write_h5ad(snakemake.output.integrated)

# UMAP plot colored by cluster and by sample, side by side
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color="leiden", ax=axes[0], show=False, title="Leiden clusters")
sc.pl.umap(adata, color="sample_id", ax=axes[1], show=False, title="Sample (batch check)")
plt.tight_layout()
plt.savefig(snakemake.output.umap_png, dpi=150)
plt.close(fig)

# Template annotation TSV — cell_type column left blank for manual fill-in
template = pd.DataFrame({
    "leiden": sorted(adata.obs["leiden"].cat.categories, key=int),
    "n_cells": adata.obs["leiden"].value_counts().reindex(
        sorted(adata.obs["leiden"].cat.categories, key=int)
    ).values,
    "cell_type": "",
})
template.to_csv(snakemake.output.cluster_template, sep="\t", index=False)
print(f"Wrote cluster_template.tsv — fill in the cell_type column using "
      f"cluster_marker_scores.csv + {snakemake.output.umap_png}, save as "
      f"{cfg['cluster_annotation_tsv']}, then run apply_cell_type_labels.")
