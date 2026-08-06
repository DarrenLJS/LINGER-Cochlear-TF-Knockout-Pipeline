"""
apply_cell_type_labels.py — Module 3 final step: map manually-assigned
cell_type labels (cluster_annotation.tsv, edited from cluster_template.tsv)
onto the integrated AnnData's obs.
"""
import sys
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad(snakemake.input.integrated)
ann = pd.read_csv(snakemake.input.annotation, sep="\t", dtype={"leiden": str})

unlabeled = ann[ann["cell_type"].isna() | (ann["cell_type"].str.strip() == "")]
if len(unlabeled):
    sys.exit(
        f"cluster_annotation.tsv has {len(unlabeled)} cluster(s) with no cell_type "
        f"assigned yet: {unlabeled['leiden'].tolist()}. Fill in every row before rerunning."
    )

mapping = dict(zip(ann["leiden"], ann["cell_type"]))
missing_clusters = set(adata.obs["leiden"].astype(str)) - set(mapping.keys())
if missing_clusters:
    sys.exit(
        f"cluster_annotation.tsv is missing row(s) for cluster(s) present in the data: "
        f"{sorted(missing_clusters)}. Every leiden cluster in integrated.h5ad needs a row."
    )

adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(mapping).astype("category")
adata.write_h5ad(snakemake.output.labeled)

print("Cell type counts:")
print(adata.obs["cell_type"].value_counts().to_string())
