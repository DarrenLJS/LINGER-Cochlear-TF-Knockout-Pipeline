import scanpy as sc
adata = sc.read_h5ad("/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/module3_integration/integrated.h5ad")
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
top = sc.get.rank_genes_groups_df(adata, group=None).groupby("group").head(15)
top.to_csv("/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/module3_integration/top_de_genes_per_cluster.csv", index=False)

# also useful — sample composition per cluster, to actually check batch overlap
import pandas as pd
pd.crosstab(adata.obs["leiden"], adata.obs["sample_id"]).to_csv(
    "/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/module3_integration/leiden_by_sample.csv"
)
