import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/module3_integration/integrated.h5ad")

candidates = ["Prph", "Tubb3", "Snap25", "Calb2", "Cntn3", "Cntn4", "Erbb4", "Grid2"]
raw_adata = adata.raw.to_adata()
raw_adata.obs["leiden"] = adata.obs["leiden"].values

present = [g for g in candidates if g in raw_adata.var_names]
missing = sorted(set(candidates) - set(present))
if missing:
    print("Not found in var_names:", missing)

mean_expr = sc.get.obs_df(raw_adata, keys=present + ["leiden"]).groupby("leiden").mean()
mean_expr.loc[["6", "10", "16"]].to_csv(
    "/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/module3_integration/sgn_check_6_10_16.csv"
)
print(mean_expr.loc[["6", "10", "16"]])
