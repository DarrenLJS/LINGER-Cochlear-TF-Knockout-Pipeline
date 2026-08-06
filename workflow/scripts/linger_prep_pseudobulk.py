"""
linger_prep_pseudobulk.py — Module 4 step 1.

Adapts our labeled.h5ad (Module 3 output: integrated RNA with a manual
`cell_type` column) plus the per-sample atac_consensus.h5ad files into the
exact shape LingerGRN.pseudo_bulk.pseudo_bulk() expects, then writes the
files LingerGRN.LINGER_tr.get_TSS / RE_TG_dis / load_data_scNN read back in
via HARDCODED relative paths (confirmed against LingerGRN==1.110 source,
2026-07-24 — pseudo_bulk() and RE_TG_dis() read/write `./data/...`, not
`outdir`, regardless of what outdir is passed elsewhere in the API).

Because of that hardcoding, this script's --workdir MUST be the current
working directory when get_TSS/RE_TG_dis run in the next rule
(linger_get_tss) — see 04_linger_init.smk, which chdir's into it.

Real API confirmed by inspecting the installed package, not assumed:
  pseudo_bulk.pseudo_bulk(adata_RNA, adata_ATAC, singlepseudobulk)
    -> requires adata_RNA.obs['barcode'] and adata_RNA.obs['label'] to
       already be populated (cell type labels) BEFORE pseudobulking —
       this is why this rule depends on labeled.h5ad, not integrated.h5ad.

CLI, not `script:` — 2026-07-25 fix. Snakemake's `conda:` directive
misclassified the pre-built LINGER env's directory path as a NAME (not a
DIR) under real SGE execution, and `conda env export --name '<path>'`
errored since names can't contain '/'. Rather than debug that further,
04_linger_init.smk now calls this env's own python binary directly by
absolute path via `shell:`, sidestepping Snakemake's conda resolution
entirely — which also means no injected `snakemake` object, hence argparse.
"""
import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

import LingerGRN.pseudo_bulk as pseudo_bulk

p = argparse.ArgumentParser()
p.add_argument("--labeled", required=True)
p.add_argument("--atac-consensus", required=True, nargs="+")
p.add_argument("--sample-ids", required=True, nargs="+")
p.add_argument("--workdir", required=True)
p.add_argument("--singlepseudobulk", type=int, default=0)
p.add_argument("--output-done", required=True)
args = p.parse_args()

WORKDIR = args.workdir
os.makedirs(f"{WORKDIR}/data", exist_ok=True)
os.chdir(WORKDIR)

labeled_path = args.labeled
atac_paths = args.atac_consensus
sample_ids = args.sample_ids
singlepseudobulk = args.singlepseudobulk

print(f"Loading labeled RNA: {labeled_path}")
adata_RNA = sc.read_h5ad(labeled_path)
assert "cell_type" in adata_RNA.obs.columns, (
    "labeled.h5ad has no 'cell_type' column — Module 3's manual annotation "
    "step (apply_cell_type_labels) must run before this rule."
)
adata_RNA.obs["barcode"] = adata_RNA.obs_names
adata_RNA.obs["label"] = adata_RNA.obs["cell_type"].values

print("Loading + concatenating per-sample atac_consensus.h5ad...")
atac_list = []
for sid, p in zip(sample_ids, atac_paths):
    a = sc.read_h5ad(p)
    a.obs_names = [f"{bc}-{sid}" for bc in a.obs_names]
    atac_list.append(a)
adata_ATAC = ad.concat(atac_list, join="outer", index_unique=None)
adata_ATAC.obs["barcode"] = adata_RNA.obs_names.intersection(adata_ATAC.obs_names)
# Restrict both to the shared barcode set actually used in labeled.h5ad
# (labeled.h5ad was built from rna_synced.h5ad, which is already
# barcode-matched to atac_consensus.h5ad per-sample in Module 2 — this is
# just re-establishing that intersection after the two concat/reload steps).
shared = adata_RNA.obs_names.intersection(adata_ATAC.obs_names)
print(f"  {len(shared)} barcodes shared between labeled RNA and pooled ATAC")
adata_RNA = adata_RNA[shared].copy()
adata_ATAC = adata_ATAC[shared].copy()
adata_ATAC.obs["barcode"] = adata_ATAC.obs_names
adata_ATAC.var["gene_ids"] = adata_ATAC.var_names  # LINGER expects peak coords here
adata_RNA.var["gene_ids"] = adata_RNA.var_names
adata_RNA.raw = None  # pseudo_bulk() sets .raw itself

print(f"Running LINGER pseudo_bulk() (singlepseudobulk={singlepseudobulk})...")
TG_pseudobulk, RE_pseudobulk = pseudo_bulk.pseudo_bulk(adata_RNA, adata_ATAC, singlepseudobulk)

# ---------------------------------------------------------------------------
# Write the exact filenames LINGER_tr.get_TSS/RE_TG_dis/load_data_scNN read
# back via relative `./data/...` paths (NOT snakemake.output — this mirrors
# the reference tutorial's working-directory convention).
# ---------------------------------------------------------------------------
TG_pseudobulk.to_csv("data/TG_pseudobulk.tsv", sep=",")
RE_pseudobulk.to_csv("data/RE_pseudobulk.tsv", sep=",")
pd.Series(RE_pseudobulk.index).to_csv("data/Peaks.txt", header=None, index=None)
pd.Series(TG_pseudobulk.index).to_csv("data/Symbol.txt", header=None, index=None)

print(f"Pseudobulk: {TG_pseudobulk.shape[1]} pseudo-samples, "
      f"{TG_pseudobulk.shape[0]} genes, {RE_pseudobulk.shape[0]} REs")

# Output marker (the real payload is under WORKDIR/data/, per above)
with open(args.output_done, "w") as f:
    f.write(f"TG_pseudobulk: {TG_pseudobulk.shape}\nRE_pseudobulk: {RE_pseudobulk.shape}\n")

print(f"\nWrote data/*.tsv + data/Peaks.txt + data/Symbol.txt under {WORKDIR}/data/")
