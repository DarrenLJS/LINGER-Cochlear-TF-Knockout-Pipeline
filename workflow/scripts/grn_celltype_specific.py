"""
grn_celltype_specific.py — Module 6 step 2 (per-cell-type GRN).

Real LingerGRN==1.110 calls (LL_net.py), one cell type at a time, matching
docs/GRN_infer.md's documented usage pattern:

  LL_net.cell_type_specific_TF_RE_binding(GRNdir, adata_RNA, adata_ATAC,
                                           genome, celltype, outdir, method)
  LL_net.cell_type_specific_cis_reg(GRNdir, adata_RNA, adata_ATAC, genome,
                                     method, outdir)      # NOTE: no celltype
                                                           # arg in this one
                                                           # per source — it
                                                           # reads adata's
                                                           # .obs['label']
                                                           # internally
  LL_net.cell_type_specific_trans_reg(GRNdir, adata_RNA, celltype, outdir)

Outputs per cell type (per source, filenames LINGER writes itself):
  cell_type_specific_cis_regulatory_{celltype}.txt
  cell_type_specific_trans_regulatory_{celltype}.txt

CLI, not `script:` — 2026-07-25 fix, see linger_prep_pseudobulk.py docstring
for why (conda: directory-path misclassification under real SGE execution).
"""
import argparse
import os
import anndata as ad
import scanpy as sc

import LingerGRN.LL_net as LL_net

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True)
p.add_argument("--grn-dir", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--celltype", required=True)
p.add_argument("--labeled", required=True)
p.add_argument("--atac-consensus", required=True, nargs="+")
p.add_argument("--sample-ids", required=True, nargs="+")
p.add_argument("--output-done", required=True)
args = p.parse_args()

WORKDIR = args.workdir
GRN_DIR = args.grn_dir
GENOME = args.genome
CELLTYPE = args.celltype

os.chdir(WORKDIR)

adata_RNA = sc.read_h5ad(args.labeled)
adata_RNA.obs["barcode"] = adata_RNA.obs_names
adata_RNA.obs["label"] = adata_RNA.obs["cell_type"].values
adata_RNA.var["gene_ids"] = adata_RNA.var_names

atac_list = []
for sid, path in zip(args.sample_ids, args.atac_consensus):
    a = sc.read_h5ad(path)
    a.obs_names = [f"{bc}-{sid}" for bc in a.obs_names]
    atac_list.append(a)
adata_ATAC = ad.concat(atac_list, join="outer")
adata_ATAC.var["gene_ids"] = adata_ATAC.var_names
shared = adata_RNA.obs_names.intersection(adata_ATAC.obs_names)
adata_RNA = adata_RNA[shared].copy()
adata_ATAC = adata_ATAC[shared].copy()
adata_ATAC.obs["barcode"] = adata_ATAC.obs_names
adata_ATAC.obs["label"] = adata_RNA.obs["label"].values

print(f"--- {CELLTYPE} ---")
print("cell_type_specific_TF_RE_binding()...")
LL_net.cell_type_specific_TF_RE_binding(
    GRN_DIR, adata_RNA, adata_ATAC, GENOME, CELLTYPE, WORKDIR + "/", "scNN"
)

print("cell_type_specific_cis_reg()...")
LL_net.cell_type_specific_cis_reg(GRN_DIR, adata_RNA, adata_ATAC, GENOME, "scNN", WORKDIR + "/")

print("cell_type_specific_trans_reg()...")
LL_net.cell_type_specific_trans_reg(GRN_DIR, adata_RNA, CELLTYPE, WORKDIR + "/")

with open(args.output_done, "w") as f:
    f.write(f"celltype: {CELLTYPE}\n")

print(f"\nDone: {CELLTYPE} -> cell_type_specific_cis/trans_regulatory_{CELLTYPE}.txt")
