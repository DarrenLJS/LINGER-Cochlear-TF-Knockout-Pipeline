"""
grn_celltype_specific.py — Module 6 step 2 (per-cell-type GRN).

Real LingerGRN==1.110 calls (LL_net.py), one cell type at a time, matching
docs/GRN_infer.md's documented usage pattern:

  LL_net.cell_type_specific_TF_RE_binding(GRNdir, adata_RNA, adata_ATAC,
                                           genome, celltype, outdir, method)
  LL_net.cell_type_specific_cis_reg(GRNdir, adata_RNA, adata_ATAC, genome,
                                     celltype, outdir, method)
                                     # FIXED 2026-09-07: this WAS documented
                                     # (wrongly) as taking no celltype arg,
                                     # with method/outdir in the last two
                                     # slots. Real source (confirmed via
                                     # TypeError on Eddie: "missing 1
                                     # required positional argument:
                                     # 'method'") is celltype at position 5,
                                     # outdir at 6, method LAST at 7 — the
                                     # original call silently passed "scNN"
                                     # as celltype and never supplied method
                                     # at all.
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

import torch
# FIX 2026-09-05, v2: same fix as grn_population_training.py — allowlisting
# LINGER_tr.Net alone wasn't enough (weights_only=True unpickling walks the
# whole object graph inside a saved Net, not just the outermost class; hit
# torch.nn.modules.linear.Linear next in the real Eddie run). Monkeypatching
# torch.load to default weights_only=False covers the whole graph in one
# shot rather than allowlisting classes one at a time as they surface — see
# grn_population_training.py's identical comment for the full explanation.
_orig_torch_load = torch.load
def _torch_load_default_unsafe(*a, **kw):
    kw.setdefault("weights_only", False)
    return _orig_torch_load(*a, **kw)
torch.load = _torch_load_default_unsafe

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
GRN_DIR = args.grn_dir.rstrip("/") + "/"
# FIX 2026-09-05: same bug as grn_population_training.py — LingerGRN's
# LL_net functions (cell_type_specific_TF_RE_binding/_cis_reg/_trans_reg)
# build GRNdir-relative paths via raw string concatenation with no
# separator inserted (confirmed pattern: LINGER_tr.load_data_scNN crashed
# on this exact class of bug at population-training time). Normalizing once
# here rather than trusting every call site downstream.
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
LL_net.cell_type_specific_cis_reg(GRN_DIR, adata_RNA, adata_ATAC, GENOME, CELLTYPE, WORKDIR + "/", "scNN")

print("cell_type_specific_trans_reg()...")
LL_net.cell_type_specific_trans_reg(GRN_DIR, adata_RNA, CELLTYPE, WORKDIR + "/")

with open(args.output_done, "w") as f:
    f.write(f"celltype: {CELLTYPE}\n")

print(f"\nDone: {CELLTYPE} -> cell_type_specific_cis/trans_regulatory_{CELLTYPE}.txt")
