"""
prep_pseudobulk_target_celltype.py — Module 8b step 1 (cell-type-resolved
perturbation input).

WHY THIS SCRIPT EXISTS: linger_perturbation.py's forward pass currently only
ever runs against ONE population-pooled pseudobulk (Module 4's
data/TG_pseudobulk.tsv / data/RE_pseudobulk.tsv, built once from ALL cells
regardless of type). The net's weights are necessarily shared/population-
trained (that's what population_training produces), but nothing about the
net or its normalization requires the INPUT fed through it at inference
time to be population-pooled too — normalization is recomputed fresh from
whatever input is given at call time (see linger_perturbation.py's
docstring). This script builds the cell-type-restricted pseudobulk input
that lets Module 8b exploit that.

NOT a repeat of grn_celltype_specific.py (Module 6 step 2): that script
calls LINGER's cell_type_specific_TF_RE_binding()/_cis_reg()/_trans_reg(),
which compute correlation-based regulatory scores per cell type directly
from single-cell data — it never forward-passes anything through
{chr}_net.pt, and it doesn't produce a pseudobulk expression/accessibility
matrix. The cell-type-restricted TG_pseudobulk/RE_pseudobulk this script
produces doesn't exist anywhere else in the pipeline; it has to be built
fresh by calling the SAME LingerGRN.pseudo_bulk.pseudo_bulk() function
Module 4 uses, just on a cell-type-subsetted adata_RNA/adata_ATAC — same
pattern as linger_prep_pseudobulk.py, not prep_pseudobulk_target.py (which
only reads Module 4's already-written population files, it doesn't call
pseudo_bulk() itself).

MIN-CELL GUARD: pseudobulking a handful of cells produces noisy per-gene
means that don't represent a real cell-type expression profile. Below
--min-cells, this script exits with a clear message and writes no output
rather than silently producing a low-N pseudobulk that looks the same
shape as a reliable one downstream. The Snakemake rule checks for this
script's actual output files, so a skipped cell type simply doesn't
produce Exp_{celltype}.tsv and the corresponding knockout rule won't run
for it — not a crash, a deliberate gap.
"""
import argparse
import os
import sys

import anndata as ad
import pandas as pd
import scanpy as sc

import LingerGRN.pseudo_bulk as pseudo_bulk

p = argparse.ArgumentParser()
p.add_argument("--labeled", required=True)
p.add_argument("--atac-consensus", required=True, nargs="+")
p.add_argument("--sample-ids", required=True, nargs="+")
p.add_argument("--celltype", required=True)
p.add_argument("--min-cells", type=int, default=100,
                help="Skip (no output written) if fewer than this many cells match --celltype "
                     "after the shared-barcode intersection below.")
p.add_argument("--singlepseudobulk", type=int, default=0)
p.add_argument("--output-dir", required=True)
args = p.parse_args()

CELLTYPE = args.celltype
os.makedirs(args.output_dir, exist_ok=True)

print(f"Loading labeled RNA: {args.labeled}")
adata_RNA = sc.read_h5ad(args.labeled)
assert "cell_type" in adata_RNA.obs.columns, (
    "labeled.h5ad has no 'cell_type' column — Module 3's manual annotation step must run first."
)

print("Loading + concatenating per-sample atac_consensus.h5ad...")
atac_list = []
for sid, path in zip(args.sample_ids, args.atac_consensus):
    a = sc.read_h5ad(path)
    a.obs_names = [f"{bc}-{sid}" for bc in a.obs_names]
    atac_list.append(a)
adata_ATAC = ad.concat(atac_list, join="outer", index_unique=None)

shared = adata_RNA.obs_names.intersection(adata_ATAC.obs_names)
adata_RNA = adata_RNA[shared].copy()
adata_ATAC = adata_ATAC[shared].copy()

# --- restrict to this cell type only, BEFORE pseudobulking ---
mask = adata_RNA.obs["cell_type"] == CELLTYPE
n_cells = int(mask.sum())
print(f"--- {CELLTYPE}: {n_cells} cells (of {adata_RNA.n_obs} total shared cells) ---")
if n_cells < args.min_cells:
    print(f"SKIPPING {CELLTYPE}: {n_cells} cells < --min-cells {args.min_cells}. "
          f"No output written — this cell type's pseudobulk would not be reliable. "
          f"See this script's docstring (MIN-CELL GUARD).")
    sys.exit(0)

adata_RNA = adata_RNA[mask].copy()
adata_ATAC = adata_ATAC[mask.values].copy()

adata_RNA.obs["barcode"] = adata_RNA.obs_names
adata_RNA.obs["label"] = adata_RNA.obs["cell_type"].values
adata_RNA.var["gene_ids"] = adata_RNA.var_names
adata_RNA.raw = None  # pseudo_bulk() sets .raw itself, same as linger_prep_pseudobulk.py

adata_ATAC.obs["barcode"] = adata_ATAC.obs_names
adata_ATAC.obs["label"] = adata_RNA.obs["label"].values
adata_ATAC.var["gene_ids"] = adata_ATAC.var_names

print(f"Running LingerGRN.pseudo_bulk.pseudo_bulk() on {CELLTYPE}-only cells "
      f"(singlepseudobulk={args.singlepseudobulk})...")
TG_pseudobulk, RE_pseudobulk = pseudo_bulk.pseudo_bulk(adata_RNA, adata_ATAC, args.singlepseudobulk)

tg_out = os.path.join(args.output_dir, f"TG_pseudobulk_{CELLTYPE}.tsv")
re_out = os.path.join(args.output_dir, f"RE_pseudobulk_{CELLTYPE}.tsv")
TG_pseudobulk.to_csv(tg_out, sep=",")
RE_pseudobulk.to_csv(re_out, sep=",")

print(f"Wrote {tg_out} ({TG_pseudobulk.shape}) and {re_out} ({RE_pseudobulk.shape})")
