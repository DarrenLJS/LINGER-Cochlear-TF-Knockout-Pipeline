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

REWRITTEN 2026-09-30 — root-cause fix for the "every cell type fails its own
sanity check" incident (population sanity rho=0.805; all 14 cell types
0.09-0.33). Root cause, confirmed against the real LingerGRN==1.110 source
(LingerGRN/pseudo_bulk.py): pseudo_bulk() does not just average expression
— every call re-derives HVG selection, per-gene scaling stats, PCA, and the
KNN neighbor graph from whatever cells it's given (see its source, lines
39-59). The OLD version of this script called pseudo_bulk() a second time,
fresh, on a cell-type-ISOLATED subset of cells — so its HVGs/scaling/PCA/
neighbor-graph were fit purely within that one cell type, a structurally
different feature space than the one Module 4's population call fit (across
all cell types jointly) and that population_training's {chr}_net.pt files
were actually trained against. That mismatch, not sample size, is why every
single cell type failed uniformly regardless of how many cells it had.

THE FIX: don't call pseudo_bulk() again at all. pseudo_bulk() names each of
its output pseudobulk columns after a real, individual cell barcode (see
source: `TG_filter1 = pd.DataFrame(TG_filter1.T, columns=allindex, ...)`,
where `allindex` is built by sampling barcodes out of each cell-type
cluster during the POPULATION run — pseudo_bulk.py lines 64-76). That means
Module 4's data/TG_pseudobulk.tsv / data/RE_pseudobulk.tsv already contain
cell-type-resolved pseudobulk columns; we just weren't reading them that
way. This script now:
  1. Loads labeled.h5ad purely to map each population pseudobulk column's
     barcode back to its real cell_type (no pseudobulking here at all).
  2. Loads Module 4's already-written population TG_pseudobulk.tsv /
     RE_pseudobulk.tsv directly (same files prep_pseudobulk_target.py's
     Target/Opn come from).
  3. Selects the subset of population columns whose barcode maps to
     --celltype, and writes those out unchanged as this cell type's
     TG_pseudobulk_{celltype}.tsv / RE_pseudobulk_{celltype}.tsv.
This guarantees exact feature-space/normalization-basis compatibility with
the trained net (it's a literal column subset of the exact matrix the net
was trained against), and needs no LingerGRN/Module 6 changes. Verified
empirically before this fix (2026-09-30): all 394 population pseudobulk
column barcodes resolved against labeled.h5ad with zero ambiguity, and the
resulting per-celltype column counts (e.g. Chondrocyte 29, Non-sensory
epithelium 47, Outer hair cell 11) reproduce pseudo_bulk()'s own
floor(sqrt(n_cells))+1 formula exactly — confirming this was always latent
in Module 4's output, just never read out per cell type.

Downstream contract is UNCHANGED: this script still writes
TG_pseudobulk_{celltype}.tsv / RE_pseudobulk_{celltype}.tsv in
--output-dir, and linger_perturbation.py's --target-path/--opn-path
consumption of them (including its missing-TF/missing-RE fill-with-0
logic — some genes/REs can still legitimately be absent from a cell
type's column subset) needs no changes at all.

MIN-CELL GUARD: unchanged in spirit — below --min-cells REAL cells of this
type in labeled.h5ad, this script exits with a clear message and writes no
output, rather than silently producing an unreliable low-N pseudobulk. The
Snakemake rule checks for this script's actual output files, so a skipped
cell type simply doesn't produce output and the corresponding knockout
rule won't run for it — not a crash, a deliberate gap.
"""
import argparse
import os
import sys

import pandas as pd
import scanpy as sc

p = argparse.ArgumentParser()
p.add_argument("--labeled", required=True,
               help="Module 3's labeled.h5ad — used ONLY to map population pseudobulk "
                    "column barcodes back to their real cell_type. No pseudobulking happens here.")
p.add_argument("--population-pseudobulk-dir", required=True,
               help="Dir holding Module 4's already-written data/TG_pseudobulk.tsv and "
                    "data/RE_pseudobulk.tsv (i.e. the module4_linger_init workdir).")
p.add_argument("--celltype", required=True)
p.add_argument("--min-cells", type=int, default=100,
               help="Skip (no output written) if fewer than this many REAL cells of "
                    "--celltype exist in labeled.h5ad.")
p.add_argument("--output-dir", required=True)
args = p.parse_args()

CELLTYPE = args.celltype
os.makedirs(args.output_dir, exist_ok=True)

print(f"Loading labeled RNA (barcode -> cell_type lookup only): {args.labeled}")
adata_RNA = sc.read_h5ad(args.labeled)
assert "cell_type" in adata_RNA.obs.columns, (
    "labeled.h5ad has no 'cell_type' column — Module 3's manual annotation step must run first."
)

n_cells = int((adata_RNA.obs["cell_type"] == CELLTYPE).sum())
print(f"--- {CELLTYPE}: {n_cells} real cells in labeled.h5ad ---")
if n_cells < args.min_cells:
    print(f"SKIPPING {CELLTYPE}: {n_cells} cells < --min-cells {args.min_cells}. "
          f"No output written — this cell type's pseudobulk would not be reliable. "
          f"See this script's docstring (MIN-CELL GUARD).")
    sys.exit(0)

barcode_to_celltype = adata_RNA.obs["cell_type"]

tg_pop_path = os.path.join(args.population_pseudobulk_dir, "data", "TG_pseudobulk.tsv")
re_pop_path = os.path.join(args.population_pseudobulk_dir, "data", "RE_pseudobulk.tsv")
print(f"Loading population pseudobulk (Module 4 output, not recomputed): {tg_pop_path}, {re_pop_path}")
TG_pop = pd.read_csv(tg_pop_path, sep=",", header=0, index_col=0)
RE_pop = pd.read_csv(re_pop_path, sep=",", header=0, index_col=0)
print(f"Population pseudobulk: {TG_pop.shape[1]} columns total")

pop_cols = pd.Index(TG_pop.columns)
resolved = pop_cols.intersection(barcode_to_celltype.index)
unresolved = pop_cols.difference(barcode_to_celltype.index)
if len(unresolved):
    print(f"WARNING: {len(unresolved)} of {len(pop_cols)} population pseudobulk column "
          f"barcodes were NOT found in labeled.h5ad's obs_names — these columns cannot be "
          f"attributed to any cell type and are excluded from every cell type's output. "
          f"First few: {list(unresolved[:5])}")

selected_cols = [c for c in resolved if barcode_to_celltype.loc[c] == CELLTYPE]
print(f"{CELLTYPE}: {len(selected_cols)} population pseudobulk columns resolve to this cell type")

if not selected_cols:
    print(f"SKIPPING {CELLTYPE}: 0 population pseudobulk columns resolved to this cell type "
          f"despite {n_cells} real cells passing --min-cells — labeled.h5ad's barcodes may not "
          f"match Module 4's pseudobulk column barcodes (e.g. different sample-suffixing). "
          f"No output written; inspect this mismatch before re-running.")
    sys.exit(0)

TG_celltype = TG_pop[selected_cols]
RE_celltype = RE_pop[selected_cols]

tg_out = os.path.join(args.output_dir, f"TG_pseudobulk_{CELLTYPE}.tsv")
re_out = os.path.join(args.output_dir, f"RE_pseudobulk_{CELLTYPE}.tsv")
TG_celltype.to_csv(tg_out, sep=",")
RE_celltype.to_csv(re_out, sep=",")

print(f"Wrote {tg_out} ({TG_celltype.shape}) and {re_out} ({RE_celltype.shape}) — "
      f"selected from Module 4's population pseudobulk, same feature space/normalization "
      f"basis the trained nets saw, not re-derived from scratch.")
