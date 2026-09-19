"""
compute_tf_activity_shift.py — Module 11 step 1 (CAPS), TF-activity-space
build for `sign_adjusted_shift`.

Runs ONE forward-pass-output file (either the pre-KO sanity-check baseline
or a post-KO `{ko_id}_predicted_expression.tsv`, from Module 8 or 8b) through
the exact same `TF_activity.regulon()` call Module 9/9b already trust, and
writes a single collapsed TF-activity vector. compute_caps_score.py then
diffs two of this script's outputs (baseline vs a knockout) to get the
network-wide TF-activity shift a knockout causes.

WHY THE SANITY-CHECK FILE, NOT REAL DATA, AS THE "PRE-KO" BASELINE: both
`_sanity_check_baseline_predicted.tsv` and `{ko_id}_predicted_expression.tsv`
are outputs of the IDENTICAL forward-pass machinery (linger_perturbation.py,
same net, same live per-call normalization) — the only difference between
them is whether a TF row was zeroed first. Diffing predicted-vs-predicted
isolates the knockout's effect. Diffing predicted-vs-real (TG_pseudobulk)
would also fold the forward pass's own predicted-vs-real gap (the sanity_rho
numbers, 0.077-0.797 across scopes) into `sign_adjusted_shift` a second time
— that gap is already captured once, on purpose, via `confidence_weight`
(Module 9b's own sanity_rho column). Folding it into both terms would
double-count the same limitation instead of separating "how big is the
knockout's effect" from "how much do we trust this scope's forward pass".

ANNDATA WRAPPING: regulon()'s real signature (confirmed by reading LingerGRN
1.110's TF_activity.py source directly, same source Module 9/9b's aging
checks were designed against) needs adata_RNA.X (dense-able), .var["gene_ids"]
and .obs["barcode"]. `{ko_id}_predicted_expression.tsv` /
`_sanity_check_baseline_predicted.tsv` are genes(rows) x samples(cols) TSVs
(same shape aggregate_validation_celltype.py's _sanity_rho already reads) —
wrapped here with EXACTLY the same convention bulk_rna_loader.py's own
_adata() helper and Module 6's grn_celltype_specific.py already use:
adata.var["gene_ids"] = adata.var_names, adata.obs["barcode"] = adata.obs_names.

NETWORK SELECTION: population scope -> network="cell population" (reads
cell_population_trans_regulatory.txt); a cell-type scope -> network=<celltype>
(reads Module 6's cell_type_specific_trans_regulatory_{celltype}.txt) —
identical selection already exercised by Module 9b's aging rules.

COLLAPSE TO ONE VECTOR: regulon() returns TFs x samples. Row-mean across
samples, same convention validate_held_out.py's expression_shift mode uses
for the analogous gene-space collapse (adata.X.mean(axis=0)).
"""
import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import LingerGRN.TF_activity as TF_activity

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True, help="Module 4 workdir — also regulon()'s `outdir` (trans-regulatory files live here)")
p.add_argument("--grn-dir", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--network", required=True, help='"cell population" for the population scope, or a cell-type name')
p.add_argument("--expression-tsv", required=True, help="genes(rows) x samples(cols) TSV — sanity-check baseline or a knockout's predicted_expression.tsv")
p.add_argument("--output-tsv", required=True, help="TF-indexed, single 'tf_activity' column")
args = p.parse_args()

expr = pd.read_csv(args.expression_tsv, sep="\t", index_col=0)
print(f"read {args.expression_tsv}: {expr.shape[0]} genes x {expr.shape[1]} samples")

adata = ad.AnnData(X=sp.csr_matrix(expr.values.T))
adata.var_names = [str(g) for g in expr.index]
adata.obs_names = [str(s) for s in expr.columns]
adata.var["gene_ids"] = adata.var_names
adata.obs["barcode"] = adata.obs_names

os.chdir(args.workdir)
regulon_scores = TF_activity.regulon(args.workdir + "/", adata, args.grn_dir, args.network, args.genome)
print(f"regulon network: {args.network!r} -> {regulon_scores.shape[0]} TFs x {regulon_scores.shape[1]} samples")

tf_activity = regulon_scores.mean(axis=1)
tf_activity.index.name = "TF"
tf_activity.rename("tf_activity").to_frame().to_csv(args.output_tsv, sep="\t")
print(f"Wrote {args.output_tsv}: {len(tf_activity)} TFs")
