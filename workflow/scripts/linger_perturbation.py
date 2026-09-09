"""
linger_perturbation.py — Module 8 step 2 (bypass design, not perturb.py).

Reimplements the forward pass directly against Module 6's trained
{chr}_net.pt files, modeled EXACTLY on LINGER_tr.sc_nn_NN() — the real
function that TRAINED those nets (confirmed from source) — rather than on
perturb.py's LINGER_simulation()/get_simulation(), which were built for a
different (human/hg19/hg38 'LINGER' method) input-index scheme and are
confirmed structurally incompatible with scNN's real training contract in
three ways (see prep_pseudobulk_target.py's docstring and the pipeline
README's Module 8 design notes for the full comparison).

PER-GENE FORWARD PASS, replicated exactly from LINGER_tr.sc_nn_NN():
  TFtemp = Exp.drop([gene]).values   if gene is itself a TF, else Exp.values
  REtemp = Opn.loc[RE_TGlink[gene]].values
  inputs = vstack(TFtemp, REtemp)
  inputs = per-row z-score normalized (mean/std computed live from THIS
           call's inputs, across samples — dim=1) — this is also exactly
           how perturb.py's own LINGER_simulation() re-normalizes at
           inference time, confirmed from source, and it's the actual
           mechanism by which a knockout propagates: zeroing a TF's row
           makes its own normalized value ~0 (mean=0, std~eps) without
           touching any other row's normalization, since normalization is
           per-row/per-feature, not global.
  y_pred = net(inputs.T)             where net = netall[gene_row_index],
                                      netall = torch.load({chr}_net.pt)

KNOCKOUT MECHANICS: --tf-list zeroes the named TF row(s) in a COPY of the
persisted Exp (from prep_pseudobulk_target.py) before the forward pass —
nothing else changes. Multiple TFs (triple knockout) zero multiple rows
in the same copy.

--sanity-check MODE (run this FIRST, before trusting any real knockout):
runs the SAME forward pass with an UNMODIFIED Exp (no knockout) and
reports Spearman correlation between predicted and real (Target) values
for every gene that got a trained net. This is the only real check
available for the TF-ordering risk flagged in prep_pseudobulk_target.py's
docstring — a low correlation here means the reconstructed Exp ordering
almost certainly does NOT match what population_training's nets were
trained against, and knockout output from this pipeline should not be
trusted until that's resolved (e.g. checking whether PYTHONHASHSEED was
fixed in the LINGER conda env at training time).
"""
import argparse
import ast
import json
import os

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

# Same PyTorch 2.6+ weights_only fix as grn_population_training.py /
# grn_celltype_specific.py — these are our own freshly-trained checkpoints,
# not untrusted downloads, matching the trust judgment LingerGRN's own
# later releases already made for these exact call sites.
_orig_torch_load = torch.load
def _torch_load_default_unsafe(*a, **kw):
    kw.setdefault("weights_only", False)
    return _orig_torch_load(*a, **kw)
torch.load = _torch_load_default_unsafe

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True, help="Module 4/6 workdir, where {chr}_net.pt files live")
p.add_argument("--module8-dir", required=True, help="dir with Exp.tsv/RE_TGlink_resolved.tsv from prep_pseudobulk_target.py")
p.add_argument("--ko-id", required=True)
p.add_argument("--tf-list", nargs="*", default=[], help="TFs to zero. Empty list == sanity-check baseline (no knockout).")
p.add_argument("--sanity-check", action="store_true",
                help="Also compute predicted-vs-real correlation against Target — forces --tf-list to be ignored for scoring purposes (still applied to the forward pass if given).")
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

WORKDIR = args.workdir

Exp = pd.read_csv(os.path.join(args.module8_dir, "Exp.tsv"), sep="\t", index_col=0)
Opn = pd.read_csv(os.path.join(WORKDIR, "data", "RE_pseudobulk.tsv"), sep=",", header=0, index_col=0)
Target = pd.read_csv(os.path.join(WORKDIR, "data", "TG_pseudobulk.tsv"), sep=",", header=0, index_col=0)
RE_TGlink = pd.read_csv(os.path.join(args.module8_dir, "RE_TGlink_resolved.tsv"), sep="\t")
re_col = [c for c in RE_TGlink.columns if c not in ("gene", "chr")][0]
RE_TGlink[re_col] = RE_TGlink[re_col].apply(ast.literal_eval)

print(f"--- {args.ko_id} (knockout TFs: {args.tf_list or '[none — baseline]'}) ---")

# Apply knockout to a COPY — Exp on disk is shared across every ko_id run.
Exp_ko = Exp.copy()
missing = [tf for tf in args.tf_list if tf not in Exp_ko.index]
if missing:
    print(f"WARNING: {missing} not found in Exp's TF list — cannot knock out, skipping these")
for tf in args.tf_list:
    if tf in Exp_ko.index:
        Exp_ko.loc[tf] = 0

predicted = {}  # gene -> np.array of predicted values across samples
missing_nets = 0
eps = 1e-6

for chrom in sorted(RE_TGlink["chr"].unique()):
    net_path = os.path.join(WORKDIR, f"{chrom}_net.pt")
    if not os.path.exists(net_path):
        print(f"WARNING: {net_path} not found — skipping chromosome {chrom} entirely "
              f"({(RE_TGlink['chr'] == chrom).sum()} genes affected)")
        continue
    netall = torch.load(net_path)
    chr_rows = RE_TGlink[RE_TGlink["chr"] == chrom].reset_index(drop=True)

    for idx in range(chr_rows.shape[0]):
        if idx not in netall:
            missing_nets += 1
            continue  # gene had no net trained (sc_nn_NN's own good==0 path, if it ever applies here)
        gene = chr_rows.loc[idx, "gene"]
        re_list = chr_rows.loc[idx, re_col]

        if gene in Exp_ko.index:
            TFtemp = Exp_ko.drop([gene]).values
        else:
            TFtemp = Exp_ko.values
        REtemp = Opn.loc[re_list].values
        inputs = np.vstack((TFtemp, REtemp))
        inputs = torch.tensor(inputs, dtype=torch.float32)
        mean = inputs.mean(dim=1)
        std = inputs.std(dim=1)
        inputs = ((inputs.T - mean) / (std + eps)).T

        net = netall[idx]
        net.eval()
        with torch.no_grad():
            y_pred = net(inputs.T)
        predicted[gene] = y_pred.detach().numpy().reshape(-1)

if missing_nets:
    print(f"NOTE: {missing_nets} genes had no trained net for their chromosome "
          f"(sc_nn's own per-gene training-failure path) — excluded from output, not an error here.")

pred_df = pd.DataFrame(predicted).T
pred_df.columns = Target.columns
pred_df.to_csv(args.output_tsv, sep="\t")
print(f"Wrote {args.output_tsv}: {pred_df.shape[0]} genes x {pred_df.shape[1]} samples")

if args.sanity_check:
    common = [g for g in pred_df.index if g in Target.index]
    rhos = []
    for g in common:
        rho, _ = spearmanr(pred_df.loc[g], Target.loc[g])
        if not np.isnan(rho):
            rhos.append(rho)
    if rhos:
        rhos = np.array(rhos)
        print(f"\nSANITY CHECK — predicted vs real Target, {len(rhos)} genes with a trained net:")
        print(f"  median Spearman rho = {np.median(rhos):.3f}   "
              f"mean = {rhos.mean():.3f}   "
              f"frac rho>0.3 = {(rhos > 0.3).mean():.2f}")
        print(f"  INTERPRETATION: this should be clearly positive and not close to 0 for a "
              f"large majority of genes — these nets were fit to reproduce Target from this "
              f"exact Exp/Opn input, so a low correlation here most likely means the "
              f"reconstructed Exp/TF ordering (see prep_pseudobulk_target.py's docstring) "
              f"does NOT match what {{chr}}_net.pt was actually trained against. Do not trust "
              f"real knockout output from this pipeline until this reads clearly positive.")
    else:
        print("\nSANITY CHECK: no overlapping genes scored — check Target/predicted gene name alignment.")
