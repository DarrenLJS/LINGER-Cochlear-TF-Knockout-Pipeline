"""
grn_population_training.py — Module 6 step 1 (population-level).

Two real LingerGRN==1.110 calls, confirmed from source:

  LINGER_tr.training(GRNdir, method='scNN', outdir, activef, species)
    -> per-chromosome neural net (SHAP-explained) trained per RE-TG link,
       parallelised across chromosomes via joblib. Writes
       {chr}_net.pt / {chr}_shap.pt / {chr}_Loss.txt / RE_TGlink.txt to
       outdir. `species` comes from GRNdir/genome_map_homer.txt, keyed by
       genome — this is the exact lookup the scNN.md tutorial does, not an
       assumption.

  LL_net.TF_RE_binding(GRNdir, adata_RNA, adata_ATAC, genome, method='scNN', outdir)
    -> reads MotifTarget.bed (Module 4's HOMER scan output) via
       load_TFbinding_scNN(), writes cell_population_TF_RE_binding.txt.

Both read/write relative to `outdir` (unlike get_TSS/RE_TG_dis in Module 4,
which use hardcoded `./data/` — see linger_get_tss.py docstring). outdir
here IS --workdir, kept the same directory as Module 4's pseudobulk/TSS
work so RE_TGlink.txt and data/*.tsv are both visible.

CLI, not `script:` — 2026-07-25 fix, see linger_prep_pseudobulk.py docstring
for why (conda: directory-path misclassification under real SGE execution).
"""
import argparse
import os
import shutil
import pandas as pd
import scanpy as sc
import anndata as ad

import LingerGRN.LINGER_tr as LINGER_tr
import LingerGRN.LL_net as LL_net

import torch
# FIX 2026-09-05, v2: add_safe_globals([LINGER_tr.Net]) (first attempt) got
# further but wasn't enough — PyTorch's weights_only=True unpickler walks
# the ENTIRE object graph inside a saved Net instance, not just the
# outermost class, and hit torch.nn.modules.linear.Linear next (real error,
# confirmed on Eddie). Allowlisting classes one at a time as they surface
# means as many retries as there are distinct layer types inside Net
# (Linear now, plausibly ReLU/Sequential/etc. next).
#
# Root cause: the installed LingerGRN version's torch.load() calls don't
# pass weights_only=False at all (confirmed: the current PyPI 1.110 release
# does, at every call site — this Eddie install predates that fix). Rather
# than replicate that fix call-by-call across LL_net.py/LINGER_tr.py/
# perturb.py (files we don't control, live in a shared conda env, would
# need reapplying after any env rebuild), monkeypatch torch.load itself so
# any call that doesn't explicitly pass weights_only defaults to False —
# matching the trust judgment the package's own later release already
# made for exactly these files (our own pipeline's freshly-generated
# checkpoints, not untrusted downloads).
_orig_torch_load = torch.load
def _torch_load_default_unsafe(*a, **kw):
    kw.setdefault("weights_only", False)
    return _orig_torch_load(*a, **kw)
torch.load = _torch_load_default_unsafe

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True)
p.add_argument("--grn-dir", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--activef", required=True)
p.add_argument("--motif-bed", required=True)
p.add_argument("--labeled", required=True)
p.add_argument("--atac-consensus", required=True, nargs="+")
p.add_argument("--sample-ids", required=True, nargs="+")
p.add_argument("--output-done", required=True)
args = p.parse_args()

WORKDIR = args.workdir
GRN_DIR = args.grn_dir.rstrip("/") + "/"
# FIX 2026-09-05: LingerGRN's own internal functions (LINGER_tr.load_data_scNN,
# confirmed at LINGER_tr.py line 323) build paths as raw string concatenation
# — `GRNdir + 'Match_TF_motif_' + species + '.txt'` — with NO separator
# inserted. They assume the caller already passed GRNdir ending in '/'.
# Passing --grn-dir without a trailing slash (as config_eddie.yaml's
# linger.grn_dir currently does) produced a real crash on the first actual
# run: FileNotFoundError on '.../provide_dataMatch_TF_motif_Mus_musculus.txt'
# (note the missing '/' between 'provide_data' and 'Match_TF_motif').
# os.path.join(GRN_DIR, "genome_map_homer.txt") below tolerates a missing
# trailing slash and is NOT the affected call — only LINGER's own internal
# concatenation is. Normalizing once here, rather than at every call site,
# so this can't regress if another LINGER_tr/LL_net call is added later.
GENOME = args.genome
ACTIVEF = args.activef

os.makedirs(WORKDIR, exist_ok=True)
os.chdir(WORKDIR)

# MotifTarget.bed was written by Module 4's linger_motif_scan rule into
# MODULE4_DIR — LL_net.TF_RE_binding expects it at outdir+'MotifTarget.bed',
# i.e. right here, so copy it in rather than duplicating the HOMER scan.
motif_src = args.motif_bed
motif_dst = os.path.join(WORKDIR, "MotifTarget.bed")
if not os.path.exists(motif_dst):
    shutil.copy(motif_src, motif_dst)

genome_map = pd.read_csv(os.path.join(GRN_DIR, "genome_map_homer.txt"), sep="\t")
genome_map.index = genome_map["genome_short"]
species = genome_map.loc[GENOME]["species_ensembl"]
print(f"genome={GENOME} -> species_ensembl={species} (from GRNdir/genome_map_homer.txt)")

print(f"LINGER_tr.training(GRNdir, method='scNN', outdir={WORKDIR}, activef={ACTIVEF}, species={species})")
# RESUME GUARD, added 2026-09-05: training() has no resume logic of its own
# and always retrains every chromosome from scratch — expensive (4h08m on
# the last real run). If a prior run got all the way through training and
# only failed on the *next* step (TF_RE_binding, the weights_only bug this
# same patch fixes), re-deriving the same chromosome list training() itself
# would compute and checking whether every {chr}_net.pt/{chr}_shap.pt
# already exists lets us skip straight past it instead of burning another
# 4+ hours to re-reach a step that's already fixed.
re_tglink_path = os.path.join(WORKDIR, "data", "RE_gene_distance.txt")
chrlist = sorted(pd.read_csv(re_tglink_path, sep="\t")["RE"].str.split(":").str[0].unique())
already_done = all(
    os.path.exists(os.path.join(WORKDIR, f"{c}_net.pt")) and os.path.exists(os.path.join(WORKDIR, f"{c}_shap.pt"))
    for c in chrlist
)
if already_done:
    print(f"RESUME: all {len(chrlist)} chromosome {{net,shap}}.pt files already present in {WORKDIR} — "
          f"skipping LINGER_tr.training(), which has no partial-completion detection of its own.")
else:
    LINGER_tr.training(GRN_DIR, "scNN", WORKDIR + "/", ACTIVEF, species)

print("Loading labeled RNA + pooled ATAC for TF_RE_binding()...")
adata_RNA = sc.read_h5ad(args.labeled)
adata_RNA.obs["barcode"] = adata_RNA.obs_names
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

print("LL_net.TF_RE_binding(method='scNN')...")
LL_net.TF_RE_binding(GRN_DIR, adata_RNA, adata_ATAC, GENOME, "scNN", WORKDIR + "/")

with open(args.output_done, "w") as f:
    f.write(f"population training complete\nworkdir: {WORKDIR}\nspecies: {species}\n")

print(f"\nDone. cell_population_TF_RE_binding.txt + RE_TGlink.txt + per-chr .pt files under {WORKDIR}/")
