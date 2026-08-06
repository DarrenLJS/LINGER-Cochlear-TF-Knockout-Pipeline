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
GRN_DIR = args.grn_dir
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
