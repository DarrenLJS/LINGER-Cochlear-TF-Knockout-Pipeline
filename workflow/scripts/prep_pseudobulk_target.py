"""
prep_pseudobulk_target.py — Module 8 step 1 (bypass design).

WHY THIS SCRIPT EXISTS (bypass over shim — see README's Module 8 design
notes for the full comparison against perturb.py's real, broken API):

  perturb.py's load_data_ptb()/get_simulation() were written for the
  human/hg19/hg38 'LINGER' training method's index/file scheme, confirmed
  from the real LingerGRN==1.110 source to be structurally incompatible
  with method='scNN' in three independent ways (hardcoded human chromosome
  list, a Target/data_merge row-count assumption scNN's filtered gene set
  violates, and an explicit-index TF/RE gather scheme that doesn't match
  scNN training's implicit "all TFs except self" scheme). Reimplementing
  the forward pass directly against {chr}_net.pt — modeled on scNN's own
  real training function LINGER_tr.sc_nn_NN(), not perturb.py's — sidesteps
  all three rather than reverse-engineering scNN's real contract into a
  file format built for a different one.

WHAT THIS SCRIPT DOES, and why each piece matters for correctness:

  1. Loads RE_TGlink FROM DISK (outdir/RE_TGlink.txt), NOT by recomputing
     it via LINGER_tr.load_data_scNN(). This is a deliberate, important
     choice: LINGER_tr.training()'s real scNN branch (confirmed from
     source) writes this exact file at the end of its Parallel() loop —
     it IS the real gene/chromosome ordering training used to key each
     {chr}_net.pt dict (by row position within that chromosome's slice,
     post reset_index). Recomputing it fresh here would only be needed if
     the file didn't exist; reading the real artifact instead means the
     per-gene network lookup below is exactly right, not "should match."

  2. Reconstructs Exp (the TF-only pseudobulk expression matrix) via the
     SAME set-intersection LINGER_tr.load_data_scNN() itself uses
     (TFlist = list(set(Target.index) & set(TFName[0].values))).
     FLAGGED RISK, not silently assumed safe: Python's set() iteration
     order depends on per-process string-hash randomization unless
     PYTHONHASHSEED is fixed. The ORIGINAL population_training run that
     produced {chr}_net.pt used some TFlist order that was never persisted
     to disk anywhere (scNN mode writes no TFName.txt/index.txt, unlike
     the human 'LINGER' path). This script's reconstructed order is
     therefore NOT guaranteed to match the order the net's fc1 weights
     were actually trained against — if it doesn't, forward-pass
     predictions will be silently wrong (right shape, meaningless values,
     no error thrown). This is exactly what linger_perturbation.py's
     mandatory --sanity-check pass is for (see that script's docstring) —
     run it and read its output before trusting any knockout result.

  3. Materializes Exp ONCE to outdir_module8/Exp.tsv, so every downstream
     knockout run (and the sanity-check baseline run) reads the SAME
     persisted order rather than re-deriving it per-process — this
     doesn't fix the risk in point 2, but it does guarantee every Module 8
     output is at least internally self-consistent with each other.

  Target and Opn are NOT re-derived here — they're the untouched
  population-level pseudobulk files Module 4 already wrote
  (data/TG_pseudobulk.tsv, data/RE_pseudobulk.tsv), read and indexed by
  label (gene/RE string), which is order-independent and therefore safe.
"""
import argparse
import os

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True, help="Module 4/6 workdir — same dir RE_TGlink.txt/TG_pseudobulk.tsv live in")
p.add_argument("--grn-dir", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--output-dir", required=True, help="Module 8 dir to write Exp.tsv/RE_TGlink_resolved.tsv into")
args = p.parse_args()

WORKDIR = args.workdir
GRN_DIR = args.grn_dir.rstrip("/") + "/"  # see grn_population_training.py's identical note on why this matters
GENOME = args.genome

os.makedirs(args.output_dir, exist_ok=True)
os.chdir(WORKDIR)

# --- species lookup, same pattern as grn_population_training.py ---
genome_map = pd.read_csv(os.path.join(GRN_DIR, "genome_map_homer.txt"), sep="\t")
genome_map.index = genome_map["genome_short"]
species = genome_map.loc[GENOME]["species_ensembl"]
print(f"genome={GENOME} -> species_ensembl={species}")

# --- Target: population pseudobulk, unfiltered, label-indexed (safe) ---
Target = pd.read_csv("data/TG_pseudobulk.tsv", sep=",", header=0, index_col=0)
print(f"Target (TG pseudobulk): {Target.shape}")

# --- Opn: population pseudobulk RE accessibility, unfiltered, label-indexed (safe) ---
Opn = pd.read_csv("data/RE_pseudobulk.tsv", sep=",", header=0, index_col=0)
print(f"Opn (RE pseudobulk): {Opn.shape}")

# --- RE_TGlink: read the REAL training artifact, not recomputed ---
re_tglink_path = os.path.join(WORKDIR, "RE_TGlink.txt")
if not os.path.exists(re_tglink_path):
    raise FileNotFoundError(
        f"{re_tglink_path} not found. This file is written by "
        f"LINGER_tr.training(method='scNN') at the end of Module 6's "
        f"linger_population_training rule — confirm that rule completed "
        f"successfully before running Module 8."
    )
re_tglink_raw = pd.read_csv(re_tglink_path, sep="\t", header=0)
# Columns confirmed from LINGER_tr.load_data_scNN()'s construction:
#   RE_TGlink = RE_TGlink.groupby('gene').apply(lambda x: x['RE'].values.tolist()).reset_index()
#   RE_TGlink['chr'] = [...]
# i.e. ['gene', '0', 'chr'] on disk — the middle column (pandas' default
# name for an unnamed .apply() aggregation result) holds each gene's list
# of RE strings, but to_csv() serializes that Python list as its str()
# repr ("['chr1:100-200', 'chr1:300-400']"), not real CSV-list syntax.
# NOT VALIDATED against your actual on-disk RE_TGlink.txt — inspect a few
# rows (`head -3 RE_TGlink.txt`) before trusting this parse if it errors.
import ast
re_col = [c for c in re_tglink_raw.columns if c not in ("gene", "chr")]
if len(re_col) != 1:
    raise ValueError(
        f"Expected exactly one RE-list column in RE_TGlink.txt besides "
        f"'gene'/'chr', found {re_col}. Inspect the file's real header "
        f"before proceeding — the groupby().apply() column-naming "
        f"convention this script assumes may not hold for your pandas "
        f"version."
    )
re_col = re_col[0]
re_tglink_raw[re_col] = re_tglink_raw[re_col].apply(ast.literal_eval)
print(f"RE_TGlink (real training artifact): {re_tglink_raw.shape[0]} genes across "
      f"{re_tglink_raw['chr'].nunique()} chromosomes")

# --- Exp: TF-only pseudobulk subset, reconstructed via LINGER's own logic ---
# CLI, not the `LINGER_tr.load_data_scNN()` convenience function, since we
# already have Target/RE_TGlink from the safer sources above — only the
# TFlist intersection logic is reused here, matching source exactly:
#   TFlist = list(set(Target.index) & set(TFName[0].values))
#   Exp = Target.loc[TFlist]
if species == "New":
    match2 = pd.read_csv(os.path.join(GRN_DIR, "MotifMatch.txt"), header=0, sep="\t")
else:
    match2 = pd.read_csv(os.path.join(GRN_DIR, f"Match_TF_motif_{species}.txt"), header=None, sep="\t")
    match2.columns = ["Motif", "TF"]
TFName = pd.DataFrame(match2["TF"].unique())
TFlist = list(set(Target.index) & set(TFName[0].values))
print(f"TFlist: {len(TFlist)} TFs — RISK: this order is a fresh set() "
      f"intersection, not necessarily the order population_training's "
      f"{{chr}}_net.pt weights were trained against. See this script's "
      f"docstring point 2 and run linger_perturbation.py --sanity-check "
      f"before trusting knockout output.")
Exp = Target.loc[TFlist]

# --- Persist, once, for every downstream Module 8 run to share ---
exp_out = os.path.join(args.output_dir, "Exp.tsv")
re_tglink_out = os.path.join(args.output_dir, "RE_TGlink_resolved.tsv")
Exp.to_csv(exp_out, sep="\t")
re_tglink_raw.to_csv(re_tglink_out, sep="\t", index=False)

print(f"\nWrote {exp_out} ({Exp.shape}) and {re_tglink_out} ({re_tglink_raw.shape[0]} genes) "
      f"— both consumed as-is by every linger_perturbation.py run, so all "
      f"Module 8 outputs share one consistent TF ordering.")
