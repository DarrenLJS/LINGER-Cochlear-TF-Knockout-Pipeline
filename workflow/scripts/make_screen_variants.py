"""
make_screen_variants.py — Module 8c, step 1: the variant list for one chunk of
the all-TF screen.

TFs are taken from the persisted population Exp.tsv (the exact 586-TF list the
nets were trained against) at run time, so the DAG needs no TF list: chunk i of
n gets the TFs at positions i, i+n, i+2n, ... (strided, so every chunk has a
similar mix). Every TF in the chunk gets one single-TF variant per requested
mode; variant ids are "<mode>__<TF>". Written as the JSON that
linger_perturbation.py --variants-json reads.
"""
import argparse
import json

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--exp", required=True, help="module8_perturbation/Exp.tsv")
p.add_argument("--chunk", type=int, required=True)
p.add_argument("--n-chunks", type=int, required=True)
p.add_argument("--modes", nargs="+", required=True, choices=["ko", "oe"])
p.add_argument("--out", required=True)
args = p.parse_args()

if not (0 <= args.chunk < args.n_chunks):
    raise ValueError(f"chunk {args.chunk} outside [0, {args.n_chunks})")
tfs = pd.read_csv(args.exp, sep="\t", index_col=0, usecols=[0]).index.tolist()
mine = [t for i, t in enumerate(tfs) if i % args.n_chunks == args.chunk]
if not mine:
    raise ValueError(f"chunk {args.chunk} of {args.n_chunks} is empty ({len(tfs)} TFs) — lower n_chunks")
variants = [{"id": f"{m}__{t}", "tfs": [t], "mode": m} for t in mine for m in args.modes]
with open(args.out, "w") as fh:
    json.dump(variants, fh)
print(f"chunk {args.chunk}/{args.n_chunks}: {len(mine)} TFs x {len(args.modes)} mode(s) = {len(variants)} variants -> {args.out}")
