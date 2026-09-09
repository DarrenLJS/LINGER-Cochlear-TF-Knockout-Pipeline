"""
aggregate_tf_activity.py — Module 7 step 2, stack per-dataset regulon scores.

Pure pandas, deliberately not using LINGER_PYTHON's heavier env (runs in
snakemake_eddie's own python via PATH — see resources block; no torch/scanpy
needed here, just pandas).
"""
import argparse
import glob
import json
import os

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--module7-dir", required=True)
p.add_argument("--roles-json", required=True)
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

roles = json.loads(args.roles_json)

series = {}
for path in sorted(glob.glob(os.path.join(args.module7_dir, "*.regulon.tsv"))):
    entry_id = os.path.basename(path).replace(".regulon.tsv", "")
    df = pd.read_csv(path, sep="\t", index_col=0)
    series[entry_id] = df.iloc[:, 0]

wide = pd.DataFrame(series)
wide = wide.T
wide.insert(0, "role", [roles.get(i, "unknown") for i in wide.index])
wide.to_csv(args.output_tsv, sep="\t")

n_bulk = (wide["role"] == "bulk").sum()
n_baseline = (wide["role"] == "baseline").sum()
print(f"Wrote {args.output_tsv}: {n_bulk} bulk + {n_baseline} baseline datasets x {wide.shape[1]-1} TFs")
