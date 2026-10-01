"""
summarise_tf_screen.py — Module 8c, step 3: one status table per screened scope.

Concatenates every chunk's status TSV (variant, mode, tfs, n_samples_changed,
shift_l2, shift_max_abs, noop) and adds, per single-TF variant:
  tf               the TF name
  tf_in_network    the TF has a non-zero column in this scope's trans-regulatory
                   matrix (population: cell_population_trans_regulatory.txt; a
                   cell type: cell_type_specific_trans_regulatory_<ct>.txt).
                   False => its TF activity is undefined (e.g. the 8 all-zero
                   population columns, Thap11 in the cell-type matrices).
  usable           (not noop) and tf_in_network — the variants the empirical
                   null and any ranking should use.
Fail-loud: the number of chunk files must equal --n-chunks and every chunk must
have the same variant-id prefix scheme.
"""
import argparse
import glob
import os

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--screen-dir", required=True, help="module8c_tf_screen/<scope>")
p.add_argument("--n-chunks", type=int, required=True)
p.add_argument("--trans-regulatory-tsv", required=True)
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

files = sorted(glob.glob(os.path.join(args.screen_dir, "chunk_*_status.tsv")))
if len(files) != args.n_chunks:
    raise ValueError(f"{args.screen_dir}: found {len(files)} chunk status files, expected {args.n_chunks}")
st = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
if st["variant"].duplicated().any():
    raise ValueError("duplicate variant ids across chunks")
st["tf"] = st["tfs"].astype(str)
tr = pd.read_csv(args.trans_regulatory_tsv, sep="\t", index_col=0)
colsum = tr.sum(axis=0)
st["tf_in_network"] = st["tf"].map(lambda t: bool(t in colsum.index and float(colsum[t]) > 0))
st["usable"] = (~st["noop"].astype(bool)) & st["tf_in_network"]
st.to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: {len(st)} variants; "
      f"{int(st['noop'].sum())} no-op, {int((~st['tf_in_network']).sum())} not in network, "
      f"{int(st['usable'].sum())} usable")
print(st.groupby("mode")[["noop", "tf_in_network", "usable"]].sum().to_string())
