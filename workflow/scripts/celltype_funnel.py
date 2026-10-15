"""
celltype_funnel.py — where every annotated cluster / cell type drops out of the pipeline (ADDED 2026-10-15).

README and figures used to say "22 cell types" / "20 cell types" while different stages used different
subsets. This reads the pipeline's own files and writes one row per distinct cell-type label with the last
stage it reached, plus a stage-count summary:

  clusters_annotated   rows of cluster_annotation.tsv (Leiden clusters, ALL of them)
  labelled             distinct non-empty `cell_type` labels (several clusters can share one label)
  with_network         labels that have a cell-type-specific trans-regulatory matrix (Module 6)
  benchmarked          labels with rows in Module 10's benchmark_report.tsv (cis_RE_TG)
  perturbed            labels with Module 8b sanity-check output (min_cells_for_pseudobulk survivors)
  ko_evaluated         labels with at least one CAPS row of status not ko_tf_no_network (Module 11)

Every input is optional (a stage whose files are absent is reported as `NA`, not as 0).
"""
import argparse
import os
import re

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--annotation", required=True)
p.add_argument("--module4-dir", required=True, help="module4_linger_init (cell_type_specific_trans_regulatory_*.txt)")
p.add_argument("--benchmark-tsv", default=None)
p.add_argument("--module8b-dir", default=None)
p.add_argument("--caps-tsv", default=None)
p.add_argument("--output-tsv", required=True)
p.add_argument("--output-summary-tsv", required=True)
args = p.parse_args()

ann = pd.read_csv(args.annotation, sep="\t")
if "cell_type" not in ann.columns:
    raise ValueError(f"{args.annotation} has no cell_type column; columns: {list(ann.columns)}")
labels = sorted(ann["cell_type"].dropna().astype(str).str.strip().loc[lambda s: s != ""].unique())

nets = None
if os.path.isdir(args.module4_dir):
    nets = {re.sub(r"^cell_type_specific_trans_regulatory_|\.txt$", "", f)
            for f in os.listdir(args.module4_dir) if f.startswith("cell_type_specific_trans_regulatory_")}
bench = None
if args.benchmark_tsv and os.path.exists(args.benchmark_tsv):
    b = pd.read_csv(args.benchmark_tsv, sep="\t")
    bench = set(b["celltype"].astype(str)) if "celltype" in b.columns else None
pert = None
if args.module8b_dir and os.path.isdir(args.module8b_dir):
    pert = {d for d in os.listdir(args.module8b_dir)
            if os.path.exists(os.path.join(args.module8b_dir, d, "_sanity_check_baseline_predicted.tsv"))}
kos = None
if args.caps_tsv and os.path.exists(args.caps_tsv):
    c = pd.read_csv(args.caps_tsv, sep="\t")
    kos = set(c.loc[(c["scope"] != "population") & (c["status"] != "ko_tf_no_network"), "scope"].astype(str))


def _flag(s, x):
    return "NA" if s is None else bool(x in s)


rows = []
for lab in labels:
    r = {"cell_type": lab, "n_clusters": int((ann["cell_type"].astype(str).str.strip() == lab).sum()),
         "with_network": _flag(nets, lab), "benchmarked": _flag(bench, lab),
         "perturbed": _flag(pert, lab), "ko_evaluated": _flag(kos, lab)}
    last = "labelled"
    for st in ("with_network", "benchmarked", "perturbed", "ko_evaluated"):
        if r[st] is True:
            last = st
    r["last_stage"] = last
    rows.append(r)
out = pd.DataFrame(rows)
out.to_csv(args.output_tsv, sep="\t", index=False)


def _n(col):
    return "NA" if (out[col] == "NA").all() else int((out[col] == True).sum())  # noqa: E712


summ = pd.DataFrame([
    {"stage": "clusters_annotated", "n": int(len(ann))},
    {"stage": "labelled", "n": len(labels)},
    {"stage": "with_network", "n": _n("with_network")},
    {"stage": "benchmarked", "n": _n("benchmarked")},
    {"stage": "perturbed", "n": _n("perturbed")},
    {"stage": "ko_evaluated", "n": _n("ko_evaluated")},
])
summ.to_csv(args.output_summary_tsv, sep="\t", index=False)
print(summ.to_string(index=False))
dropped = out[out["last_stage"] != "ko_evaluated"]
if len(dropped):
    print("labels that did not reach ko_evaluated:\n" + dropped[["cell_type", "last_stage"]].to_string(index=False))
