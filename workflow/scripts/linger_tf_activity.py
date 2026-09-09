"""
linger_tf_activity.py — Module 7, one dataset's regulon/TF-activity score.

Real LingerGRN==1.110 call (TF_activity.py, confirmed against the published
wheel):

  TF_activity.regulon(outdir, adata_RNA, GRNdir, network, genome)

network="cell population" for every dataset here — NOT "general" (that path
is hardcoded to hg19/hg38, see 07_bulk_tf_activity.smk header). Requires
Module 6 to have written cell_population_trans_regulatory.txt (a small
addition to grn_population_training.py — see README_MODULES_7_10.md).

CLI, not `script:` — same reason as every Module 4/6 script: LINGER_PYTHON is
invoked by absolute path, bypassing Snakemake's `conda:` directive.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from bulk_rna_loader import load_rna_only  # noqa: E402

import LingerGRN.TF_activity as TF_activity

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True)
p.add_argument("--grn-dir", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--entry-json", required=True, help="JSON-encoded config entry for this dataset")
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

os.chdir(args.workdir)

entry = json.loads(args.entry_json)
entry_id = entry.get("sample_id") or entry.get("ref_id")
print(f"--- {entry_id} (role={entry.get('role')}) ---")

adata_RNA = load_rna_only(entry)
print(f"Loaded {adata_RNA.n_obs} samples/cells x {adata_RNA.n_vars} genes")

score = TF_activity.regulon(
    args.workdir + "/", adata_RNA, args.grn_dir, "cell population", args.genome
)
# score: TFs (rows) x samples/cells (cols). Collapse to one value per TF by
# taking the mean across columns — appropriate for bulk RNA-seq (few/no
# replicates per condition to preserve) and for the baseline references
# (population-level "unperturbed" TF activity is exactly what Module 9 wants,
# not a per-cell breakdown these datasets don't have anyway).
mean_score = score.mean(axis=1)
mean_score.name = entry_id
mean_score.to_frame().to_csv(args.output_tsv, sep="\t")

print(f"Wrote {args.output_tsv} ({len(mean_score)} TFs)")
