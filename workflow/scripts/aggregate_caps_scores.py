"""
aggregate_caps_scores.py — Module 11 step 3, combine every per-(ko_id,
scope) CAPS row into one report.

EQUAL TREATMENT, same convention as aggregate_validation_celltype.py: every
row (population's own ko_id x population, and every cell type's ko_id x
that cell type) is concatenated under one identical schema with a `scope`
column, sorted by CAPS so the ranking Module 11 exists to produce is
directly visible — not split into a population table and a cell-type table.
"""
import argparse
import glob
import os

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--module11-dir", required=True, help="MODULE11_DIR — holds {scope}/{ko_id}_caps.tsv")
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

rows = []
for path in sorted(glob.glob(os.path.join(args.module11_dir, "*", "*_caps.tsv"))):
    rows.append(pd.read_csv(path, sep="\t"))

if rows:
    report = pd.concat(rows, ignore_index=True)
    report = report.sort_values("CAPS", ascending=False).reset_index(drop=True)
else:
    report = pd.DataFrame(columns=[
        "scope", "ko_id", "tfs", "sign_adjusted_shift", "n_common_aging_tfs",
        "regulatory_centrality", "sanity_rho", "confidence_weight", "CAPS",
    ])
    print("WARNING: no Module 11 per-(ko_id, scope) CAPS files found — writing empty report")

report.to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: {len(report)} rows, "
      f"{report['scope'].nunique() if len(report) else 0} scopes "
      f"(population + cell types), {report['ko_id'].nunique() if len(report) else 0} ko_ids")
print(report[["scope", "ko_id", "sign_adjusted_shift", "regulatory_centrality",
              "confidence_weight", "CAPS"]].to_string(index=False))
