"""
aggregate_validation.py — Module 9 step 2, combine per-dataset scores into
one report with a pass/fail column.

Pure pandas — same reasoning as aggregate_tf_activity.py for not needing
the heavier LINGER_PYTHON env here.

REVISED 2026-10-01 — PASS/FAIL now comes from validate_held_out.py's
`verdict` column (bootstrap CI of the excess Spearman rho over a level-only
null; see that script's REVISED docstring block). The old rule
(rho > 0 AND auroc > 0.5) is gone: it was satisfied by gene-level structure
alone (v1's predicted shift was ~99.9% a knockout-independent level term),
which is how the Tmie negative control "passed" at population level with an
AUROC of 0.346 while failing in all 14 cell types.
  verdict values: PASS | FAIL | INCONCLUSIVE | NOT_EVALUABLE | not_scored
  (not_scored: aging rows, and descriptive checks such as tbx2_conversion).
  For the negative control the verdict is already inverted inside
  validate_held_out.py (PASS = no reprogramming signal beyond the null).
AUROC and AUPR remain as reported columns only.
"""
import argparse
import glob
import os

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--module9-dir", required=True)
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

rows = []
for path in sorted(glob.glob(os.path.join(args.module9_dir, "*_score.tsv"))):
    rows.append(pd.read_csv(path, sep="\t"))
report = pd.concat(rows, ignore_index=True)


if "verdict" not in report.columns:
    raise ValueError("score files carry no 'verdict' column — they were written by the pre-2026-10-01 "
                     "validate_held_out.py; delete module9_validation/ and rerun module9_all.")
report["pass_fail"] = report["verdict"]
report.to_csv(args.output_tsv, sep="\t", index=False)

print(f"Wrote {args.output_tsv}:")
cols = [c for c in ["sample_id", "check", "spearman_rho", "level_null_rho", "excess_rho",
                    "excess_ci_lo", "excess_ci_hi", "auroc", "pass_fail"] if c in report.columns]
print(report[cols].to_string(index=False))
