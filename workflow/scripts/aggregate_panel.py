"""
aggregate_panel.py — Module 9c: combine the generalisation-panel score rows into panel_report.tsv.

One row per scored panel unit (a dataset, or one arm of GSE86204) with the usual validate_held_out.py columns
(spearman_rho, level_null_rho, excess_rho + CI, verdict, provenance md5s), joined to the panel metadata
(dataset, gene, arm, eligibility, exploratory). Entries whose gene is not in the trained network get a
`descriptive_only` row (no perturbation, no score) so the panel's coverage is visible in one table.
Adds `pass_fail` = verdict (PASS / FAIL / INCONCLUSIVE / NOT_EVALUABLE / descriptive_only).
"""
import argparse
import glob
import json
import os

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--panel-dir", required=True)
p.add_argument("--meta-json", required=True)
p.add_argument("--descriptive-json", required=True)
p.add_argument("--output-tsv", required=True)
args = p.parse_args()

meta = json.loads(args.meta_json)
rows = []
for path in sorted(glob.glob(os.path.join(args.panel_dir, "*_score.tsv"))):
    r = pd.read_csv(path, sep="\t")
    pid = os.path.basename(path)[: -len("_score.tsv")]
    m = meta.get(pid)
    if m is None:
        raise ValueError(f"{path}: panel id {pid!r} is not in the config-derived metadata")
    for k, v in m.items():
        r[k] = v
    r["panel_id"] = pid
    rows.append(r)
for d in json.loads(args.descriptive_json):
    d = dict(d)
    d.update({"panel_id": d["test_id"], "verdict": "descriptive_only",
              "note": "gene not in the trained network: not perturbed, not scored"})
    rows.append(pd.DataFrame([d]))
if not rows:
    raise ValueError(f"no panel rows found under {args.panel_dir}")
rep = pd.concat(rows, ignore_index=True)
rep["pass_fail"] = rep["verdict"]
first = ["panel_id", "test_id", "dataset", "gene", "arm", "eligibility", "exploratory", "ko_id", "verdict", "pass_fail",
         "spearman_rho", "level_null_rho", "excess_rho", "excess_ci_lo", "excess_ci_hi", "rho_ci_lo", "rho_ci_hi",
         "n_real_target", "n_real_reference", "pred_shift_sd"]
rep = rep[[c for c in first if c in rep.columns] + [c for c in rep.columns if c not in first]]
rep.to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: {len(rep)} rows")
print(rep[[c for c in ["panel_id", "gene", "arm", "ko_id", "verdict", "spearman_rho", "excess_rho"] if c in rep.columns]].to_string(index=False))
