"""
build_coverage_table.py — Module 9c: coverage of the supervisor's xlsx test units by LINGER (ADDED 2026-10-17).

The supervisor's decks require coverage to be reported openly (a model that cannot legally take a test records N/A with the reason,
rather than silently dropping it). LINGER perturbs transcription factors only, so most xlsx units cannot be scored. This table lists
EVERY xlsx test unit once, with its LINGER status, the reason, and (where scored) the Module 9c verdict.

INPUTS
  resources/xlsx_test_units.tsv        one row per xlsx unit: gene, OE/KD class (as in the xlsx), xlsx dimensions, LINGER status + reason,
                                       linger_panel_ids (';'-separated panel_id values of Module 9c that score it)
  module9c_panel/panel_report.tsv      Module 9c report (panel_id, verdict, excess_rho, ...)

OUTPUTS (--out-dir)
  coverage_table.tsv     one row per xlsx unit: status, reason, verdict(s), excess rho of the scored arm(s)
  coverage_summary.tsv   (dimension x class) counts: units in the xlsx vs units LINGER can score, and the verdict tally
"""
import argparse
import os

import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--units-tsv", required=True)
p.add_argument("--panel-report", required=True)
p.add_argument("--out-dir", required=True)
args = p.parse_args()
os.makedirs(args.out_dir, exist_ok=True)

units = pd.read_csv(args.units_tsv, sep="\t", dtype=str).fillna("")
rep = pd.read_csv(args.panel_report, sep="\t")
rep = rep.set_index("panel_id") if "panel_id" in rep.columns else rep.assign(panel_id=rep.index).set_index("panel_id")

rows = []
for _, u in units.iterrows():
    pids = [x for x in u["linger_panel_ids"].split(";") if x]
    verd, rho, lo, hi, missing = [], [], [], [], []
    for pid in pids:
        if pid in rep.index:
            r = rep.loc[pid]
            verd.append(str(r.get("verdict", "")))
            rho.append(r.get("excess_rho", np.nan))
            lo.append(r.get("excess_ci_lo", np.nan))
            hi.append(r.get("excess_ci_hi", np.nan))
        else:
            missing.append(pid)
    status = u["linger_status"]
    if status.startswith("scored") and not verd:
        final = f"{status}: panel row missing ({','.join(missing)}) — rerun module9c_all"
    elif status.startswith("scored"):
        final = f"{status}: {'/'.join(verd)}"
    elif status == "epsilon_exploratory":
        final = "epsilon only: see metrics_epsilon.tsv"
    elif status == "tf_not_in_network":
        final = "N/A: TF not in trained network (descriptive only)"
    elif status == "not_a_tf":
        final = "N/A: not a transcription factor"
    elif status == "not_run_array":
        final = "N/A: array data not supported"
    else:
        final = status
    d = u.to_dict()
    d.update(final_status=final, verdicts=";".join(verd), excess_rho=";".join("" if pd.isna(x) else f"{x:.3f}" for x in rho),
             excess_ci_lo=";".join("" if pd.isna(x) else f"{x:.3f}" for x in lo),
             excess_ci_hi=";".join("" if pd.isna(x) else f"{x:.3f}" for x in hi))
    rows.append(d)
cov = pd.DataFrame(rows)
cov.to_csv(os.path.join(args.out_dir, "coverage_table.tsv"), sep="\t", index=False)

# (dimension x class) summary
srows = []
for dim in "1234567":
    for cls in ("OE", "KD"):
        sub = cov[cov["xlsx_dimensions"].str.split(",").apply(lambda L: dim in L) & (cov["perturb_class_xlsx"] == cls)]
        scored = sub[sub["linger_status"].str.startswith("scored")]
        vt = ";".join(f"{k}={v}" for k, v in
                      pd.Series([x.split(":")[-1].strip().split("/")[0] for x in scored["final_status"]]).value_counts().items())
        srows.append({"dimension": int(dim), "class_xlsx": cls, "n_xlsx_units": int(len(sub)), "n_linger_scored": int(len(scored)),
                      "n_tf_not_in_network": int((sub["linger_status"] == "tf_not_in_network").sum()),
                      "n_not_a_tf": int((sub["linger_status"] == "not_a_tf").sum()),
                      "n_other_na": int(sub["linger_status"].isin(["not_run_array", "epsilon_exploratory"]).sum()),
                      "verdicts_of_scored": vt})
pd.DataFrame(srows).to_csv(os.path.join(args.out_dir, "coverage_summary.tsv"), sep="\t", index=False)
print(f"coverage_table.tsv: {len(cov)} units; status counts:\n{cov['linger_status'].value_counts().to_string()}")
print(pd.DataFrame(srows).to_string(index=False))
