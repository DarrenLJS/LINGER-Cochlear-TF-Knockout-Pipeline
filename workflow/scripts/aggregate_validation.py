"""
aggregate_validation.py — Module 9 step 2, combine per-dataset scores into
one report with a pass/fail column.

Pure pandas — same reasoning as aggregate_tf_activity.py for not needing
the heavier LINGER_PYTHON env here.

PASS/FAIL LOGIC, stated explicitly rather than buried in a threshold:
  - Positive checks (atoh1_gfi1_pou4f3_overexpression, tbx2_conversion):
    pass if spearman_rho > 0 (predicted and real shift agree in
    direction, after this check's sign adjustment) AND auroc > 0.5.
  - Negative control (negative_control_must_not_reprogram): pass if
    the SAME statistics do NOT show a reprogramming-consistent signal —
    spearman_rho <= 0 OR auroc <= 0.5. This is the row's invert_pass=True.
  - Aging checks: no pass/fail — these validate the baseline reference
    itself (see validate_held_out.py's docstring), reported for
    completeness, not scored.

THRESHOLDS ARE A DEFAULT, NOT A FINAL ANSWER: rho>0 / auroc>0.5 is the
weakest meaningful bar (any positive signal at all). Worth tightening
once real numbers exist and can be judged against how confident the
result needs to be to trust Module 8's extrapolation into the aging
context, per the HTML plan doc's own framing of what this validation is
for.
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


def _pass_fail(row):
    if row["mode"] == "aging_tf_activity":
        return "not_scored"
    signal = (row["spearman_rho"] > 0) and (pd.notna(row["auroc"]) and row["auroc"] > 0.5)
    if row["invert_pass"]:
        return "PASS" if not signal else "FAIL"
    return "PASS" if signal else "FAIL"


report["pass_fail"] = report.apply(_pass_fail, axis=1)
report.to_csv(args.output_tsv, sep="\t", index=False)

print(f"Wrote {args.output_tsv}:")
print(report[["sample_id", "check", "spearman_rho", "auroc", "aupr", "pass_fail"]].to_string(index=False))
