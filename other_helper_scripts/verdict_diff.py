#!/usr/bin/env python3
"""
verdict_diff.py — regression guard: did a code/config change alter any validation verdict?

Compares two validation_report_combined.tsv files (the BEFORE you saved, the AFTER the rerun) row by row on
(scope, sample_id) and reports every change in `pass_fail` (and, with --rho-tol, any rho that moved by more than
the tolerance). Rows present in only one file are listed. Exit status 0 = no verdict changed, 1 = at least one did,
so it can gate a rerun.

Stages 1-2 of the 2026-10-15 build (figures, config, provenance, funnel, Module 9c) must change NO verdicts:
    cp $SCRATCH/module9b_validation_celltype/validation_report_combined.tsv  before_combined.tsv     # BEFORE rerunning
    ... rerun ...
    python3 other_helper_scripts/verdict_diff.py before_combined.tsv \
        $SCRATCH/module9b_validation_celltype/validation_report_combined.tsv

Stage 5 (screen-based null for the headline scopes) is EXPECTED to change some negative-control verdicts; run the
same command to get that as a labelled, separate diff.

Usage: verdict_diff.py BEFORE AFTER [--rho-tol 0.005] [--out diff.tsv]
"""
import argparse
import sys

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("before")
p.add_argument("after")
p.add_argument("--rho-tol", type=float, default=0.005, help="report |delta spearman_rho| above this (default 0.005)")
p.add_argument("--out", default=None, help="write the full row-level comparison here")
a = p.parse_args()

KEY = ["scope", "sample_id"]
b, f = pd.read_csv(a.before, sep="\t"), pd.read_csv(a.after, sep="\t")
for name, d in (("BEFORE", b), ("AFTER", f)):
    miss = [c for c in KEY + ["pass_fail"] if c not in d.columns]
    if miss:
        sys.exit(f"{name} lacks column(s) {miss}")
m = b.merge(f, on=KEY, how="outer", suffixes=("_before", "_after"), indicator=True)
only_b, only_a = m[m["_merge"] == "left_only"], m[m["_merge"] == "right_only"]
both = m[m["_merge"] == "both"].copy()
both["verdict_changed"] = both["pass_fail_before"].astype(str) != both["pass_fail_after"].astype(str)
if "spearman_rho_before" in both and "spearman_rho_after" in both:
    both["d_rho"] = both["spearman_rho_after"] - both["spearman_rho_before"]
    both["rho_moved"] = both["d_rho"].abs() > a.rho_tol
else:
    both["d_rho"], both["rho_moved"] = float("nan"), False
chg = both[both["verdict_changed"]]
print(f"{len(both)} rows compared; {len(chg)} verdict change(s); {int(both['rho_moved'].sum())} rho shift(s) > {a.rho_tol}; "
      f"{len(only_b)} row(s) only in BEFORE, {len(only_a)} only in AFTER")
if len(chg):
    print("\nVERDICT CHANGES:")
    print(chg[KEY + ["pass_fail_before", "pass_fail_after", "d_rho"]].to_string(index=False))
mv = both[both["rho_moved"] & ~both["verdict_changed"]]
if len(mv):
    print("\nrho moved but verdict unchanged:")
    print(mv[KEY + ["pass_fail_after", "d_rho"]].to_string(index=False))
if len(only_b) or len(only_a):
    print("\nrows only in one file:")
    print(pd.concat([only_b[KEY].assign(side="BEFORE only"), only_a[KEY].assign(side="AFTER only")]).to_string(index=False))
if a.out:
    both.to_csv(a.out, sep="\t", index=False)
sys.exit(1 if len(chg) else 0)
