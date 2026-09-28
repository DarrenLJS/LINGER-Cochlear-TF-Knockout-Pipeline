#!/usr/bin/env python3
"""
check_tf_eligibility.py — is a candidate gene a row LINGER can actually
knock out, before it goes anywhere near perturbation.knockouts?

WHY THIS EXISTS: linger_perturbation.py's knockout mechanism is "zero this
TF's row in Exp before the forward pass" — it only works for genes that are
rows in Exp, i.e. genes that are BOTH (a) present in the population
pseudobulk (Target.index) AND (b) recognized as a TF by LINGER's own
reference TF list (TFName, HOMER-motif-derived — see
prep_pseudobulk_target.py's docstring: `TFlist = list(set(Target.index) &
set(TFName[0].values))`). prep_pseudobulk_target.py already computes this
exact intersection and materializes it to {module8_dir}/Exp.tsv — so rather
than recompute the intersection ourselves (a second, potentially-drifting
implementation of the same logic), this script reads Exp.tsv directly and
checks candidate gene membership in its index. That is the real,
already-verified answer: no guessing, no re-deriving.

USE CASE (2026-09-27 pipeline extension): the Cochlea-model-testing-set-2.xlsx
generalization panel proposes knocking out Pknox2, Casz1, Zbtb20, Prdm16,
and possibly Notch1/Ctnnb1 (via GSE86204's combinatorial arm), none of which
are in the current perturbation.knockouts (Atoh1/Gfi1/Pou4f3/Tbx2). Several
biologically plausible genes in that same panel (Kcnj2, Trim71, Neu4, Dusp1,
Myo7a, Fgfr3, Ptch1) are NOT DNA-binding TFs at all (a channel, an RBP/E3
ligase, an enzyme, a phosphatase, a motor protein, two receptors) and were
already excluded from scope on that basis without needing this check.
Notch1 and Ctnnb1 are the genuinely ambiguous case: they are signaling
proteins, not TFs themselves, but their downstream effectors (Rbpj for
Notch, Tcf7/Lef1 for Wnt/beta-catenin) ARE canonical DNA-binding TFs and
may be what HOMER's motif database actually carries. This script checks
all of them rather than assuming.

Run on Eddie, in the LINGER conda env is NOT required — this only reads a
plain TSV, so any Python with pandas works. Needs Module 8's Exp.tsv to
already exist, which needs Modules 4/6 (LINGER init + GRN training) to have
run at least once first.

Usage:
    python check_tf_eligibility.py --module8-dir /path/to/module8_output_dir
    # writes tf_eligibility_report.tsv and prints a per-gene verdict

    # to check a different/extra gene list:
    python check_tf_eligibility.py --module8-dir /path/to/module8 \\
        --genes Pknox2 Casz1 Zbtb20 Prdm16 Notch1 Ctnnb1 Rbpj Tcf7 Lef1
"""
import argparse
import os

import pandas as pd

# Default candidate list: the 4 new genes the generalization panel wants to
# add (Pknox2/Casz1/Zbtb20/Prdm16, all backed by a confirmed-eligible dataset
# in the xlsx panel: GSE171921/GSE300215+GSE279618/GSE196199/GSE193046), the
# 2 ambiguous signaling genes from GSE86204 (Notch1, Ctnnb1), their
# canonical DNA-binding effectors as a fallback check (Rbpj, Tcf7, Lef1),
# and the 4 genes already in production as a sanity check that this script's
# method agrees with what the pipeline already knows works.
DEFAULT_GENES = [
    # already in perturbation.knockouts — sanity check, expect ELIGIBLE
    "Atoh1", "Gfi1", "Pou4f3", "Tbx2",
    # generalization panel, confirmed dataset available — expect to check
    "Pknox2", "Casz1", "Zbtb20", "Prdm16",
    # generalization panel, pending (GSE86204) — the actual open question
    "Notch1", "Ctnnb1",
    # canonical effectors of the two genes above, in case HOMER carries
    # these motifs instead of the receptor/co-activator name itself
    "Rbpj", "Tcf7", "Lef1",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module8-dir", required=True,
                     help="Module 8 output dir containing Exp.tsv "
                          "(written by prep_pseudobulk_target.py)")
    ap.add_argument("--genes", nargs="*", default=DEFAULT_GENES,
                     help="Candidate gene symbols to check (default: see script header)")
    ap.add_argument("--out", default="tf_eligibility_report.tsv")
    args = ap.parse_args()

    exp_path = os.path.join(args.module8_dir, "Exp.tsv")
    if not os.path.exists(exp_path):
        raise SystemExit(
            f"ERROR: {exp_path} not found. This file is written by "
            f"prep_pseudobulk_target.py (Module 8 step 1) — run that first, "
            f"or point --module8-dir at the directory it wrote into."
        )

    Exp = pd.read_csv(exp_path, sep="\t", index_col=0)
    real_tf_index = list(Exp.index)
    # Exact match set, plus a case-insensitive lookup for reporting near-misses
    # (e.g. if HOMER's reference used all-caps or human-style casing somewhere
    # rather than mouse's Title-case convention) — a near-miss is NOT treated
    # as eligible, only surfaced so a human can check whether it's a real
    # casing mismatch worth fixing or a genuinely different gene.
    exact_set = set(real_tf_index)
    lower_map = {}
    for tf in real_tf_index:
        lower_map.setdefault(tf.lower(), []).append(tf)

    print(f"Loaded {exp_path}: {len(real_tf_index)} TFs in Exp.tsv's index "
          f"(the real, already-materialized intersection of this pipeline's "
          f"population pseudobulk genes and LINGER's own reference TF list).\n")

    rows = []
    for gene in args.genes:
        exact_hit = gene in exact_set
        ci_hits = lower_map.get(gene.lower(), [])
        ci_only = (not exact_hit) and bool(ci_hits)
        if exact_hit:
            verdict = "ELIGIBLE"
        elif ci_only:
            verdict = "CASING_MISMATCH_CHECK_MANUALLY"
        else:
            verdict = "NOT_FOUND_NOT_ELIGIBLE"
        rows.append({
            "gene": gene, "verdict": verdict,
            "exact_match_in_Exp_index": exact_hit,
            "case_insensitive_matches": "; ".join(ci_hits) if ci_hits else "",
        })
        print(f"[{verdict:>28}] {gene}"
              + (f"  (found as: {', '.join(ci_hits)})" if ci_only else ""))

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)

    n_eligible = sum(r["verdict"] == "ELIGIBLE" for r in rows)
    n_not = sum(r["verdict"] == "NOT_FOUND_NOT_ELIGIBLE" for r in rows)
    n_check = sum(r["verdict"] == "CASING_MISMATCH_CHECK_MANUALLY" for r in rows)
    print(f"\n{n_eligible}/{len(rows)} ELIGIBLE (real row in Exp.tsv — safe to add to "
          f"perturbation.knockouts as a knockout, or use as an OE target once Module 8 "
          f"gets an OE mode).")
    print(f"{n_not}/{len(rows)} NOT_FOUND_NOT_ELIGIBLE (not a TF row in this pipeline's "
          f"trained network — Module 8's forward pass has nothing to zero/boost for these; "
          f"route their validation data through an expression-only check instead, never "
          f"perturbation.knockouts).")
    print(f"{n_check}/{len(rows)} CASING_MISMATCH_CHECK_MANUALLY (found under a different "
          f"case — open {args.out} and confirm by eye before treating as eligible).")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
