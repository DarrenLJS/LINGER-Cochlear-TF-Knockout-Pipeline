"""
validate_held_out.py — Module 9, one held-out/negative-control dataset's
validation score.

Methodology per the HTML plan doc's Module 9 spec, taken at face value
rather than invented from scratch:
  - Output: "Correlation metrics, AUROC/AUPR" (stated explicitly).
  - Positive checks compare a Module 8 KNOCKOUT prediction directly
    against real overexpression/conversion data ("Atoh1/Gfi1/Pou4f3
    prediction vs. real overexpression data") — NOT a separately-built
    overexpression simulation. This script therefore reuses whichever
    Module 8 `*_predicted_expression.tsv` the check config names, applying
    a sign flip for checks where the held-out ground truth is the
    OPPOSITE direction of what Module 8 simulated (TF gain vs TF loss).
  - Aging checks use Module 7's role="baseline" regulon scores as the
    "unperturbed" reference point, per the plan doc and last session's
    revised-plan table, diffed against each aging dataset's own
    TF-activity profile (computed here the same way Module 7 does, via
    TF_activity.regulon(network="cell population")).
    REVISED 2026-09-13 per linger_cochlear_pipeline_progress02.html: the
    plan doc names GSE274279 the "ground truth for aging vector" and the
    other four aging datasets "aging-vector-support" — i.e. the four
    support datasets are meant to be checked AGAINST GSE274279's shift
    vector, not merely diffed against baseline in isolation. Two roles,
    set via check_spec["aging_role"]:
      - "reference" (GSE274279 / check=aging_vector): computes its own
        real_shift vs Module 7 baseline and writes it to
        --output-shift-tsv as "the aging vector". No external ground
        truth exists for this row, so its own rho/auroc/aupr are NaN —
        it's the reference, not a result to score.
      - "support" (check=aging_vector_support): computes its own
        real_shift the same way, then Spearman-correlates it against the
        reference vector (read from --aging-vector-tsv) and computes
        AUROC/AUPR treating the reference's top-K |shift| TFs as the
        "hit" ground truth and this dataset's |shift| as the score. This
        is the actual "does this independent aging dataset agree with
        the aging vector's direction" check the plan doc describes.
  - Negative control: same shift-comparison machinery as the positive
    reprogramming checks, but the PASS condition is inverted — this
    dataset's real profile should NOT correlate with the reprogramming
    direction the way GSE224627/GSE233559's real profiles do.

TWO METRICS COMPUTED FOR EVERY DATASET, matching the HTML's stated output
exactly (not picking one over the other):
  1. Spearman correlation between the predicted shift vector (from Module
     8's knockout, sign-adjusted per check) and the real shift vector
     (held-out expression minus Module 4's population pseudobulk
     baseline), over the shared gene set.
  2. AUROC/AUPR, reusing Module 10's existing metric choice rather than
     inventing a new one: the top-K real-shifted genes (|real_shift|,
     K configurable) are treated as the binary ground-truth "hit" set,
     ranked against |predicted_shift| as the score.

FLAGGED, NOT GUESSED: exactly which cells/samples in each held-out
dataset represent the "converted"/"reprogrammed" state vs a baseline
state isn't specified anywhere upstream (these are whole downloaded
datasets, not pre-split into before/after). This script uses the
dataset's OVERALL mean expression profile against Module 4's population
pseudobulk baseline as the real shift — a reasonable default, but WORTH
CONFIRMING against each dataset's actual structure (e.g. if GSE224627
has explicit condition labels, a proper per-condition diff would be more
rigorous than an overall-mean shift). Uses bulk_rna_loader's same
best-effort format auto-detector as Module 7, with the same caveat about
unverified per-dataset layouts.

ADDED 2026-09-19 — --baseline-tsv / --network, for Module 9b (cell-type-
resolved validation, 09b_validation_celltype.smk). Both default to the
exact population behavior above (Module 9 itself calls this script with
neither flag set, so it is completely unaffected), so this file has ONE
implementation shared by both modules rather than a forked duplicate:
  - --baseline-tsv overrides the "before" reference used by
    mode=="expression_shift" (default: Module 4's population
    TG_pseudobulk.tsv). Module 9b points this at a cell type's own
    TG_pseudobulk_{celltype}.tsv (already written by Module 8b's
    prep_pseudobulk_target_celltype) instead — diffing a cell type's
    knockout prediction against the POPULATION baseline would conflate
    genuine cross-cell-type expression differences with the knockout's
    actual effect; diffing against that cell type's own pre-knockout
    profile isolates the knockout effect the same way Module 9 already
    does at the population level.
  - --network overrides the LingerGRN `network` argument to
    TF_activity.regulon() for mode=="aging_tf_activity" (default: "cell
    population"). Passing a cell-type name here instead is a real,
    supported LingerGRN 1.110 code path — confirmed by reading
    TF_activity.py directly: any value other than "cell population"/
    "general" is read as `cell_type_specific_trans_regulatory_{network}.
    txt`, which Module 6 (grn_celltype_specific.py) already writes, one
    per cell type, into the same --workdir this script already reads.

CAVEAT that applies to EVERY Module 9b check, not just the aging ones,
worth restating here since it's easy to lose sight of once the plumbing
works: every held-out/negative-control/aging dataset loaded by this
script (via bulk_rna_loader.load_rna_only) is BULK RNA-seq, not single-
cell/multiome — there is no per-cell-type label to subset by on the real
side of any comparison. Pointing --network or --baseline-tsv at a
specific cell type changes which MODEL (regulatory network / pre-
knockout reference) the same bulk sample is scored against, not which
CELLS from that sample are used — it answers "does this bulk sample look
consistent with cell type X's regulatory activity / knockout response",
not "what did cell type X's cells in this sample actually do" (no such
ground truth exists in this dataset collection). Report results from
this mode with that distinction stated, not implied.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from bulk_rna_loader import load_rna_only  # noqa: E402

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True, help="Module 4 workdir — data/TG_pseudobulk.tsv baseline lives here")
p.add_argument("--grn-dir", required=True)
p.add_argument("--genome", required=True)
p.add_argument("--module8-dir", required=True)
p.add_argument("--module7-summary", required=True, help="module7_tf_activity/tf_activity_summary.tsv")
p.add_argument("--entry-json", required=True, help="held-out/negative-control config entry")
p.add_argument("--check-spec-json", required=True,
                help='{"mode": "expression_shift"|"aging_tf_activity", '
                     '"ko_id": "<module8 ko_id>"|null, "sign": 1|-1, '
                     '"invert_pass": true|false, "top_k": 200}')
p.add_argument("--output-tsv", required=True)
p.add_argument("--output-shift-tsv", default=None,
               help="Only used when check_spec['aging_role']=='reference': writes this "
                    "dataset's real TF-activity shift vector so aging_vector_support "
                    "checks can be correlated against it.")
p.add_argument("--aging-vector-tsv", default=None,
               help="Only used when check_spec['aging_role']=='support': path to the "
                    "reference ('aging_vector') shift vector written via --output-shift-tsv.")
p.add_argument("--baseline-tsv", default=None,
               help="ADDED for Module 9b. Overrides the 'before' reference for "
                    "mode=='expression_shift' (default: WORKDIR/data/TG_pseudobulk.tsv, "
                    "population). Pass a cell type's own TG_pseudobulk_{celltype}.tsv for "
                    "cell-type-resolved validation.")
p.add_argument("--network", default="cell population",
               help="ADDED for Module 9b. Overrides the LingerGRN `network` arg to "
                    "TF_activity.regulon() for mode=='aging_tf_activity' (default: 'cell "
                    "population'). Pass a cell type name to use Module 6's "
                    "cell_type_specific_trans_regulatory_{network}.txt instead — a real "
                    "LingerGRN 1.110 code path, not a population-only special case.")
args = p.parse_args()

entry = json.loads(args.entry_json)
spec = json.loads(args.check_spec_json)
sample_id = entry.get("sample_id")
check = entry.get("check")
top_k = spec.get("top_k", 200)

print(f"--- {sample_id} (check={check}, mode={spec['mode']}) ---")


def _auroc_aupr(real_shift, predicted_shift, top_k):
    common = real_shift.index.intersection(predicted_shift.index)
    real = real_shift.loc[common]
    pred = predicted_shift.loc[common]
    finite = real.notna() & pred.notna()
    n_dropped = int((~finite).sum())
    if n_dropped:
        print(f"_auroc_aupr: dropping {n_dropped}/{len(finite)} entries with NaN "
              f"(undefined TF-activity/shift, e.g. zero-variance regulon in a small "
              f"dataset) before scoring: {list(real.index[~finite])}")
    real, pred = real[finite], pred[finite]
    if len(real) < 10:
        return np.nan, np.nan, len(real)
    k = min(top_k, len(real) // 2)
    if k < 2:
        return np.nan, np.nan, len(real)
    threshold = real.abs().sort_values(ascending=False).iloc[k - 1]
    labels = (real.abs() >= threshold).astype(int).values
    scores = pred.abs().values
    if len(set(labels)) < 2:
        return np.nan, np.nan, len(real)
    return roc_auc_score(labels, scores), average_precision_score(labels, scores), len(real)


if spec["mode"] == "expression_shift":
    # --- real shift: held-out dataset's overall mean expression vs Module 4's population baseline ---
    adata = load_rna_only(entry)
    real_mean = pd.Series(
        np.asarray(adata.X.mean(axis=0)).reshape(-1), index=adata.var_names
    )
    real_mean = np.log2(1 + real_mean.clip(lower=0))

    baseline_path = args.baseline_tsv or os.path.join(args.workdir, "data", "TG_pseudobulk.tsv")
    baseline = pd.read_csv(baseline_path, sep=",", header=0, index_col=0)
    baseline_mean = np.log2(1 + baseline.mean(axis=1).clip(lower=0))
    print(f"baseline: {baseline_path}")

    common_genes = real_mean.index.intersection(baseline_mean.index)
    real_shift = (real_mean.loc[common_genes] - baseline_mean.loc[common_genes])
    print(f"real_shift: {len(real_shift)} genes (held-out {sample_id} vs population baseline)")

    ko_id = spec["ko_id"]
    ko_path = os.path.join(args.module8_dir, f"{ko_id}_predicted_expression.tsv")
    ko_pred = pd.read_csv(ko_path, sep="\t", index_col=0)
    ko_pred_mean = ko_pred.mean(axis=1)
    predicted_shift = (ko_pred_mean - baseline_mean.reindex(ko_pred_mean.index)) * spec["sign"]
    print(f"predicted_shift: {len(predicted_shift)} genes (Module 8 {ko_id}, sign={spec['sign']})")

    common = real_shift.index.intersection(predicted_shift.index)
    real_c, pred_c = real_shift.loc[common], predicted_shift.loc[common]
    finite = real_c.notna() & pred_c.notna()
    if (~finite).sum():
        print(f"spearman: dropping {int((~finite).sum())}/{len(finite)} NaN entries before correlating")
    real_c, pred_c = real_c[finite], pred_c[finite]
    if len(real_c) < 10:
        rho, pval = np.nan, np.nan
    else:
        rho, pval = spearmanr(real_c, pred_c)
    auroc, aupr, n_common = _auroc_aupr(real_shift, predicted_shift, top_k)

elif spec["mode"] == "aging_tf_activity":
    # --- real: held-out dataset's own TF activity, via the same LingerGRN call Module 7 uses ---
    import LingerGRN.TF_activity as TF_activity
    adata = load_rna_only(entry)
    os.chdir(args.workdir)
    real_tf_activity = TF_activity.regulon(args.workdir + "/", adata, args.grn_dir, args.network, args.genome)
    print(f"regulon network: {args.network!r}")
    real_mean = real_tf_activity.mean(axis=1)

    # --- baseline: Module 7's role="baseline" mean TF activity ---
    m7 = pd.read_csv(args.module7_summary, sep="\t", index_col=0)
    baseline_rows = m7[m7["role"] == "baseline"].drop(columns=["role"])
    baseline_mean = baseline_rows.mean(axis=0)  # mean over the 6 baseline datasets, per TF
    baseline_mean.index.name = None

    common_tfs = real_mean.index.intersection(baseline_mean.index)
    real_shift = real_mean.loc[common_tfs] - baseline_mean.loc[common_tfs]

    aging_role = spec.get("aging_role")
    if aging_role == "reference":
        # This dataset IS "the aging vector" (GSE274279 per the plan doc) —
        # there's no external ground truth to score it against, so it's
        # reported descriptively (NaN rho/auroc/aupr) and written out for
        # the "support" datasets to correlate against.
        if args.output_shift_tsv:
            real_shift.rename("shift").to_frame().to_csv(args.output_shift_tsv, sep="\t")
            print(f"aging vector: wrote {len(real_shift)}-TF reference shift to {args.output_shift_tsv}")
        else:
            print("WARNING: aging_role=='reference' but --output-shift-tsv not given — "
                  "aging_vector_support checks will have nothing to correlate against")
        rho, pval = np.nan, np.nan
        auroc, aupr, n_common = np.nan, np.nan, len(real_shift)

    elif aging_role == "support":
        if not args.aging_vector_tsv:
            raise ValueError("check_spec['aging_role']=='support' requires --aging-vector-tsv "
                              "(the reference aging_vector dataset's shift file)")
        aging_vector = pd.read_csv(args.aging_vector_tsv, sep="\t", index_col=0)["shift"]
        common = real_shift.index.intersection(aging_vector.index)
        print(f"aging TF-activity shift computed for {len(real_shift)} TFs "
              f"({len(common)} shared with the reference aging vector)")
        av_c, rs_c = aging_vector.loc[common], real_shift.loc[common]
        finite = av_c.notna() & rs_c.notna()
        if (~finite).sum():
            print(f"spearman: dropping {int((~finite).sum())}/{len(finite)} NaN entries "
                  f"before correlating (undefined TF-activity, likely zero-variance "
                  f"regulon genes in this dataset): {list(av_c.index[~finite])}")
        av_c, rs_c = av_c[finite], rs_c[finite]
        if len(av_c) < 10:
            rho, pval = np.nan, np.nan
        else:
            rho, pval = spearmanr(av_c, rs_c)
        # Reference vector's top-K |shift| TFs = ground-truth "hit" set;
        # this dataset's |shift| = score — same convention as _auroc_aupr
        # uses elsewhere (real_shift=ground truth, predicted_shift=score).
        auroc, aupr, n_common = _auroc_aupr(aging_vector, real_shift, top_k)

    else:
        raise ValueError("check_spec['mode']=='aging_tf_activity' requires "
                          "check_spec['aging_role'] to be 'reference' or 'support', "
                          f"got: {aging_role!r}")

else:
    raise ValueError(f"Unknown check spec mode: {spec['mode']}")

result = {
    "sample_id": sample_id,
    "check": check,
    "mode": spec["mode"],
    "ko_id": spec.get("ko_id"),
    "sign": spec.get("sign"),
    "invert_pass": spec.get("invert_pass", False),
    "n_genes_or_tfs": n_common,
    "spearman_rho": rho,
    "spearman_pval": pval,
    "auroc": auroc,
    "aupr": aupr,
}
pd.DataFrame([result]).to_csv(args.output_tsv, sep="\t", index=False)
print(f"Wrote {args.output_tsv}: rho={rho}, auroc={auroc}, aupr={aupr}")
