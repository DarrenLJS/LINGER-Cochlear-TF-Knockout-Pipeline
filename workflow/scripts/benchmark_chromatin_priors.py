"""
benchmark_chromatin_priors.py — Module 10, Option A benchmarking.

KNOWN GAP, flagged rather than guessed: I don't have visibility into the
actual downloaded file formats for GSE305205 (Hi-C loops) or the three
CUT&RUN datasets. config_eddie.yaml's comment says these are "pre-processed
peaks/loops (no realignment needed)", which most commonly means:
  - Hi-C loops : a .bedpe-style file (two anchor regions + score per row),
                 or two separate BED files of anchors.
  - CUT&RUN    : a standard narrowPeak/BED file of called peaks.
This script's loaders auto-detect the common shapes (see _load_hic_loops /
_load_cutrun_peaks below) but WERE NOT VALIDATED against your actual
downloaded files. Before trusting the benchmark numbers: `ls` the
GSE305205_hic and GSE150386_cutrun_atoh1 directories on Eddie and confirm
the file this script picks is the right one (printed to the log either way).

TF identity for CUT&RUN ground truth is inferred from each entry's
`prior_id` string (looks for a known TF name substring: atoh1, pou4f3,
gfi1). GSE181307_cutrun has no TF in its prior_id/dataset name in the
current config — CONFIRM which TF it targets before trusting its numbers;
until then it's benchmarked against every TF found in the TF-RE binding
matrix with a loud warning, which is not a rigorous comparison.

Method: for each cell type x edge-type x ground-truth-dataset combination,
label each inferred edge as a "hit" (1) if its RE overlaps a ground-truth
interval (bedtools intersect, via pybedtools) else "miss" (0), then compute
ROC-AUC and average precision (AUPR) against the inferred score — same
sklearn functions LingerGRN.Benchmk.bm_trans() itself uses, just fed a
chromatin-interval ground truth instead of a ranked-gene-list one.
"""
import argparse
import glob
import gzip
import json
import os
import re

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

TF_NAME_PATTERNS = ["atoh1", "pou4f3", "gfi1", "tbx2"]


def _parse_re_region(re_string):
    """RE strings from LingerGRN's scNN cis/TF-RE outputs are region
    identifiers like 'chr1:100000-100500' (confirmed format: N_overlap /
    O_overlap in LL_net.load_region-adjacent code, and the peaks that feed
    consensus_peaks.bed upstream). Returns (chrom, start, end) or None."""
    m = re.match(r"^(chr[\w]+):(\d+)-(\d+)$", str(re_string))
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def _is_gzipped(path):
    """Sniff the gzip magic number (0x1f 0x8b) rather than trusting the
    file extension — some GEO downloads are gzip without a .gz suffix, and
    occasionally the reverse (a .gz-named file that's already plain text)."""
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def _load_bed_like(path):
    """Best-effort loader for a BED/narrowPeak/bedpe-ish file: assumes the
    first 3 whitespace/tab-delimited columns of any non-header row are
    chrom/start/end. For bedpe (Hi-C loops), ALSO extracts columns 4-6 as a
    second anchor if present, and treats both anchors as valid overlap
    targets (a cis edge landing in either loop anchor counts as supported).
    Transparently handles gzip-compressed input (most GEO narrowPeak/tsv
    downloads are shipped as .gz)."""
    rows = []
    opener = gzip.open if _is_gzipped(path) else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                rows.append((parts[0], int(parts[1]), int(parts[2])))
                if len(parts) >= 6:
                    rows.append((parts[3], int(parts[4]), int(parts[5])))
            except ValueError:
                continue  # header row that slipped through
    return rows


def _find_ground_truth_files(directory, role=None, file_glob=None):
    """Returns a LIST of matched files (not just one), so multiple
    replicates / multiple named contexts (e.g. two file_glob-selected
    POU4F3 CUT&RUN conditions) can be unioned by the caller.

    If file_glob is given, it takes priority and is matched literally
    against files directly under `directory` (not recursively — these
    directories are flat GEO dumps, and a literal glob keeps selection
    explicit and auditable from config rather than "whatever sorted-first
    matched a generic pattern").

    Otherwise falls back to auto-detection, using role-specific patterns:
    Hi-C loop files are commonly .tsv/.bedpe/loop-named, not "peak"-named,
    so a single peak-oriented pattern list (the old behaviour) silently
    missed real Hi-C files that don't happen to have "bed"/"peak" in the
    name — as it did for GSE305205's *_10kb_01.tsv.gz."""
    if file_glob:
        hits = sorted(glob.glob(os.path.join(directory, file_glob)))
        return [h for h in hits if os.path.isfile(h)]

    if role == "hic_loops":
        patterns = ["*.bedpe*", "*loop*", "*.tsv*", "*.bed*"]
    else:  # cutrun and anything else peak-shaped
        patterns = ["*.narrowPeak*", "*peak*", "*.bed*", "*.bedpe*"]

    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(directory, "**", pat), recursive=True))
        hits = [h for h in hits if os.path.isfile(h)]
        if hits:
            return hits[:1]  # keep prior single-file-auto-detect behaviour
    return []


def _overlaps_any(chrom, start, end, intervals_by_chrom):
    for ic, istart, iend in intervals_by_chrom.get(chrom, []):
        if start < iend and end > istart:
            return True
    return False


def _tf_from_prior_id(prior_id, dataset_desc=""):
    text = (prior_id + " " + dataset_desc).lower()
    for tf in TF_NAME_PATTERNS:
        if tf in text:
            return tf.capitalize() if tf != "pou4f3" else "Pou4f3"
    return None


def benchmark_cis(module4_dir, celltypes, hic_intervals, outdir):
    rows = []
    intervals_by_chrom = {}
    for c, s, e in hic_intervals:
        intervals_by_chrom.setdefault(c, []).append((c, s, e))

    for ct in celltypes:
        path = os.path.join(module4_dir, f"cell_type_specific_cis_regulatory_{ct}.txt")
        if not os.path.exists(path):
            print(f"[cis:{ct}] missing {path}, skipping")
            continue
        df = pd.read_csv(path, sep="\t", header=None, names=["RE", "TG", "score"])
        labels, scores = [], []
        for _, row in df.iterrows():
            parsed = _parse_re_region(row["RE"])
            if parsed is None:
                continue
            chrom, start, end = parsed
            labels.append(int(_overlaps_any(chrom, start, end, intervals_by_chrom)))
            scores.append(row["score"])
        if len(set(labels)) < 2:
            print(f"[cis:{ct}] only one class present ({sum(labels)}/{len(labels)} hits) — AUC undefined, skipping")
            continue
        auc = roc_auc_score(labels, scores)
        aupr = average_precision_score(labels, scores)
        rows.append({"celltype": ct, "edge_type": "cis_RE_TG", "ground_truth": "hic_loops",
                      "n_edges": len(labels), "n_hits": sum(labels), "auroc": auc, "aupr": aupr})
        print(f"[cis:{ct}] n={len(labels)} hits={sum(labels)} AUROC={auc:.3f} AUPR={aupr:.3f}")
    return rows


def benchmark_tf_re(module4_dir, celltypes, cutrun_datasets, outdir):
    rows = []
    for prior_id, intervals, tf_name in cutrun_datasets:
        intervals_by_chrom = {}
        for c, s, e in intervals:
            intervals_by_chrom.setdefault(c, []).append((c, s, e))

        for ct in celltypes:
            path = os.path.join(module4_dir, f"cell_type_specific_TF_RE_binding_{ct}.txt")
            if not os.path.exists(path):
                print(f"[tfre:{ct}:{prior_id}] missing {path}, skipping")
                continue
            df = pd.read_csv(path, sep="\t", index_col=0)
            tf_cols = [tf_name] if (tf_name and tf_name in df.columns) else list(df.columns)
            if tf_name and tf_name not in df.columns:
                print(f"[tfre:{ct}:{prior_id}] WARNING: TF '{tf_name}' not in binding matrix columns, "
                      f"falling back to all {len(df.columns)} TFs — treat this row cautiously")
            for tf_col in tf_cols:
                labels, scores = [], []
                for re_string, score in df[tf_col].items():
                    parsed = _parse_re_region(re_string)
                    if parsed is None:
                        continue
                    chrom, start, end = parsed
                    labels.append(int(_overlaps_any(chrom, start, end, intervals_by_chrom)))
                    scores.append(score)
                if len(set(labels)) < 2:
                    continue
                auc = roc_auc_score(labels, scores)
                aupr = average_precision_score(labels, scores)
                rows.append({"celltype": ct, "edge_type": f"TF_RE_binding[{tf_col}]",
                             "ground_truth": prior_id, "n_edges": len(labels),
                             "n_hits": sum(labels), "auroc": auc, "aupr": aupr})
                print(f"[tfre:{ct}:{prior_id}:{tf_col}] n={len(labels)} hits={sum(labels)} "
                      f"AUROC={auc:.3f} AUPR={aupr:.3f}")
    return rows


p = argparse.ArgumentParser()
p.add_argument("--module4-dir", required=True)
p.add_argument("--celltypes", required=True, nargs="+")
p.add_argument("--chrom-priors-json", required=True)
p.add_argument("--outdir", required=True)
p.add_argument("--output-report", required=True)
args = p.parse_args()

os.makedirs(args.outdir, exist_ok=True)
priors = json.loads(args.chrom_priors_json)

hic_intervals = []
cutrun_datasets = []
for entry in priors:
    role = entry.get("role")
    if role == "skip":
        reason = entry.get("skip_reason", "no reason given")
        print(f"[{entry['prior_id']}] SKIPPED by config (role: skip) — {reason}")
        continue
    directory = entry["path"]
    gt_files = _find_ground_truth_files(directory, role=role, file_glob=entry.get("file_glob"))
    if not gt_files:
        print(f"[{entry['prior_id']}] WARNING: no ground-truth file found under {directory}"
              f"{' matching file_glob=' + entry['file_glob'] if entry.get('file_glob') else ''}, skipping")
        continue
    print(f"[{entry['prior_id']}] using {gt_files} (role={role})")
    intervals = []
    for gt_file in gt_files:
        intervals.extend(_load_bed_like(gt_file))
    if role == "hic_loops":
        hic_intervals.extend(intervals)
    elif role == "cutrun":
        tf_name = _tf_from_prior_id(entry["prior_id"], entry.get("dataset", ""))
        if tf_name is None:
            print(f"[{entry['prior_id']}] WARNING: could not infer TF identity from prior_id — "
                  f"CONFIRM manually and hardcode it here before trusting results")
        cutrun_datasets.append((entry["prior_id"], intervals, tf_name))
    # atac_prior / subtype_rna / subtype_atac / h3k4me3 roles: not used as
    # cis/TF-RE ground truth here — out of scope for Option A as decided.

all_rows = []
if hic_intervals:
    all_rows += benchmark_cis(args.module4_dir, args.celltypes, hic_intervals, args.outdir)
else:
    print("WARNING: no Hi-C loop intervals loaded — cis benchmarking skipped entirely")

if cutrun_datasets:
    all_rows += benchmark_tf_re(args.module4_dir, args.celltypes, cutrun_datasets, args.outdir)
else:
    print("WARNING: no CUT&RUN datasets loaded — TF-RE benchmarking skipped entirely")

report = pd.DataFrame(all_rows)
report.to_csv(args.output_report, sep="\t", index=False)
print(f"\nWrote {args.output_report}: {len(report)} benchmark rows")
