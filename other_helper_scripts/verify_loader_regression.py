"""
One-off regression check: rerun the 20 datasets that already succeeded
against the ORIGINAL loader, now through the REWRITTEN bulk_rna_loader.py,
and diff numerically against the real existing output. Doesn't touch or
overwrite the originals.

Run on Eddie from the pipeline's code directory:
    cd /exports/eddie/scratch/s2906787/linger_pipeline/code
    /exports/csce/eddie/biology/groups/bioinfmsc/anaconda/envs/s2906787/LINGER/bin/python \
        verify_loader_regression.py
"""
import json
import os
import subprocess
import sys

import yaml
import pandas as pd
import numpy as np

CODE_DIR = "/exports/eddie/scratch/s2906787/linger_pipeline/code"
LINGER_PYTHON = "/exports/csce/eddie/biology/groups/bioinfmsc/anaconda/envs/s2906787/LINGER/bin/python"
WORKDIR = "/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/module4_linger_init"
GRN_DIR = "/exports/eddie/scratch/s2906787/linger_pipeline/provide_data/provide_data"
GENOME = "mm10"
ORIG_DIR = "/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/module7_tf_activity"
VERIFY_DIR = "/exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/module7_verify_20"
os.makedirs(VERIFY_DIR, exist_ok=True)

KNOWN_GOOD = [
    "GSE116702", "GSE135703_adult_SC", "GSE136196_stria", "GSE157398_Base_RNA",
    "GSE163798", "GSE168041_sgn_lateralwall", "GSE193158", "GSE198167",
    "GSE202920_P14_P28", "GSE213784", "GSE240605", "GSE242669", "GSE275243",
    "GSE283534_mature_hc", "GSE294736", "GSE299064", "GSE300428", "GSE312253",
    "GSE329565", "GSE83599",
]

with open(os.path.join(CODE_DIR, "config/config_eddie.yaml")) as f:
    config = yaml.safe_load(f)

path_roots = {
    "multiome_root": config["multiome_root"],
    "chrom_root": config["chrom_root"],
    "scratch": config["scratch"],
}


def resolve(obj):
    if isinstance(obj, str):
        return obj.format(**path_roots)
    if isinstance(obj, dict):
        return {k: resolve(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [resolve(v) for v in obj]
    return obj


bulk = [dict(resolve(e), role="bulk") for e in config.get("extra_bulk_rna_inputs", [])]
baseline = [dict(resolve(e), role="baseline") for e in config.get("model_construction_refs", [])]
entries = {e.get("sample_id", e.get("ref_id")): e for e in (bulk + baseline)}

results = []
for entry_id in KNOWN_GOOD:
    if entry_id not in entries:
        print(f"SKIP {entry_id}: not found in config (check spelling / config drift)")
        results.append((entry_id, "NOT_IN_CONFIG", None))
        continue
    entry = entries[entry_id]
    out_path = os.path.join(VERIFY_DIR, f"{entry_id}.regulon.tsv")
    cmd = [
        LINGER_PYTHON, "workflow/scripts/linger_tf_activity.py",
        "--workdir", WORKDIR, "--grn-dir", GRN_DIR, "--genome", GENOME,
        "--entry-json", json.dumps(entry), "--output-tsv", out_path,
    ]
    print(f"--- {entry_id} ---")
    r = subprocess.run(cmd, cwd=CODE_DIR, capture_output=True, text=True)
    if r.returncode != 0:
        # MemoryError inside regulon()'s quantile_normalize() is very
        # likely this SCRIPT's fault, not the loader's: the real pipeline
        # rule requests vmem_mb: 64000 via SGE, but this script runs via
        # plain subprocess.run() with no SGE submission, so it inherits
        # whatever memory limit the interactive/login shell actually has —
        # almost certainly far less. Labelled separately so the summary
        # doesn't conflate this with a real loader regression.
        if "MemoryError" in r.stderr and "bulk_rna_loader.py" not in r.stderr:
            print(f"OOM (likely this script's interactive memory limit, not the loader — see comment): {entry_id}")
            print(r.stderr[-500:])
            results.append((entry_id, "OOM_LIKELY_SCRIPT_MEMORY_LIMIT", None))
        else:
            print(f"FAIL (rewritten loader crashed on a previously-successful dataset!): {entry_id}")
            print(r.stderr[-2000:])
            results.append((entry_id, "CRASHED", None))
        continue

    orig_path = os.path.join(ORIG_DIR, f"{entry_id}.regulon.tsv")
    if not os.path.exists(orig_path):
        results.append((entry_id, "NO_ORIGINAL_TO_COMPARE", None))
        continue

    orig = pd.read_csv(orig_path, sep="\t", index_col=0)
    new = pd.read_csv(out_path, sep="\t", index_col=0)
    if orig.shape != new.shape:
        results.append((entry_id, "SHAPE_MISMATCH", f"{orig.shape} vs {new.shape}"))
        continue
    if not (orig.index.equals(new.index) and orig.columns.equals(new.columns)):
        results.append((entry_id, "LABEL_MISMATCH", None))
        continue

    # NaN-aware comparison: regulon() can legitimately produce NaN (e.g. a
    # TF/RE pair with zero variance), and plain np.abs(a - b).max() returns
    # NaN if EITHER array has any NaN anywhere, regardless of whether the
    # two are otherwise identical — that was silently producing the
    # inconclusive "VALUE_DIFF, nan" results in the first run. Compare NaN
    # positions explicitly, then compare only the non-NaN cells.
    orig_nan_mask = np.isnan(orig.values)
    new_nan_mask = np.isnan(new.values)
    if not np.array_equal(orig_nan_mask, new_nan_mask):
        results.append((entry_id, "NAN_POSITION_MISMATCH",
                         f"{orig_nan_mask.sum()} NaNs (orig) vs {new_nan_mask.sum()} (new)"))
        continue
    non_nan = ~orig_nan_mask
    if non_nan.sum() == 0:
        max_abs_diff = 0.0  # everything NaN in both — vacuously identical
    else:
        max_abs_diff = np.abs(orig.values[non_nan] - new.values[non_nan]).max()
    results.append((entry_id, "OK" if max_abs_diff < 1e-8 else "VALUE_DIFF", max_abs_diff))

print("\n=== SUMMARY ===")
for entry_id, status, detail in results:
    print(f"{entry_id:35s} {status:20s} {detail if detail is not None else ''}")

n_ok = sum(1 for _, s, _ in results if s == "OK")
print(f"\n{n_ok}/{len(results)} datasets confirmed identical output.")
