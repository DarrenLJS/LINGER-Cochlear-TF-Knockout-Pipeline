"""
qc_preprocess_sample.py — Module 1 entrypoint, ONE sample per invocation.

Invoked by Snakemake's `script:` directive (see workflow/rules/01_qc_preprocessing.smk),
which injects the `snakemake` object with `.wildcards`, `.output`, `.params`,
`.log`. Loads RNA+ATAC, QC-filters, runs Scrublet, intersects barcodes, writes
rna_filtered.h5ad / atac_filtered.h5ad, and a per-sample summary JSON that
workflow/scripts/aggregate_module1_summary.py later collects into
module1_summary.csv (replicating the original single-script run's combined
summary table, now assembled from N independent job outputs instead of one
in-process loop).
"""
import json
import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_lib import LOADERS, qc_filter_rna, run_scrublet, qc_filter_atac

sample = snakemake.params.sample_entry
qc_params = snakemake.params.qc_params
sample_id = sample["sample_id"]

log_path = snakemake.log[0]
log_f = open(log_path, "w")


def logprint(*a):
    msg = " ".join(str(x) for x in a)
    print(msg)
    print(msg, file=log_f)


try:
    logprint(f"=== {sample_id} ({sample['type']}) ===")

    loader = LOADERS[sample["type"]]
    rna, atac = loader(sample)
    logprint(f"  Loaded: RNA {rna.shape}, ATAC {atac.shape}")

    rna = qc_filter_rna(rna, qc_params)
    logprint(f"  RNA after gene/mito QC: {rna.shape}")

    rna = run_scrublet(rna, qc_params, sample)
    logprint(f"  RNA after doublet removal: {rna.shape}")

    atac = qc_filter_atac(atac, qc_params)
    logprint(f"  ATAC after peak QC: {atac.shape}")

    shared_barcodes = rna.obs_names.intersection(atac.obs_names)
    n_shared = len(shared_barcodes)
    logprint(f"  Shared barcodes between RNA and ATAC: {n_shared}")
    if n_shared < 50:
        logprint(
            f"  WARNING: only {n_shared} shared barcodes — RNA and ATAC may not be "
            "true paired multiome for this sample. Flagging for manual review rather than dropping."
        )

    if n_shared > 0:
        rna_paired = rna[shared_barcodes].copy()
        atac_paired = atac[shared_barcodes].copy()
    else:
        rna_paired, atac_paired = rna, atac  # keep unpaired data for manual inspection

    rna_paired.write_h5ad(snakemake.output.rna)
    atac_paired.write_h5ad(snakemake.output.atac)

    summary = {
        "sample_id": sample_id,
        "dataset": sample["dataset"],
        "n_cells_final": rna_paired.n_obs,
        "n_genes_final": rna_paired.n_vars,
        "n_peaks_final": atac_paired.n_vars,
        "n_shared_barcodes": n_shared,
        "status": "OK",
    }
    logprint(f"  Saved rna_filtered.h5ad + atac_filtered.h5ad for {sample_id}")

except Exception as e:
    logprint(f"  FAILED: {sample_id} — {e}")
    traceback.print_exc(file=log_f)
    summary = {"sample_id": sample_id, "dataset": sample["dataset"], "status": "FAILED", "error": str(e)}
    with open(snakemake.output.summary, "w") as f:
        json.dump(summary, f, indent=2)
    log_f.close()
    raise  # let Snakemake mark this job as failed; --keep-going still lets other samples proceed

with open(snakemake.output.summary, "w") as f:
    json.dump(summary, f, indent=2)
log_f.close()
