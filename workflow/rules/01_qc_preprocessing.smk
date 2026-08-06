# =============================================================================
# workflow/rules/01_qc_preprocessing.smk
# Module 1 — QC & preprocessing
#
# One job per sample (wildcard `sample`), replacing the original script's
# internal `for sample in SAMPLES` loop with try/except-per-sample. Snakemake
# gives the same "one failure doesn't lose the rest" property natively via
# --keep-going, plus proper output-based staleness tracking instead of the
# original's manual `if os.path.exists(rna_out) and os.path.exists(atac_out)`
# skip check.
#
# Conda env: envs/linger_preproc.yaml (scanpy, anndata, scrublet)
# =============================================================================

rule qc_preprocess_sample:
    """
    Load RNA+ATAC for one sample, QC-filter, run Scrublet, intersect
    barcodes, write rna_filtered.h5ad + atac_filtered.h5ad.
    Logic in workflow/scripts/qc_lib.py + qc_preprocess_sample.py, ported
    unchanged from the original 01_qc_preprocessing.py.
    """
    input:
        raw = raw_inputs_for_sample,
    output:
        rna     = f"{SCRATCH}/{{sample}}/rna_filtered.h5ad",
        atac    = f"{SCRATCH}/{{sample}}/atac_filtered.h5ad",
        summary = f"{SCRATCH}/{{sample}}/module1_summary.json",
    params:
        sample_entry = lambda wc: sample_entry(wc.sample),
        qc_params    = config["qc_params"],
        sge_extra    = sge_extra("qc_preprocess_sample"),
    log:
        f"{SCRATCH}/logs/01_qc_preprocess_{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/01_qc_preprocess_{{sample}}.tsv",
    resources:
        runtime = config["resources"]["qc_preprocess_sample"]["runtime_min"],
        sge_extra = sge_extra("qc_preprocess_sample"),
    threads: config["resources"]["qc_preprocess_sample"].get("slots", 1)
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/qc_preprocess_sample.py"


rule aggregate_module1_summary:
    """Combine every sample's per-sample summary JSON into one CSV."""
    input:
        summaries = expand(f"{SCRATCH}/{{sample}}/module1_summary.json", sample=SAMPLES),
    output:
        f"{SCRATCH}/module1_summary.csv",
    log:
        f"{SCRATCH}/logs/01_aggregate_summary.log",
    resources:
        runtime = 15,
        sge_extra = "-V -l h_vmem=2000M",
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/aggregate_module1_summary.py"
