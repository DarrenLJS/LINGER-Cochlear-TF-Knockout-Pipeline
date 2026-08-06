# =============================================================================
# workflow/rules/02_consensus_peaks.smk
# Module 2 — consensus peak set & re-quantification (02a-02d collapsed)
#
# 02a extract_peaks       — per-sample BED from atac_filtered.h5ad
# 02b merge_consensus_peaks — bedtools merge across all samples -> consensus_peaks.bed
# 02c requantify_consensus  — per-sample ATAC re-quantified against consensus peaks
# 02d sync_rna_atac_barcodes — per-sample RNA subset to match consensus-ATAC barcodes,
#                              written to a NEW file (rna_synced.h5ad) rather than
#                              overwriting rna_filtered.h5ad in place — see
#                              sync_rna_atac_barcodes.py docstring for why.
# sanity_checks            — aggregate QC report across all samples, run after 02d
#
# Conda env: envs/linger_preproc.yaml (bedtools, snapatac2, scanpy)
# =============================================================================

rule extract_peaks:
    input:
        atac = f"{SCRATCH}/{{sample}}/atac_filtered.h5ad",
    output:
        bed = f"{SCRATCH}/_peak_bed/{{sample}}.bed",
    log:
        f"{SCRATCH}/logs/02a_extract_peaks_{{sample}}.log",
    resources:
        runtime = config["resources"]["extract_peaks"]["runtime_min"],
        sge_extra = sge_extra("extract_peaks"),
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/extract_peaks.py"


rule merge_consensus_peaks:
    """Concatenate every per-sample BED and merge overlapping intervals with
    `bedtools merge` into one consensus peak set all samples are re-quantified
    against in 02c. Lightweight text sort/merge — safe on a login node."""
    input:
        beds = expand(f"{SCRATCH}/_peak_bed/{{sample}}.bed", sample=SAMPLES),
    output:
        consensus = f"{SCRATCH}/consensus_peaks.bed",
    params:
        merge_gap = config["qc_params"]["merge_gap_bp"],
        concat = f"{SCRATCH}/_all_peaks_concat.bed",
        sorted_bed = f"{SCRATCH}/_all_peaks_sorted.bed",
    log:
        f"{SCRATCH}/logs/02b_merge_consensus_peaks.log",
    resources:
        runtime = config["resources"]["merge_consensus_peaks"]["runtime_min"],
        sge_extra = sge_extra("merge_consensus_peaks"),
    conda:
        "../../envs/linger_preproc.yaml"
    shell:
        r"""
        set -euo pipefail
        exec &> {log}

        echo "Concatenating per-sample BEDs..."
        cat {input.beds} > {params.concat}
        echo "  $(wc -l < {params.concat}) total peak intervals before merge"

        echo "Sorting..."
        sort -k1,1 -k2,2n {params.concat} > {params.sorted_bed}

        echo "Merging (gap={params.merge_gap}bp)..."
        bedtools merge -i {params.sorted_bed} -d {params.merge_gap} > {output.consensus}

        N_CONSENSUS=$(wc -l < {output.consensus})
        echo ""
        echo "Consensus peak set: ${{N_CONSENSUS}} intervals -> {output.consensus}"
        echo ""
        echo "Peak width distribution:"
        awk '{{print $3-$2}}' {output.consensus} | sort -n | awk '
            {{ a[NR]=$1; sum+=$1 }}
            END {{ print "  min="a[1]"  median="a[int(NR/2)]"  max="a[NR]"  mean="sum/NR }}'
        """


rule requantify_consensus:
    input:
        fragments      = lambda wc: sample_entry(wc.sample)["atac_fragments"],
        atac_filtered  = f"{SCRATCH}/{{sample}}/atac_filtered.h5ad",
        consensus_bed  = f"{SCRATCH}/consensus_peaks.bed",
        chrom_sizes    = config["references"]["mm10_chrom_sizes"],
    output:
        atac_consensus = f"{SCRATCH}/{{sample}}/atac_consensus.h5ad",
    log:
        f"{SCRATCH}/logs/02c_requantify_{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/02c_requantify_{{sample}}.tsv",
    resources:
        runtime = config["resources"]["requantify_consensus"]["runtime_min"],
        sge_extra = sge_extra("requantify_consensus"),
    threads: config["resources"]["requantify_consensus"].get("slots", 1)
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/requantify_consensus.py"


rule sync_rna_atac_barcodes:
    input:
        rna_filtered   = f"{SCRATCH}/{{sample}}/rna_filtered.h5ad",
        atac_consensus = f"{SCRATCH}/{{sample}}/atac_consensus.h5ad",
    output:
        rna_synced = f"{SCRATCH}/{{sample}}/rna_synced.h5ad",
    log:
        f"{SCRATCH}/logs/02d_sync_{{sample}}.log",
    resources:
        runtime = config["resources"]["sync_rna_atac_barcodes"]["runtime_min"],
        sge_extra = sge_extra("sync_rna_atac_barcodes"),
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/sync_rna_atac_barcodes.py"


rule sanity_checks:
    input:
        rna_synced     = expand(f"{SCRATCH}/{{sample}}/rna_synced.h5ad", sample=SAMPLES),
        atac_consensus = expand(f"{SCRATCH}/{{sample}}/atac_consensus.h5ad", sample=SAMPLES),
        atac_filtered  = expand(f"{SCRATCH}/{{sample}}/atac_filtered.h5ad", sample=SAMPLES),
        consensus_bed  = f"{SCRATCH}/consensus_peaks.bed",
    output:
        summary = f"{SCRATCH}/qc_sanity_checks/sanity_summary.csv",
    params:
        sample_ids = SAMPLES,
        qc_params  = config["qc_params"],
        qc_dir     = f"{SCRATCH}/qc_sanity_checks",
    log:
        f"{SCRATCH}/logs/02e_sanity_checks.log",
    resources:
        runtime = config["resources"]["sanity_checks"]["runtime_min"],
        sge_extra = sge_extra("sanity_checks"),
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/sanity_checks.py"
