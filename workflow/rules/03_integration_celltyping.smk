# =============================================================================
# workflow/rules/03_integration_celltyping.smk
# Module 3 — dataset integration & cell type labelling
#
# integrate_harmony       — pool all samples, HVG/PCA, Harmony batch correction,
#                            Leiden clustering, UMAP -> integrated.h5ad + a
#                            cluster_template.tsv (cell_type column blank)
# score_markers            — per-cluster expression of cochlear marker genes,
#                            to support the manual annotation call
# apply_cell_type_labels   — MANUAL step in between: copy cluster_template.tsv
#                            to the path in config["integration"]["cluster_annotation_tsv"],
#                            fill in the cell_type column using cluster_umap.png
#                            + cluster_marker_scores.csv, THEN run this rule.
#                            Not part of `rule all` — run explicitly:
#                              snakemake --use-conda --cores 4 \
#                                <scratch>/module3_integration/labeled.h5ad
#
# Conda env: envs/linger_preproc.yaml (scanpy, harmonypy, leidenalg)
# =============================================================================

MODULE3_DIR = f"{SCRATCH}/module3_integration"

rule integrate_harmony:
    input:
        rna_synced = expand(f"{SCRATCH}/{{sample}}/rna_synced.h5ad", sample=SAMPLES),
    output:
        integrated       = f"{MODULE3_DIR}/integrated.h5ad",
        umap_png         = f"{MODULE3_DIR}/cluster_umap.png",
        cluster_template = f"{MODULE3_DIR}/cluster_template.tsv",
    params:
        integration_cfg = INTEGRATION_CFG,
        sample_ids      = SAMPLES,
    log:
        f"{SCRATCH}/logs/03_integrate_harmony.log",
    benchmark:
        f"{SCRATCH}/benchmarks/03_integrate_harmony.tsv",
    resources:
        runtime = config["resources"]["integrate_harmony"]["runtime_min"],
        sge_extra = sge_extra("integrate_harmony"),
    threads: config["resources"]["integrate_harmony"].get("slots", 1)
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/integrate_harmony.py"


rule score_markers:
    input:
        integrated = f"{MODULE3_DIR}/integrated.h5ad",
    output:
        scores_csv  = f"{MODULE3_DIR}/cluster_marker_scores.csv",
        heatmap_png = f"{MODULE3_DIR}/cluster_marker_heatmap.png",
    params:
        marker_genes = INTEGRATION_CFG["marker_genes"],
    log:
        f"{SCRATCH}/logs/03_score_markers.log",
    resources:
        runtime = config["resources"]["score_markers"]["runtime_min"],
        sge_extra = sge_extra("score_markers"),
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/score_markers.py"


rule apply_cell_type_labels:
    """
    NOT part of `rule all` — this is the manual-in-the-loop step. Requires
    you to have already copied cluster_template.tsv to
    config["integration"]["cluster_annotation_tsv"] and filled in every
    cell_type value.
    """
    input:
        integrated = f"{MODULE3_DIR}/integrated.h5ad",
        annotation = INTEGRATION_CFG["cluster_annotation_tsv"],
    output:
        labeled = f"{MODULE3_DIR}/labeled.h5ad",
    log:
        f"{SCRATCH}/logs/03_apply_cell_type_labels.log",
    resources:
        runtime = config["resources"]["apply_cell_type_labels"]["runtime_min"],
        sge_extra = sge_extra("apply_cell_type_labels"),
    conda:
        "../../envs/linger_preproc.yaml"
    script:
        "../scripts/apply_cell_type_labels.py"
