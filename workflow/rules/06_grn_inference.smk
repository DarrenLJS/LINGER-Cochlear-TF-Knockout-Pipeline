# =============================================================================
# workflow/rules/06_grn_inference.smk
# Module 6 — GRN inference
#
# linger_population_training — LINGER_tr.training(method='scNN') +
#                               LL_net.TF_RE_binding(method='scNN')
# linger_celltype_grn         — per-cell-type cis/trans regulatory networks,
#                               wildcarded over CELL_TYPES (below)
#
# model_construction_refs (GSE135703/202920/283534/111349/136196/168041) are
# NOT yet wired into a rule here — LingerGRN's public API doesn't expose an
# obvious "extra baseline expression" argument to training()/cell_type_specific_*
# the way the HTML's Module 6 description implies. Confirm intended usage
# before I build this in: most likely they're meant to extend adata_RNA
# before pseudobulking (Module 4), not modify Module 6 itself — flagging
# rather than guessing.
#
# CELL_TYPES is read from cluster_annotation_tsv at DAG-build time. That file
# only exists after Module 3's manual apply_cell_type_labels step, so this
# whole module is unreachable (and Snakemake will say so clearly) until
# that's done — matches the HTML's dependency ordering.
#
# Conda env: NOT via Snakemake's `conda:` directive as of 2026-07-25 — see
# 04_linger_init.smk header for why (SGE misclassified the env directory
# path as a name). Both rules below call LINGER_PYTHON directly.
# =============================================================================

MODULE6_DIR = f"{SCRATCH}/module6_grn"

def _cell_types():
    ann_path = INTEGRATION_CFG["cluster_annotation_tsv"]
    if not os.path.exists(ann_path):
        return []
    df = pd.read_csv(ann_path, sep="\t")
    return sorted(df["cell_type"].dropna().unique().tolist())

CELL_TYPES = _cell_types()


rule linger_population_training:
    input:
        pseudobulk_done = f"{SCRATCH}/module4_linger_init/pseudobulk.done",
        tss_done        = f"{SCRATCH}/module4_linger_init/tss_redist.done",
        motif_bed       = f"{SCRATCH}/module4_linger_init/MotifTarget.bed",
        labeled         = f"{SCRATCH}/module3_integration/labeled.h5ad",
        atac_consensus  = expand(f"{SCRATCH}/{{sample}}/atac_consensus.h5ad", sample=SAMPLES),
    output:
        done = f"{MODULE6_DIR}/population_training.done",
    params:
        workdir    = f"{SCRATCH}/module4_linger_init",   # same workdir as Module 4 — see script docstring
        grn_dir    = LINGER_CFG["grn_dir"],
        genome     = LINGER_CFG["genome"],
        activef    = LINGER_CFG["activef"],
        sample_ids = SAMPLES,
    log:
        f"{SCRATCH}/logs/06a_population_training.log",
    benchmark:
        f"{SCRATCH}/benchmarks/06a_population_training.tsv",
    resources:
        runtime   = config["resources"]["linger_population_training"]["runtime_min"],
        sge_extra = sge_extra("linger_population_training"),
    threads: config["resources"]["linger_population_training"].get("slots", 1)
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        # Same PATH/LD_LIBRARY_PATH restoration as 04_linger_init.smk's rules —
        # LINGER_PYTHON is called by absolute path (bypassing `conda activate`),
        # so the env's own bin/ (bedtools/intersectBed for pybedtools) and
        # lib/ (conda-forge libstdc++, avoids GLIBCXX_3.4.30 not found on
        # older-libstdc++ Eddie nodes) never land on PATH/LD_LIBRARY_PATH
        # otherwise. See 04_linger_init.smk module docstring for the full story.
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/grn_population_training.py \
            --workdir "{params.workdir}" \
            --grn-dir "{params.grn_dir}" \
            --genome "{params.genome}" \
            --activef "{params.activef}" \
            --motif-bed "{input.motif_bed}" \
            --labeled "{input.labeled}" \
            --atac-consensus {input.atac_consensus} \
            --sample-ids {params.sample_ids} \
            --output-done "{output.done}"
        """


rule linger_celltype_grn:
    input:
        population_done = f"{MODULE6_DIR}/population_training.done",
        labeled          = f"{SCRATCH}/module3_integration/labeled.h5ad",
        atac_consensus   = expand(f"{SCRATCH}/{{sample}}/atac_consensus.h5ad", sample=SAMPLES),
    output:
        done = f"{MODULE6_DIR}/{{celltype}}.done",
    params:
        workdir    = f"{SCRATCH}/module4_linger_init",
        grn_dir    = LINGER_CFG["grn_dir"],
        genome     = LINGER_CFG["genome"],
        sample_ids = SAMPLES,
    log:
        f"{SCRATCH}/logs/06b_celltype_grn_{{celltype}}.log",
    resources:
        runtime   = config["resources"]["linger_celltype_grn"]["runtime_min"],
        sge_extra = sge_extra("linger_celltype_grn"),
    threads: config["resources"]["linger_celltype_grn"].get("slots", 1)
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        # Same PATH/LD_LIBRARY_PATH restoration as 04_linger_init.smk's rules —
        # see linger_population_training above and 04_linger_init.smk's
        # module docstring for why this is needed.
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/grn_celltype_specific.py \
            --workdir "{params.workdir}" \
            --grn-dir "{params.grn_dir}" \
            --genome "{params.genome}" \
            --celltype "{wildcards.celltype}" \
            --labeled "{input.labeled}" \
            --atac-consensus {input.atac_consensus} \
            --sample-ids {params.sample_ids} \
            --output-done "{output.done}"
        """


rule module6_all:
    """Convenience target: all annotated cell types' GRNs. Not (yet) in
    `rule all` — Module 6 depends on the manual Module 3 label step, so it
    can't be part of the default DAG until CELL_TYPES is non-empty."""
    input:
        expand(f"{MODULE6_DIR}/{{celltype}}.done", celltype=CELL_TYPES),
