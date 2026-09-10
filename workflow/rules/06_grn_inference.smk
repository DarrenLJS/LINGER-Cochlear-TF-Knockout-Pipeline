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

# FIX 2026-09-08 — real cyclic-dependency bug, confirmed via a Snakemake
# dry-run while integrating Module 7: linger_celltype_grn's output pattern
# `{celltype}.done` is unconstrained, so Snakemake considers it a possible
# (wildcard) producer of ANY `.done` file in MODULE6_DIR — including
# `population_training.done` itself (celltype="population_training" is a
# syntactically valid match). That was never triggered before because
# nothing outside Module 6's own internal chain named
# population_training.done as a direct target; Module 7's
# linger_tf_activity_one rule is the first to do so, which made Snakemake
# search for alternate producers and hit the self-referential route
# (linger_celltype_grn(celltype="population_training") requires
# population_training.done as its own input) -> CyclicGraphException.
# Constraining celltype to the real, known values from cluster_annotation.tsv
# closes this off without renaming any already-produced file on Eddie.
import re
wildcard_constraints:
    celltype = "|".join(re.escape(c) for c in CELL_TYPES) if CELL_TYPES else "(?!)"


rule linger_population_training:
    input:
        pseudobulk_done = f"{SCRATCH}/module4_linger_init/pseudobulk.done",
        tss_done        = f"{SCRATCH}/module4_linger_init/tss_redist.done",
        motif_bed       = f"{SCRATCH}/module4_linger_init/MotifTarget.bed",
        labeled         = f"{SCRATCH}/module3_integration/labeled.h5ad",
        atac_consensus  = expand(f"{SCRATCH}/{{sample}}/atac_consensus.h5ad", sample=SAMPLES),
    output:
        done = f"{MODULE6_DIR}/population_training.done",
        # Added 2026-09-08 — grn_population_training.py's cis_reg()/trans_reg()
        # addition (see script's own 2026-09-07 docstring note) writes these
        # two files but they were only ever tracked implicitly via `done`.
        # Declaring them explicitly lets Snakemake's DAG treat them as real
        # outputs — needed by Module 7's linger_tf_activity_one rule, which
        # depends on cell_population_trans_regulatory.txt directly.
        cis_regulatory   = f"{SCRATCH}/module4_linger_init/cell_population_cis_regulatory.txt",
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_population_trans_regulatory.txt",
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
        export PYTHONHASHSEED=0
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
        # Added 2026-09-08 — same reasoning as linger_population_training's
        # output: block above. grn_celltype_specific.py already writes these
        # three files per celltype (confirmed in its own docstring); they
        # just weren't declared. Module 10's benchmark_grn_edges rule reads
        # the cis + TF-RE files by this exact path pattern already — this
        # just makes Snakemake's DAG aware of the real producer, so reruns/
        # cache invalidation behave correctly instead of relying on the
        # files already existing on disk from a prior manual run.
        cis_regulatory    = f"{SCRATCH}/module4_linger_init/cell_type_specific_cis_regulatory_{{celltype}}.txt",
        trans_regulatory  = f"{SCRATCH}/module4_linger_init/cell_type_specific_trans_regulatory_{{celltype}}.txt",
        tf_re_binding     = f"{SCRATCH}/module4_linger_init/cell_type_specific_TF_RE_binding_{{celltype}}.txt",
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
        export PYTHONHASHSEED=0
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
