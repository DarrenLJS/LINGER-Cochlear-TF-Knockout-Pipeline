# =============================================================================
# workflow/rules/07_bulk_tf_activity.smk
# Module 7 — TF activity estimation (bulk RNA supplement)
#
# Runs LingerGRN.TF_activity.regulon() — the package's expression-only mode,
# confirmed from the real LingerGRN==1.110 wheel (TF_activity.py: regulon()
# takes only adata_RNA, no ATAC) — over every RNA-only dataset in the
# collection. Two config categories feed this, both RNA-only, both handled
# identically:
#
#   extra_bulk_rna_inputs    (38 entries) -> role="bulk"     -> HTML's
#                                             "~34 bulk/sorted RNA-seq
#                                             datasets" for Module 7 proper.
#   model_construction_refs  (6 entries)  -> role="baseline" -> mechanically
#                                             can ONLY go through regulon()
#                                             (see README's Module 6/7
#                                             design-question note — LINGER_tr
#                                             .training() and perturb.
#                                             load_data_ptb() both require
#                                             paired ATAC, these datasets
#                                             don't have any). Kept in a
#                                             separate role="baseline" output
#                                             bucket so Module 9 can use them
#                                             as the "unperturbed" reference
#                                             point rather than mixing them
#                                             into the bulk/aging TF-activity
#                                             summary.
#
# network="cell population" is used for every dataset here — NOT "general".
# "general" was my first draft and is WRONG: it routes through
# TF_activity.bulk_reg() -> LL_net.load_region(), which is hardcoded to
# hg19/hg38 (branches on `if genome=='hg19'` / `if genome=='hg38'` with no
# mm10 case — confirmed against the real LingerGRN==1.110 source). Same
# human-only hardcoding the README already flagged for the 'LINGER' training
# method, just hit again in a different function.
#
# "cell population" avoids that entirely — it just reads a precomputed
# cell_population_trans_regulatory.txt, and LL_net.trans_reg(method='scNN')
# (which writes that file) is genome-agnostic, mirroring the
# already-verified-working cell_type_specific_trans_reg(). REQUIRES A SMALL
# ADDITION TO MODULE 6: grn_population_training.py currently only calls
# LL_net.TF_RE_binding(), not cis_reg()/trans_reg() — see the patch noted in
# the accompanying README_MODULES_7_10.md. Without that addition,
# cell_population_trans_regulatory.txt won't exist and this whole module
# will fail at the first job.
# =============================================================================

MODULE7_DIR = f"{SCRATCH}/module7_tf_activity"

_BULK_ENTRIES = [dict(e, role="bulk") for e in EXTRA_BULK_RNA_CFG]
_BASELINE_ENTRIES = [dict(e, role="baseline") for e in MODEL_CONSTRUCTION_REFS_CFG]
TF_ACTIVITY_ENTRIES = {
    e.get("sample_id", e.get("ref_id")): e for e in (_BULK_ENTRIES + _BASELINE_ENTRIES)
}


def _tf_activity_entry_json(wildcards):
    import json
    return json.dumps(TF_ACTIVITY_ENTRIES[wildcards.entry_id])


rule linger_tf_activity_one:
    input:
        population_done = f"{SCRATCH}/module6_grn/population_training.done",
        # Written by the grn_population_training.py addition in
        # PATCH_grn_population_training.py.diff — add this path to that
        # rule's own `output:` block in 06_grn_inference.smk too, so
        # Snakemake's DAG actually tracks it instead of just trusting the
        # .done marker.
        trans_regulatory = f"{SCRATCH}/module4_linger_init/cell_population_trans_regulatory.txt",
    output:
        regulon = f"{MODULE7_DIR}/{{entry_id}}.regulon.tsv",
    params:
        workdir  = f"{SCRATCH}/module4_linger_init",   # regulon() reads GRNdir-relative files written by Module 4/6
        grn_dir  = LINGER_CFG["grn_dir"],
        genome   = LINGER_CFG["genome"],
        entry_json = _tf_activity_entry_json,
    log:
        f"{SCRATCH}/logs/07a_tf_activity_{{entry_id}}.log",
    resources:
        runtime   = config["resources"]["linger_tf_activity_one"]["runtime_min"],
        sge_extra = sge_extra("linger_tf_activity_one"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/linger_tf_activity.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --entry-json {params.entry_json:q} \
            --output-tsv {output.regulon}
        """


rule module7_aggregate:
    """Stack every entry's regulon score vector into one wide summary table,
    tagged by role (bulk vs baseline) so Module 9 can pull out the baseline
    subset without re-deriving it."""
    input:
        expand(f"{MODULE7_DIR}/{{entry_id}}.regulon.tsv", entry_id=list(TF_ACTIVITY_ENTRIES.keys())),
    output:
        summary = f"{MODULE7_DIR}/tf_activity_summary.tsv",
    params:
        entries_json = lambda wc: __import__("json").dumps(
            {k: v.get("role") for k, v in TF_ACTIVITY_ENTRIES.items()}
        ),
    log:
        f"{SCRATCH}/logs/07b_aggregate.log",
    resources:
        runtime   = config["resources"]["module7_aggregate"]["runtime_min"],
        sge_extra = sge_extra("module7_aggregate"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        {LINGER_PYTHON} workflow/scripts/aggregate_tf_activity.py \
            --module7-dir {MODULE7_DIR} \
            --roles-json {params.entries_json:q} \
            --output-tsv {output.summary}
        """


rule module7_all:
    input:
        f"{MODULE7_DIR}/tf_activity_summary.tsv",
