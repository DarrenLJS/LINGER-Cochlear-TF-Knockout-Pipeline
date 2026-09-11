# =============================================================================
# workflow/rules/10_grn_benchmarking.smk
# Module 10 — GRN benchmarking (Module 5's Option A, decided)
#
# Hi-C loops (GSE305205) and CUT&RUN peaks (GSE150386/150391/181307) are used
# purely as independent ground truth here, against Module 6's ALREADY-WRITTEN
# per-cell-type edge files — no changes to Modules 4/6, no separate Module 5
# preprocessing stage. This is deliberately NOT calling LingerGRN.Benchmk.
# bm_trans() — that function expects a ChIP-seq ranked-gene-list ground
# truth format (skiprows=5, 'symbol'/'score' columns), not chromatin
# interaction data, so a bespoke bedtools-intersect + sklearn AUC/AUPR
# script is the right tool here, not a package function that doesn't fit.
#
# Two edge types benchmarked, matching the HTML's Module 10 description:
#   cis   : cell_type_specific_cis_regulatory_{celltype}.txt (RE,TG,score)
#           vs GSE305205 Hi-C loop anchors
#   TF-RE : cell_type_specific_TF_RE_binding_{celltype}.txt (RE x TF matrix)
#           vs GSE150386/150391/181307 CUT&RUN peaks, matched by TF identity
#           (Atoh1 for GSE150386, Pou4f3 for GSE150391; GSE181307 has no TF
#           named in chromatin_prior_extra's role field — CONFIRM which TF
#           it's for before trusting that comparison; benchmarked against
#           all TFs as a fallback with a loud warning if unconfirmed).
# =============================================================================

MODULE10_DIR = f"{SCRATCH}/module10_benchmarking"


rule benchmark_grn_edges:
    input:
        # Every annotated cell type's cis + TF-RE binding files from Module 6.
        cis_files = expand(
            f"{SCRATCH}/module4_linger_init/cell_type_specific_cis_regulatory_{{celltype}}.txt",
            celltype=CELL_TYPES,
        ),
        tfre_files = expand(
            f"{SCRATCH}/module4_linger_init/cell_type_specific_TF_RE_binding_{{celltype}}.txt",
            celltype=CELL_TYPES,
        ),
    output:
        report = f"{MODULE10_DIR}/benchmark_report.tsv",
    params:
        celltypes = CELL_TYPES,
        chrom_priors_json = lambda wc: __import__("json").dumps(EXTRA_CHROMATIN_PRIOR_CFG_FULL),
        module4_dir = f"{SCRATCH}/module4_linger_init",
        outdir = MODULE10_DIR,
    log:
        f"{SCRATCH}/logs/10_benchmark.log",
    resources:
        runtime   = config["resources"]["benchmark_grn_edges"]["runtime_min"],
        sge_extra = sge_extra("benchmark_grn_edges"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{LINGER_ENV_BIN}:$PATH"
        # FIX 2026-09-11 — real bug, confirmed via disk inspection + per-log
        # trace: celltype names here contain spaces and hyphens (e.g.
        # "Non-sensory epithelium - unresolved Gata3+"). Without :q,
        # Snakemake stringifies params.celltypes by joining elements with
        # spaces and NO per-element quoting, so the shell then word-splits
        # each multi-word celltype name into several separate --celltypes
        # tokens (argparse nargs="+" happily swallows all of them). This
        # is exactly why benchmark_chromatin_priors.py's own log showed
        # nonsense fragment "celltypes" like "-", "unresolved", "Gata3+"
        # instead of the real names, and why every file lookup inside it
        # missed (real per-celltype files DO exist on disk under the
        # correct cell_type_specific_TF_RE_binding_{{celltype}}.txt name —
        # confirmed via `find`). :q shell-quotes each element individually,
        # so a name containing spaces survives as one argument.
        {LINGER_PYTHON} workflow/scripts/benchmark_chromatin_priors.py \
            --module4-dir {params.module4_dir} \
            --celltypes {params.celltypes:q} \
            --chrom-priors-json {params.chrom_priors_json:q} \
            --outdir {params.outdir} \
            --output-report {output.report}
        """
