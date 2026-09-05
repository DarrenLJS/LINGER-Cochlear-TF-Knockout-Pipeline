# =============================================================================
# workflow/rules/04_linger_init.smk
# Module 4 — LINGER initialisation (mouse mode)
#
# linger_prep_pseudobulk — labeled.h5ad + atac_consensus.h5ad -> LINGER's
#                           pseudobulk inputs (data/TG_pseudobulk.tsv etc)
# linger_get_tss          — LINGER_tr.get_TSS + RE_TG_dis against GRNdir
# linger_motif_scan       — HOMER motif scan over consensus_peaks.bed ->
#                            MotifTarget.bed (external tool, not a LingerGRN
#                            function — required input to Module 6's
#                            LL_net.TF_RE_binding(method='scNN'))
#
# Depends on Module 3's MANUAL labeled.h5ad (apply_cell_type_labels), since
# LingerGRN.pseudo_bulk.pseudo_bulk() requires cell-type labels before it can
# pseudobulk — confirmed from source, see linger_prep_pseudobulk.py docstring.
#
# Conda env: NOT via Snakemake's `conda:` directive as of 2026-07-25 — that
# directive misclassified the pre-built LINGER env's directory path as a
# NAME under real SGE execution ('conda env export --name' errored on the
# path's slashes). linger_prep_pseudobulk/linger_get_tss now call that
# env's own python binary directly by absolute path via `shell:` instead,
# sidestepping Snakemake's conda resolution entirely. linger_motif_scan
# still uses no conda env at all — HOMER is a standalone Perl install.
# =============================================================================

MODULE4_DIR = f"{SCRATCH}/module4_linger_init"
LINGER_CFG = config["linger"]
LINGER_PYTHON = f"{LINGER_CFG['conda_env_linger']}/bin/python"
LINGER_ENV_BIN = f"{LINGER_CFG['conda_env_linger']}/bin"
# Calling LINGER_PYTHON directly (see module docstring above) also skips the
# LD_LIBRARY_PATH setup `conda activate` normally does, not just PATH.
# Conda-forge builds (e.g. scipy's compiled _highs_wrapper.so, pulled in via
# scipy.optimize <- scipy.stats <- LingerGRN.LINGER_tr) link against the
# env's own bundled libstdc++, which is newer than some Eddie compute nodes'
# system /lib64/libstdc++.so.6. Without this, the rule fails with
# "GLIBCXX_3.4.30 not found" — but only on older-libstdc++ nodes, so it
# doesn't reproduce every run. Confirmed 2026-09-05 (node1f16 failed this
# way; node1n14 didn't, same env, same code).
LINGER_ENV_LIB = f"{LINGER_CFG['conda_env_linger']}/lib"

rule linger_prep_pseudobulk:
    input:
        labeled        = f"{SCRATCH}/module3_integration/labeled.h5ad",
        atac_consensus = expand(f"{SCRATCH}/{{sample}}/atac_consensus.h5ad", sample=SAMPLES),
    output:
        done = f"{MODULE4_DIR}/pseudobulk.done",
    params:
        workdir    = MODULE4_DIR,
        sample_ids = SAMPLES,
    log:
        f"{SCRATCH}/logs/04a_prep_pseudobulk.log",
    resources:
        runtime   = config["resources"]["linger_prep_pseudobulk"]["runtime_min"],
        sge_extra = sge_extra("linger_prep_pseudobulk"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        # linger_get_tss/linger_prep_pseudobulk call LINGER_PYTHON by absolute
        # path instead of via Snakemake's `conda:` directive (see module
        # docstring above — SGE misclassified the env dir as an env name
        # under `conda:`). That sidesteps `conda activate`, which normally
        # also puts the env's bin/ on PATH — so binaries the env installs
        # (bedtools/intersectBed for pybedtools, etc) silently vanish from
        # PATH even though the env has them. Restore that explicitly, same
        # pattern linger_motif_scan already uses for homer_bin below.
        export PATH="{LINGER_ENV_BIN}:$PATH"
        # See LINGER_ENV_LIB comment above — restores conda-forge's compiled
        # extensions' expected libstdc++, same reasoning as PATH above.
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/linger_prep_pseudobulk.py \
            --labeled {input.labeled} \
            --atac-consensus {input.atac_consensus} \
            --sample-ids {params.sample_ids} \
            --workdir {params.workdir} \
            --output-done {output.done}
        """


rule linger_get_tss:
    input:
        pseudobulk_done = f"{MODULE4_DIR}/pseudobulk.done",
    output:
        done = f"{MODULE4_DIR}/tss_redist.done",
    params:
        workdir        = MODULE4_DIR,
        grn_dir        = LINGER_CFG["grn_dir"],
        genome         = LINGER_CFG["genome"],
        tss_distance_bp = LINGER_CFG["tss_distance_bp"],
    log:
        f"{SCRATCH}/logs/04b_get_tss.log",
    resources:
        runtime   = config["resources"]["linger_get_tss"]["runtime_min"],
        sge_extra = sge_extra("linger_get_tss"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        # RE_TG_dis() shells out to `intersectBed` via pybedtools. bedtools
        # is a real dependency of the LINGER env (envs/linger.yaml), but
        # since we invoke LINGER_PYTHON by absolute path instead of via
        # `conda:` (see module docstring above), conda activate never
        # runs and the env's bin/ never lands on PATH — so intersectBed
        # can't be found even though it's installed right where it should
        # be. Restore it explicitly, same pattern linger_motif_scan
        # already uses for homer_bin below.
        export PATH="{LINGER_ENV_BIN}:$PATH"
        # See LINGER_ENV_LIB comment above — restores conda-forge's compiled
        # extensions' expected libstdc++, same reasoning as PATH above.
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/linger_get_tss.py \
            --workdir {params.workdir} \
            --grn-dir {params.grn_dir} \
            --genome {params.genome} \
            --tss-distance-bp {params.tss_distance_bp} \
            --output-done {output.done}
        """


rule linger_motif_scan:
    """
    HOMER motif scan over the consensus peak set -> MotifTarget.bed, the
    input LL_net.load_TFbinding_scNN() reads for TF-RE binding scores
    (Module 6). External tool, not a LingerGRN function, and NOT a conda
    package in this setup — setup_linger.sh installs HOMER as a standalone
    Perl tool at linger.homer_bin, added to ~/.bashrc, not via conda. No
    `conda:` directive here for that reason — this rule runs in the
    ambient shell with homer_bin prepended to PATH explicitly.

    CONFIRMED 2026-07-24 against real Eddie output (previously a guess):
    - HOMER mm10 genome package IS installed (`configureHomer.pl -list`
      shows `+  mm10  v7.0`).
    - motif_file is `all_motif_rmdup_Mammal`, not the placeholder
      `all_motif_rmdup.motif` from the earlier draft — GRNdir splits this
      by taxon (`_fly`, `_Mammal`, `_Plant`); mouse uses `_Mammal`.
    - GRNdir also has `Match_TF_motif_Mus_musculus.txt` and
      `genome_map_homer.txt` — confirms the full scNN chain resolves:
      genome='mm10' -> genome_map_homer.txt -> species='Mus_musculus'
      -> Match_TF_motif_Mus_musculus.txt (read internally by
      LL_net.load_TFbinding_scNN() in Module 6, not this rule).
    - references/jaspar/'s raw JASPAR PFMs are NOT used anywhere in this
      chain — nothing in LINGER's scNN code path reads .jaspar files
      directly. Left alone; not referenced by any rule.
    """
    input:
        consensus_bed = f"{SCRATCH}/consensus_peaks.bed",
    output:
        motif_bed = f"{MODULE4_DIR}/MotifTarget.bed",
    params:
        genome     = LINGER_CFG["genome"],
        homer_bin  = LINGER_CFG["homer_bin"],
        motif_file = f"{LINGER_CFG['grn_dir']}/all_motif_rmdup_Mammal",
    log:
        f"{SCRATCH}/logs/04c_motif_scan.log",
    resources:
        runtime   = config["resources"]["linger_motif_scan"]["runtime_min"],
        sge_extra = sge_extra("linger_motif_scan"),
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        export PATH="{params.homer_bin}:$PATH"
        annotatePeaks.pl {input.consensus_bed} {params.genome} \
            -m {params.motif_file} \
            -mbed {output.motif_bed} \
            > /dev/null
        echo "Wrote {output.motif_bed}"
        """
