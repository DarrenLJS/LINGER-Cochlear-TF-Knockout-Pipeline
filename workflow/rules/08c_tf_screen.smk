# =============================================================================
# workflow/rules/08c_tf_screen.smk
# Module 8c — all-TF in-silico screen (ADDED 2026-10-01)
#
# Every TF in the 586-TF population list is knocked out (mode ko) and/or
# overexpressed (mode oe, model-capped) ONE AT A TIME, in one or more scopes,
# using linger_perturbation.py's batch mode (each chromosome net loaded once,
# baseline + all variants of the chunk through it as one stacked batch).
#
# WHAT IT IS FOR
#   * an empirical null for Module 9/9b: "how does THIS knockout's agreement
#     with the real shift compare with a typical single-TF knockout's?" (the
#     negative control no longer rests on a fixed threshold);
#   * the candidate set for a known-aging-factor ranking (dimension 6) and for
#     matched-random comparisons.
#
# OUTPUT per scope (module8c_tf_screen/<scope>/):
#   chunk_NNN_shift.tsv.gz   genes x variants, mean over samples of
#                            [perturbed - baseline] prediction (+ a
#                            "__baseline_mean__" column); ids "<mode>__<TF>"
#   chunk_NNN_status.tsv     per variant: n_samples_changed, shift_l2, noop
#   chunk_NNN_variants.json  the variant list that chunk ran
#   screen_status.tsv        all chunks merged + tf_in_network / usable flags
#
# RESOURCES ARE A FIRST GUESS: a single-KO job measured ~55 s / 5.3 GB, but a
# chunk runs (1 + variants) forward passes per gene, so its compute scales with
# the variant count. RUN ONE CHUNK FIRST and read its [timing] lines before
# submitting the rest (see the README / delivery notes).
# =============================================================================

import re

MODULE8C_DIR = f"{SCRATCH}/module8c_tf_screen"
_SCREEN_CFG = config.get("tf_screen", {})
SCREEN_SCOPES = [] if LEGACY_FORWARD else list(_SCREEN_CFG.get("scopes", ["population"]))
SCREEN_MODES = list(_SCREEN_CFG.get("modes", ["ko", "oe"]))
SCREEN_N_CHUNKS = int(_SCREEN_CFG.get("n_chunks", 20))
assert SCREEN_MODES and all(m in ("ko", "oe") for m in SCREEN_MODES), f"tf_screen.modes={SCREEN_MODES!r}"
assert SCREEN_N_CHUNKS >= 1, "tf_screen.n_chunks must be >= 1"
_bad_scopes = [s for s in SCREEN_SCOPES if s != "population" and s not in CELLTYPE_KO_TYPES]
if _bad_scopes:
    raise ValueError(f"tf_screen.scopes {_bad_scopes} are not 'population' or an enabled Module 8b cell type "
                     f"({CELLTYPE_KO_TYPES})")
SCREEN_CHUNKS = [f"{i:03d}" for i in range(SCREEN_N_CHUNKS)]

wildcard_constraints:
    screen_scope = "|".join(re.escape(s) for s in SCREEN_SCOPES) if SCREEN_SCOPES else "(?!)",
    screen_chunk = r"\d{3}",


def _screen_inputs(wc):
    d = {
        "exp": f"{MODULE8_DIR}/Exp.tsv",
        "re_tglink": f"{MODULE8_DIR}/RE_TGlink_resolved.tsv",
    }
    if wc.screen_scope == "population":
        d["sanity_check"] = f"{MODULE8_DIR}/_sanity_check_baseline_predicted.tsv"
    else:
        ct = wc.screen_scope
        d["sanity_check"] = f"{MODULE8B_DIR}/{ct}/_sanity_check_baseline_predicted.tsv"
        d["tg"] = f"{MODULE8B_DIR}/{ct}/TG_pseudobulk_{ct}.tsv"
        d["re_"] = f"{MODULE8B_DIR}/{ct}/RE_pseudobulk_{ct}.tsv"
    return d


def _screen_scope_args(wc, input):
    if wc.screen_scope == "population":
        return ""
    return f'--target-path "{input.tg}" --opn-path "{input.re_}"'


rule tf_screen_chunk:
    """One chunk of the all-TF screen: ~SCREEN_N_CHUNKS-th of the TFs, every
    configured mode, in one scope. The sanity-check output is an input only so
    the sanity check has run (and its log exists) first — same discipline as
    Modules 8/8b."""
    input:
        unpack(_screen_inputs),
    output:
        shift = f"{MODULE8C_DIR}/{{screen_scope}}/chunk_{{screen_chunk}}_shift.tsv.gz",
        status = f"{MODULE8C_DIR}/{{screen_scope}}/chunk_{{screen_chunk}}_status.tsv",
        variants = f"{MODULE8C_DIR}/{{screen_scope}}/chunk_{{screen_chunk}}_variants.json",
    params:
        workdir = f"{SCRATCH}/module4_linger_init",
        module8_dir = MODULE8_DIR,
        scope_args = _screen_scope_args,
        modes = " ".join(SCREEN_MODES),
    log:
        f"{SCRATCH}/logs/08d_tf_screen_{{screen_scope}}_{{screen_chunk}}.log",
    resources:
        runtime   = config["resources"]["tf_screen_chunk"]["runtime_min"],
        sge_extra = sge_extra("tf_screen_chunk"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        export PYTHONHASHSEED=0
        {LINGER_PYTHON} workflow/scripts/make_screen_variants.py \
            --exp "{input.exp}" \
            --chunk {wildcards.screen_chunk} --n-chunks {SCREEN_N_CHUNKS} \
            --modes {params.modes} \
            --out "{output.variants}"
        {LINGER_PYTHON} workflow/scripts/linger_perturbation.py \
            --workdir {params.workdir} \
            --module8-dir {params.module8_dir} \
            --ko-id "_screen_{wildcards.screen_scope}_{wildcards.screen_chunk}" \
            --variants-json "{output.variants}" \
            --oe-quantile {OE_QUANTILE} \
            {PERTURB_FLAGS} \
            {params.scope_args} \
            --shift-tsv "{output.shift}" \
            --status-tsv "{output.status}"
        """


def _screen_trans_reg(wc):
    if wc.screen_scope == "population":
        return f"{SCRATCH}/module4_linger_init/cell_population_trans_regulatory.txt"
    return f"{SCRATCH}/module4_linger_init/cell_type_specific_trans_regulatory_{wc.screen_scope}.txt"


rule tf_screen_summarise:
    input:
        status = expand(f"{MODULE8C_DIR}/{{{{screen_scope}}}}/chunk_{{chunk}}_status.tsv", chunk=SCREEN_CHUNKS),
        trans_regulatory = _screen_trans_reg,
    output:
        f"{MODULE8C_DIR}/{{screen_scope}}/screen_status.tsv",
    params:
        screen_dir = f"{MODULE8C_DIR}/{{screen_scope}}",
    log:
        f"{SCRATCH}/logs/08d_tf_screen_summarise_{{screen_scope}}.log",
    resources:
        runtime   = config["resources"]["tf_screen_summarise"]["runtime_min"],
        sge_extra = sge_extra("tf_screen_summarise"),
    shell:
        r"""
        set -euo pipefail
        exec &> "{log}"
        export PATH="{LINGER_ENV_BIN}:$PATH"
        export LD_LIBRARY_PATH="{LINGER_ENV_LIB}:$LD_LIBRARY_PATH"
        {LINGER_PYTHON} workflow/scripts/summarise_tf_screen.py \
            --screen-dir "{params.screen_dir}" \
            --n-chunks {SCREEN_N_CHUNKS} \
            --trans-regulatory-tsv "{input.trans_regulatory}" \
            --output-tsv "{output}"
        """


rule module8c_all:
    input:
        expand(f"{MODULE8C_DIR}/{{screen_scope}}/screen_status.tsv", screen_scope=SCREEN_SCOPES),
