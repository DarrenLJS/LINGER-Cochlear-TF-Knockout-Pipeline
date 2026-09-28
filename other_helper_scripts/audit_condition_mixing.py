#!/usr/bin/env python3
"""
audit_condition_mixing.py — Stage B structural audit for the
condition-mixing bug found in GSE281207 (Module 9 negative control) and
confirmed in several other entries on the first run of this script
(GSE224627, GSE233559, GSE153882, GSE154833, GSE196870, GSE273939 in
Module 9; GSE150002, GSE116702, GSE139145, GSE149257, GSE198167, GSE242321,
GSE268504, GSE275243, GSE294736, GSE296324, GSE312476, GSE90821 in Module 7).

THE BUG BEING AUDITED FOR: every one of these config entries gets loaded
via bulk_rna_loader.load_rna_only() into one AnnData, then collapsed to a
SINGLE mean vector (linger_tf_activity.py's `score.mean(axis=1)`, or
validate_held_out.py's `adata.X.mean(axis=0)` / `real_tf_activity.mean(axis=1)`)
before being compared against anything. If an entry's samples actually span
more than one biological condition (genotype, treatment, age), that mean is
a blend, not a real "this dataset's profile."

THREE INDEPENDENT WAYS this shows up, all checked here:
  1. obs_names encode multiple conditions directly (bulk/pseudobulk data,
     e.g. GSE281207's "P21-HET-A"/"P21-KO-B" — the loader keeps real sample
     names as obs_names for every bulk format). Checked via TOKEN_FAMILIES.
  2. The config entry itself already declares >1 raw source file merged
     into one entry (multi_file_merge's `files:` list, or a glob_pattern
     matching >1 file with different condition tokens in their names — the
     GSE242321/GSE256069 shape). Checked via describe_config_multiplicity().
  3. Single-cell/single-nucleus entries where obs_names are just barcodes
     (e.g. "AAACCCAAG...-1") and carry NO condition info at all — the
     condition, if recoverable, lives either in adata.obs columns already
     baked into a real .h5ad file, or in a sibling per-cell metadata file
     the loader never reads (bulk_rna_loader's mtx-trio branch only ever
     reads matrix+barcodes+features, nothing else, even if a metadata file
     sits right next to them). Checked via scan_obs_columns() and
     find_possible_metadata_files() — this is what resolves the "unknown"
     cases (GSE274279 aging reference vector, and the barcode-based
     model_construction_refs entries) that TOKEN_FAMILIES alone cannot,
     since there is no obs_names to scan.
  4. A dataset's directory holds MORE raw files than the loader ever
     reads — bulk_rna_loader._find_one() globs non-recursively and returns
     only the first match when a pattern hits more than one file. Any
     entry not already using a multi-file-aware format silently drops
     every sample after the first, with no trace in the loaded AnnData at
     all (not blended, just gone). Checked via check_silent_first_file_risk().
     This is why a barcode-based entry reporting "no condition signal
     found" is not automatically safe to trust — it could mean the dataset
     genuinely has one sample, or it could mean several GSM files were on
     disk and only the first was ever loaded.

KNOWN FIXED BUGS IN EARLIER VERSIONS OF THIS SCRIPT, both found by running
it against the real data and checking the printed obs_names by eye rather
than trusting the `status` column:
  - v1: a family flagged only when >=2 *different regex patterns* matched
    (e.g. both "wt" and "ko"). Missed GSE154833 ("1mo"/"9mo"/"26mo" — three
    ages, one pattern `\d+mo`). Fixed in v2 to count distinct MATCHED
    SUBSTRINGS, not distinct patterns.
  - v2: still vocabulary-bound, so it missed GSE153882 ("9m"/"26m", not
    "9mo"), GSE196870 (7 groups mixing month/day/hour units plus two
    unrelated experiment arms with no shared vocabulary), and GSE273939
    (arbitrary treatment codes "Ele"/"Thr"/"Twe" that no keyword list could
    ever anticipate). No fixed vocabulary generalizes to free-form
    dataset-specific naming, so v3 adds a vocabulary-FREE structural check
    (cluster_by_prefix, below): strip each obs_name's trailing
    replicate/well-index suffix and cluster what's left. >=2 distinct
    prefix groups means >=2 conditions, regardless of what they're named.
    This is what actually catches Ele/Thr/Twe, and is now the PRIMARY
    signal; TOKEN_FAMILIES stays as a secondary, human-readable "why" when
    it does match a known vocabulary. Validated against every real case
    seen across two audit runs (all of the confirmed-mixed datasets above,
    plus deliberately-single-condition cases — plain WT replicates, "ctrl1
    .. ctrl5", "sample_1 .. sample_3" — to check it does NOT over-flag).
    Known remaining blind spot: it assumes a replicate suffix is numeric/
    letter/well-id-shaped; a dataset whose true replicates carry distinct
    WORD suffixes (e.g. "Sample_Alpha"/"Sample_Beta" as arbitrary IDs for
    the same condition) would still false-flag. Read the printed obs_names
    yourself for anything this flags before treating it as certain.

This script does NOT fix anything. It loads every dataset EXACTLY the way
the real pipeline does (same loader, same entries, same config) and prints
the real ground truth (obs_names, obs columns, and sibling files) before
anything gets averaged away.

Needs the real downloaded GEO files, which live only on Eddie's scratch
space, not in the code checkout — run it there.

Run on Eddie from the pipeline's code directory (same convention as
other_helper_scripts/verify_loader_regression.py):

    cd /exports/eddie/scratch/s2906787/linger_pipeline/code
    /exports/csce/eddie/biology/groups/bioinfmsc/anaconda/envs/s2906787/LINGER/bin/python \\
        other_helper_scripts/audit_condition_mixing.py

(LINGER_PYTHON, not linger_preproc — bulk_rna_loader.py's lazy imports of
scanpy/h5py/openpyxl/xlrd are already guaranteed present there since
Module 7's own script runs this exact loader under that same interpreter.)

Do NOT `conda activate` first — call the interpreter by absolute path, same
reasoning as every other Module 4/6/7/8/9 script in this repo (see
Snakefile's top-of-file note on why `--use-conda` handles activation and
manual activation is not the convention here).

Writes audit_condition_mixing_report.tsv in the current directory and also
prints everything to stdout — pipe through `tee` if you want both a saved
copy and to watch it live:

    ... audit_condition_mixing.py 2>&1 | tee audit_condition_mixing.log
"""
import argparse
import glob
import os
import re
import sys
import traceback

import yaml

# Condition axes to scan obs_names against. Each axis is a set of named
# BUCKETS, and each bucket is a set of alternate spellings for the SAME
# semantic value (e.g. bucket "wildtype" covers "wt19", "WT-1", "wild-type").
# A dataset is flagged on an axis if >=2 DIFFERENT BUCKETS actually match
# among its own obs_names — e.g. "wildtype" AND "mutant" both present.
# Bucketing by MEANING rather than literal substring matters here: a naive
# "does >=2 distinct matched strings appear" check would wrongly flag a
# dataset that's just WT-1/WT-2/WT-3 (three distinct strings, one real
# condition), and would wrongly MISS a dataset that's WT vs Mut (each only
# ever contributes one distinct string on its own side, so a same-pattern
# check never sees 2 values unless wt/mut are recognized as one axis).
# Patterns allow an optional trailing digit run (`wt19`, `mut21`, no
# separator) as well as separator-joined replicate numbers (`WT-1`).
GENOTYPE_AXIS = {
    "wildtype":     [r"\bwt\d*\b", r"wild.?type"],
    "knockout":     [r"\bko\d*\b", r"knock.?out", r"-/-", r"\bnull\b"],
    "heterozygous": [r"\bhet\d*\b", r"heterozygous"],
    "homozygous":   [r"\bhomo\d*\b", r"homozygous"],
    "mutant":       [r"\bmut\d*\b", r"mutant"],
    "control":      [r"\bctrl\d*\b", r"\bctl\d*\b", r"\bcontrol\b"],
}
TREATMENT_AXIS = {
    "treated":   [r"treat(?!ment\b)\w*"],
    "untreated": [r"untreat\w*"],
    "chir":      [r"\bchir\b"],
    "dapt":      [r"\bdapt\b"],
    "tamoxifen": [r"tamoxifen"],
    "vehicle":   [r"\bvehicle\b"],
}
# age_or_timepoint has no small fixed vocabulary — every distinct matched
# literal (e.g. "1mo" vs "9mo" vs "26mo", or "p21" vs "p70") is its own
# bucket, extracted directly rather than pre-enumerated.
AGE_PATTERNS = [r"\bp\d{1,3}\b", r"\bd\d{1,3}\b", r"\d+\s*mo(?![a-z])", r"e\d{1,2}\.\d"]

# obs_names that look like this are single-cell/nucleus barcodes, not real
# sample names — scanning them against TOKEN_FAMILIES is meaningless (no
# condition info is ever encoded in a 10x barcode). Route these to the
# obs-column / sibling-file checks instead.
_BARCODE_RE = re.compile(r"^[ACGTN]{10,20}(-\d+)?$", re.IGNORECASE)

# Filename fragments that suggest a per-cell/per-sample metadata file the
# mtx-trio loader branch would never read (it only ever opens
# matrix/barcodes/features, nothing else) — a real candidate for where a
# condition/age/genotype label is hiding when obs_names are bare barcodes.
_METADATA_FILENAME_HINTS = [
    "meta", "annotation", "cluster", "celltype", "cell_type", "sample_info",
    "coldata", "cellinfo", "cell_info", "design", "phenotype", "pdata",
    "characteristics", "colData",
]


def resolve(obj, path_roots):
    if isinstance(obj, str):
        return obj.format(**path_roots)
    if isinstance(obj, dict):
        return {k: resolve(v, path_roots) for k, v in obj.items()}
    if isinstance(obj, list):
        return [resolve(v, path_roots) for v in obj]
    return obj


def looks_like_barcodes(obs_names, threshold=0.8):
    if not obs_names:
        return False
    n_matching = sum(1 for n in obs_names if _BARCODE_RE.match(str(n)))
    return (n_matching / len(obs_names)) >= threshold


def _normalize_for_matching(obs_names):
    # NOTE: Python's \b treats '_' as a word character, so "1mo_sv_a" has
    # NO boundary between "mo" and "_" and \bmo\b-style patterns silently
    # miss it. Real sample names here are underscore-joined constantly
    # ("1mo_SV_A", "E18.5_Insm1-Het_OHC1") — normalize '_' to '-' (a real
    # separator for \b) before matching. Display copy elsewhere is unaffected.
    return [str(n).lower().replace("_", "-") for n in obs_names]


def scan_token_families(obs_names):
    """Flags an axis (genotype, treatment) when >=2 different semantic
    BUCKETS are present — not >=2 literal strings, see module docstring on
    why that distinction matters. age_or_timepoint has no fixed bucket
    vocabulary, so it flags on >=2 distinct literal matched values instead
    (e.g. "1mo" vs "9mo")."""
    joined = _normalize_for_matching(obs_names)
    hits = {}

    for axis_name, buckets in (("genotype", GENOTYPE_AXIS), ("treatment", TREATMENT_AXIS)):
        present_buckets = set()
        for bucket, patterns in buckets.items():
            rx_list = [re.compile(pat, re.IGNORECASE) for pat in patterns]
            if any(rx.search(n) for n in joined for rx in rx_list):
                present_buckets.add(bucket)
        if len(present_buckets) >= 2:
            hits[axis_name] = sorted(present_buckets)

    age_values = set()
    for pat in AGE_PATTERNS:
        rx = re.compile(pat, re.IGNORECASE)
        for n in joined:
            age_values.update(m.lower() for m in rx.findall(n))
    if len(age_values) >= 2:
        hits["age_or_timepoint"] = sorted(age_values)

    return hits


def prefix_of(name):
    """Strip a trailing replicate/well-index suffix from one sample name,
    greedily, token by token from the end. What survives is treated as the
    sample's "condition label" for clustering purposes — see
    cluster_by_prefix() and the module docstring for why this exists and
    what it does and doesn't catch."""
    tokens = [t for t in re.split(r"[_\-\s]+", str(name)) if t != ""]
    changed = True
    while changed and tokens:
        changed = False
        last = tokens[-1]
        if len(tokens) > 1 and re.fullmatch(r"\d+", last):
            tokens.pop(); changed = True; continue
        if len(tokens) > 1 and re.fullmatch(r"[A-Za-z]", last):
            tokens.pop(); changed = True; continue
        if len(tokens) > 1 and re.fullmatch(r"[A-Za-z]\d{1,3}", last):
            tokens.pop(); changed = True; continue  # well ids: A01, B12
        if len(tokens) > 1 and re.fullmatch(r"rep\d*", last, re.IGNORECASE):
            tokens.pop(); changed = True; continue
        if len(tokens) > 1 and last.lower() in ("value", "rpkm", "fpkm", "reads", "counts"):
            tokens.pop(); changed = True; continue
        # fused alpha+digit tail (e.g. "wt19", "mut21", "ele4") — strip the
        # digit suffix in place without discarding the token, so it can
        # still collapse to its alphabetic condition label even when there
        # was no separator to split on in the first place.
        m = re.fullmatch(r"([A-Za-z]+)(\d+)", last)
        if m and m.group(1):
            tokens[-1] = m.group(1)
            changed = True
            continue
    return "_".join(tokens).lower() if tokens else str(name).lower()


def cluster_by_prefix(obs_names):
    """The primary, vocabulary-free condition-mixing signal — see module
    docstring. Returns {prefix: [original names]} and the caller flags
    when len(groups) >= 2."""
    groups = {}
    for n in obs_names:
        groups.setdefault(prefix_of(n), []).append(n)
    return groups


def scan_obs_columns(adata):
    """For single-cell entries, obs_names (barcodes) carry no condition
    info — but if the real source was a .h5ad that already had a
    condition/age/genotype column baked into .obs, it survives the load
    untouched (bulk_rna_loader's h5ad branch is a bare ad.read_h5ad(), it
    doesn't strip anything). Report every obs column's cardinality, and
    flag+detail any column with a plausible "condition axis" shape (more
    than 1, but not absurdly many, distinct values — a real per-cell QC
    float column would have thousands of unique values and gets skipped).
    `barcode` is excluded — load_rna_only() always adds it itself from
    obs_names, so it is never informative here."""
    report_lines = []
    flagged = False
    for col in adata.obs.columns:
        if col == "barcode":
            continue
        try:
            n_unique = adata.obs[col].nunique(dropna=True)
        except TypeError:
            continue  # unhashable column contents, skip
        if 2 <= n_unique <= 30:
            counts = adata.obs[col].value_counts(dropna=False).to_dict()
            report_lines.append(f"obs['{col}'] ({n_unique} values): {counts}")
            flagged = True
        else:
            report_lines.append(f"obs['{col}'] ({n_unique} values, not condition-shaped)")
    return "; ".join(report_lines), flagged


_MULTI_SAMPLE_LOADER_FORMATS = {"cellranger_h5_merge", "multi_file_merge", "per_sample_merge"}
_RAW_FILE_PATTERNS = {
    "matrix.mtx": ["*matrix.mtx.gz", "*matrix.mtx"],
    "barcodes": ["*barcodes.tsv.gz", "*barcodes.tsv"],
    "features_or_genes": ["*features.tsv.gz", "*genes.tsv.gz", "*features.tsv", "*genes.tsv"],
    "h5": ["*.h5"],
    "h5ad": ["*.h5ad"],
}


def check_silent_first_file_risk(entry):
    """bulk_rna_loader._find_one() globs NON-recursively and returns only
    hits[0] when a pattern matches more than one file — confirmed by
    reading the function directly. For any entry whose format reaches
    _find_one via an AMBIGUOUS, directory-wide wildcard (i.e. NOT already
    pinned to one exact file), if its directory actually contains more
    than one raw file matching the same pattern _find_one would use, every
    sample/GSM after the first is silently dropped — not blended, not
    visible in obs_names, not detectable from the loaded AnnData at all.
    This is a different bug from condition-mixing (dropping data outright
    vs. averaging it in), and is exactly why "no condition signal found in
    the loaded data" is not the same claim as "this dataset only has one
    condition" — it could also mean "the other conditions were never
    loaded in the first place."

    THREE ways an entry can already be pinned to one exact file/pattern,
    none of which is actually ambiguous even though _find_one is still the
    function that runs — all three must be exempted or this check
    false-positives (confirmed on the first run: GSE182202_P1_RNA and both
    GSE157398 entries use `rna_prefix`, which _load_mtx_trio() resolves via
    exact _first_existing() checks, not a directory glob at all; GSE120462
    _GFP/_Ikzf2 each set their OWN entry-specific `glob_patterns` naming
    one exact filename apiece, which IS what _find_one uses for them, so
    each entry only ever matches its own file):
      1. `rna_prefix` set — _load_mtx_trio uses exact existence checks
         against `prefix + "matrix.mtx.gz"` etc., never a wildcard.
      2. `glob_pattern` (singular) set — the per_sample_merge /
         multi_file_merge-style formats that already iterate every match
         on purpose (also covered by _MULTI_SAMPLE_LOADER_FORMATS below).
      3. `glob_patterns` (PLURAL) set — _load_cellranger_h5's own
         entry-specific override, a single-file allowlist, not a
         directory-wide default.
    Returns (is_at_risk, {pattern_name: n_matching_files}) for every
    pattern with >1 match."""
    if (entry.get("format") in _MULTI_SAMPLE_LOADER_FORMATS
            or entry.get("glob_pattern")
            or entry.get("glob_patterns")
            or entry.get("rna_prefix")):
        return False, {}  # already pinned to an exact file/pattern, not ambiguous
    directory = entry.get("path") or os.path.dirname(entry.get("rna_prefix", "") or "")
    if not directory or not os.path.isdir(directory):
        return False, {}
    multi = {}
    for label, patterns in _RAW_FILE_PATTERNS.items():
        n = 0
        for pat in patterns:
            n += len(glob.glob(os.path.join(directory, pat)))
        if n > 1:
            multi[label] = n
    return bool(multi), multi


def find_possible_metadata_files(entry):
    """Filesystem-level check, independent of what the loader actually
    read: does this entry's directory contain any file that LOOKS like a
    per-cell/per-sample metadata table the mtx-trio branch would silently
    ignore? This is what turns a "we can't tell" into either "yes, and
    here it is, go read it" or "no such file exists, the loader isn't
    hiding anything.\""""
    directory = entry.get("path")
    if not directory or not os.path.isdir(directory):
        return []
    hits = []
    for fname in sorted(os.listdir(directory)):
        low = fname.lower()
        if any(hint in low for hint in _METADATA_FILENAME_HINTS):
            hits.append(fname)
    return hits


def describe_config_multiplicity(entry):
    """Independent of obs_names: does the CONFIG ITSELF already declare
    more than one raw source file for this one entry? (multi_file_merge's
    `files:` list, or a glob_pattern that resolves to >1 file on disk).
    This is exactly the shape GSE242321 (Atf6-fpkm + WT-fpkm) and
    GSE296324 (KO_summary + wildtype_summary) both have."""
    notes = []
    files = entry.get("files")
    if files:
        patterns = [f.get("pattern") for f in files if isinstance(f, dict)]
        notes.append(f"multi_file_merge with {len(files)} file patterns: {patterns}")
    glob_pattern = entry.get("glob_pattern")
    directory = entry.get("path")
    if glob_pattern and directory:
        try:
            hits = sorted(glob.glob(os.path.join(directory, glob_pattern)))
            notes.append(f"glob_pattern {glob_pattern!r} matches {len(hits)} file(s): "
                         f"{[os.path.basename(h) for h in hits]}")
        except Exception as e:
            notes.append(f"glob_pattern check failed: {e}")
    return "; ".join(notes)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config_eddie.yaml")
    ap.add_argument("--scripts-dir", default="workflow/scripts",
                     help="Where bulk_rna_loader.py lives, relative to repo root")
    ap.add_argument("--out", default="audit_condition_mixing_report.tsv")
    args = ap.parse_args()

    sys.path.insert(0, os.path.abspath(args.scripts_dir))
    from bulk_rna_loader import load_rna_only  # noqa: E402

    with open(args.config) as f:
        config = yaml.safe_load(f)

    path_roots = {
        "multiome_root": config["multiome_root"],
        "chrom_root": config["chrom_root"],
        "scratch": config["scratch"],
    }

    def section(name, wrap_single=False, role=None):
        # `role` mirrors what the real Snakemake rules inject before
        # calling load_rna_only (07_bulk_tf_activity.smk: role="bulk" for
        # extra_bulk_rna_inputs, role="baseline" for model_construction_refs)
        # — bulk_rna_loader._load_cellranger_h5 branches on entry["role"]
        # to collapse to a pseudobulk row for role=="bulk" entries
        # (GSE120462_GFP/Ikzf2 specifically). Without injecting it here,
        # this audit would load those two as raw 737k-barcode AnnData
        # instead of the real pipeline's single-row pseudobulk — confirmed
        # to matter on the second run (that mismatch was part of why they
        # first showed up as barcode-like at all).
        raw = config.get(name)
        if raw is None:
            return []
        if wrap_single:
            raw = [raw]
        entries = [dict(resolve(e, path_roots), _category=name) for e in raw]
        if role:
            for e in entries:
                e["role"] = role
        return entries

    # generalization_test_panel entries use `test_id`/`gene` rather than
    # sample_id/ref_id — load_rna_only() itself only ever looks at path/
    # format/glob_pattern-style keys, so this reuses the same loader path,
    # but the `sample_id = entry.get("sample_id") or entry.get("ref_id")`
    # line further down needs a matching id here too, or every row prints
    # "None" instead of the real test_id. Handled via a small pre-pass
    # rather than touching the shared `sample_id` lookup, so behavior for
    # every pre-existing category is unchanged.
    gen_panel = section("generalization_test_panel")
    for e in gen_panel:
        e.setdefault("sample_id", e.get("test_id"))

    entries = (
        section("extra_bulk_rna_inputs", role="bulk")
        + section("model_construction_refs", role="baseline")
        + section("held_out_inputs")
        + section("negative_control_input", wrap_single=True)
        + gen_panel
    )
    print(f"Auditing {len(entries)} config entries "
          f"(extra_bulk_rna_inputs + model_construction_refs + held_out_inputs + "
          f"negative_control_input + generalization_test_panel) for within-entry "
          f"condition mixing.\n"
          f"NOTE: generalization_test_panel entries ({len(gen_panel)} of them) point at "
          f"setup_scripts/download_datasets.sh Section 10's download paths, which must "
          f"exist on disk already (run that section first) — and their path/format/"
          f"glob_pattern in config_eddie.yaml are still placeholders pending this exact "
          f"audit, so expect to iterate on config_eddie.yaml based on what this run finds, "
          f"not to have it pass clean on the first try.\n")

    rows = []
    for entry in entries:
        sample_id = entry.get("sample_id") or entry.get("ref_id")
        category = entry["_category"]
        row = {
            "sample_id": sample_id, "dataset": entry.get("dataset"), "category": category,
            "n_samples": None, "looks_like_barcodes": None, "obs_names": "",
            "prefix_clusters": "", "token_family_flags": "", "obs_column_report": "",
            "possible_unread_metadata_files": "", "config_multiplicity_notes": "",
            "silent_first_file_risk": "", "status": "OK", "error": "",
        }
        try:
            config_notes = describe_config_multiplicity(entry)
            row["config_multiplicity_notes"] = config_notes

            metadata_files = find_possible_metadata_files(entry)
            row["possible_unread_metadata_files"] = "; ".join(metadata_files)

            silent_risk, silent_counts = check_silent_first_file_risk(entry)
            row["silent_first_file_risk"] = (
                "; ".join(f"{k}={v} files, only 1st loaded" for k, v in silent_counts.items())
                if silent_risk else ""
            )

            adata = load_rna_only(entry)
            obs_names = list(adata.obs_names)
            row["n_samples"] = len(obs_names)
            row["obs_names"] = "; ".join(obs_names)

            barcode_like = looks_like_barcodes(obs_names)
            row["looks_like_barcodes"] = barcode_like

            token_hits = scan_token_families(obs_names) if not barcode_like else {}
            row["token_family_flags"] = (
                "; ".join(f"{fam}={vals}" for fam, vals in token_hits.items())
                if token_hits else ""
            )

            # Primary signal — vocabulary-free, see module docstring. Skip
            # for barcode-like names: they have no shared "condition
            # prefix" by construction, clustering them would just flag
            # every single-cell dataset spuriously.
            prefix_groups = cluster_by_prefix(obs_names) if not barcode_like else {}
            prefix_flag = len(prefix_groups) >= 2
            row["prefix_clusters"] = (
                "; ".join(f"{p}({len(v)})" for p, v in sorted(prefix_groups.items()))
                if prefix_flag else ""
            )

            obs_col_report, obs_col_flag = scan_obs_columns(adata)
            row["obs_column_report"] = obs_col_report

            # silent_risk takes priority over everything else below: it
            # means data may be MISSING entirely (other GSM files never
            # loaded), which makes "no condition signal found" an unsafe
            # conclusion regardless of what the other checks say.
            if silent_risk:
                row["status"] = "UNRESOLVED_SILENT_FIRST_FILE_ONLY"
            elif token_hits or prefix_flag or config_notes or obs_col_flag:
                row["status"] = "FLAGGED_LIKELY_MIXED"
            elif barcode_like and metadata_files:
                row["status"] = "UNRESOLVED_CHECK_METADATA_FILE_MANUALLY"
            elif barcode_like:
                row["status"] = "barcode_based_no_condition_signal_found"
            else:
                row["status"] = "looks_single_condition"

            print(f"[{row['status']:>34}] {sample_id} ({category}, n={len(obs_names)}"
                  f"{', barcode-like' if barcode_like else ''}): "
                  f"{obs_names[:12]}{' ...' if len(obs_names) > 12 else ''}")
            if silent_risk:
                print(f"                                     SILENT-DROP RISK: {row['silent_first_file_risk']} "
                      f"in {entry.get('path')!r} — this dataset may have more samples on disk than were loaded")
            if prefix_flag:
                print(f"                                     prefix clusters: {row['prefix_clusters']}")
            if token_hits:
                print(f"                                     token flags: {row['token_family_flags']}")
            if obs_col_flag:
                print(f"                                     obs columns: {obs_col_report}")
            if config_notes:
                print(f"                                     config notes: {config_notes}")
            if metadata_files:
                print(f"                                     unread candidate metadata file(s): {metadata_files}")

        except Exception as e:
            row["status"] = "ERROR"
            row["error"] = f"{type(e).__name__}: {e}"
            print(f"[{'ERROR':>34}] {sample_id} ({category}): {row['error']}")
            traceback.print_exc(file=sys.stderr)

        rows.append(row)

    import csv
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_flagged = sum(r["status"] == "FLAGGED_LIKELY_MIXED" for r in rows)
    n_unresolved_meta = sum(r["status"] == "UNRESOLVED_CHECK_METADATA_FILE_MANUALLY" for r in rows)
    n_unresolved_silent = sum(r["status"] == "UNRESOLVED_SILENT_FIRST_FILE_ONLY" for r in rows)
    n_error = sum(r["status"] == "ERROR" for r in rows)
    print(f"\n{n_flagged}/{len(rows)} entries flagged as likely mixing conditions; "
          f"{n_unresolved_meta}/{len(rows)} unresolved (barcode-based, with a candidate "
          f"metadata file sitting unread in the directory — go look at it by hand); "
          f"{n_unresolved_silent}/{len(rows)} unresolved (directory holds more raw files than "
          f"the loader reads — samples may be silently missing, go look at it by hand); "
          f"{n_error}/{len(rows)} failed to load (see errors above/log).")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
