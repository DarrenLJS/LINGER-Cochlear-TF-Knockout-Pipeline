"""
bulk_rna_loader.py — shared loader for Module 7's RNA-only datasets.

Covers BOTH config categories that feed Module 7 (both are RNA-only, so both
go through TF_activity.regulon(), see 07_bulk_tf_activity.smk header):
  - extra_bulk_rna_inputs (38 entries)    -> role="bulk"
  - model_construction_refs (6 entries)   -> role="baseline"

REWRITTEN 2026-09-08 after a real investigation across all 24 of the original
44 datasets that failed against the first version of this loader (see
README's Module 7 development-log entry for the full per-dataset writeup).
That investigation confirmed real, on-disk layouts across FIVE distinct new
shapes beyond the original three, and confirmed — by reading TF_activity.py's
real regulon() source, not by assumption — that regulon()'s own internal
per-sample proportion + quantile normalization makes the exact input unit
(raw counts vs. TPM vs. FPKM) far less consequential than originally
assumed, as long as what's kept is genuine per-sample expression and not
metadata/statistics columns. That finding is why the fixes below prioritize
correctly EXCLUDING non-sample columns over choosing one "correct" unit.

AUDIT LOGGING: every call logs, to the calling task's own log file (already
redirected via `exec &> {log}` in the rule's shell block), which loader
branch served it and — for every branch that filters columns — exactly
which columns were kept vs. dropped. This exists specifically so a full
clean rerun is auditable against the previous run's real output, not just
trusted by reasoning about the code.

Auto-detected layouts, in try-order (unless a `format:` override is set):
  1. 10x mtx trio (matrix.mtx[.gz] + barcodes.tsv[.gz] + features/genes.tsv[.gz])
  2. A single tab/comma-delimited counts matrix file, genes-as-rows.
     Non-numeric and numeric-but-non-sample columns are now excluded rather
     than crashing the whole matrix — see _coerce_numeric_columns().
  3. A single .h5ad file.
  4. Excel (.xlsx or disguised/real .xls, optionally .gz) — tries a real
     Excel parse first, falls back to delimited-text on failure, since this
     investigation found BOTH real binary .xls (GSE262205, confirmed via
     OLE2 magic bytes) and fake .xls-labeled plain text (GSE271610,
     confirmed via hex dump) in the same batch of datasets.

FORMAT OVERRIDES (config_eddie.yaml `format:` key), for shapes no
auto-detector should guess at:
  - per_sample_merge: N sibling per-sample files (e.g. RSEM genes.results)
    merged column-wise into one matrix. Required: glob_pattern,
    value_column. Optional: id_column (default "gene_id"),
    rollup ("sum_by_gene" for isoform-level files, requires
    rollup_id_column), sample_name_regex (else derives sample name from
    the filename directly).
  - multi_file_merge: N condition-level files (e.g. GSE242321's
    Atf6-fpkm/WT-fpkm split, GSE296324's KO_summary/wildtype_summary
    split), each contributing its own subset of sample columns, merged
    column-wise. Required: files (list of {pattern, sample_column_regex
    or exclude_columns}).
  - cellranger_h5: CellRanger's own (possibly non-standard-genome-keyed)
    HDF5 export, read via scanpy.read_10x_h5(). Optional: h5_genome_key
    for older exports that don't use a literal genome name Scanpy expects.
  - cellranger_h5_merge: N sibling per-sample CellRanger .h5 files in one
    directory, each collapsed to a pseudobulk column and merged into a
    multi-sample matrix. Required: glob_pattern (specific enough to
    exclude any manifest/.tar siblings in the same directory). Optional:
    sample_name_source ("library_id" default, reads each file's own h5
    attr; "filename" uses the basename instead).

sheet_name (config key, "excel" auto-detected path only): overrides the
default of reading the workbook's first sheet, for workbooks whose first
sheet isn't the real data (e.g. a README/notes tab — confirmed real case
GSE154833, sheet order ['READ ME', 'All values', 'Means', 'Diff Exp']).

sample_column_regex / exclude_columns (config keys, usable with the
"matrix"/"excel" auto-detected paths and inside multi_file_merge's per-file
specs): when a file mixes real per-sample columns with metadata or
statistics columns that automatic dtype-coercion can't safely distinguish
(e.g. a numeric gene-length column, or FPKM/count/reads triplicated per
sample), set sample_column_regex to an anchored regex matching only the
real sample columns, or exclude_columns to a denylist of specific names to
drop before dtype-coercion runs on what's left.
"""
import ast
import glob
import gzip
import io
import os
import re

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def _log(msg):
    # The calling rule's shell already redirects stdout to the task's own
    # log file (`exec &> {log}` in 07_bulk_tf_activity.smk), so a plain
    # print is sufficient and keeps this module dependency-free.
    print(f"[bulk_rna_loader] {msg}")


def _is_gz(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def _smart_open(path):
    """Binary-mode open that checks real gzip magic bytes (\\x1f\\x8b) rather
    than trusting a '.gz' extension. Confirmed necessary because some GEO
    submitters name a file '*.gz' without ever actually running gzip on it
    (e.g. GSE157398's bulk-RNA barcodes/features files) — pandas/scipy's
    extension-based compression inference fails loudly on these with
    'Not a gzipped file'. Safe on genuinely gzipped files too."""
    return gzip.open(path, "rb") if _is_gz(path) else open(path, "rb")


def _smart_open_text(path):
    """Text-mode counterpart of _smart_open, for line-by-line scanning."""
    return gzip.open(path, "rt") if _is_gz(path) else open(path, "rt")


def _find_one(directory, patterns):
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(directory, pat)))
        if hits:
            return hits[0]
    return None


def _first_existing(paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return None


def _detect_header_row(path, id_column, lookahead=5):
    """
    Some GEO per-sample-merge files prepend a caption line (and sometimes a
    blank line after it) before the real tab-delimited header — confirmed
    real shape for GSE312476, where 10 of 11 sibling files read:
        Bulk RNA-Seq anlaysis from Figure 1H
        <blank line>
        gene_symbol\tgene_id\t<sample>
    while the 11th file (control_1) has the real header on line 1 with no
    caption at all. A fixed skiprows count would break on whichever shape
    it wasn't tuned for. Instead, scan the first `lookahead` raw lines
    (unaffected by pandas' blank-line handling) and return the index of the
    first one that actually contains id_column as a tab-separated token.
    """
    with _smart_open_text(path) as f:
        for i in range(lookahead):
            line = f.readline()
            if not line:
                break
            if id_column in line.rstrip("\n").split("\t"):
                return i
    raise ValueError(
        f"per_sample_merge: could not find a header row containing "
        f"id_column {id_column!r} within the first {lookahead} lines of "
        f"{path!r}. Inspect the file directly — it may need `header: "
        f"false` or a different `id_column` override."
    )


# --------------------------------------------------------------------------
# Branch 1: 10x mtx trio — UNCHANGED from the original loader. Confirmed via
# Snakemake's own `code` rerun-trigger (hashes the rule's shell command, not
# this module) plus this function being untouched that any of the original
# 20 successes routed through here are unaffected by this rewrite.
# --------------------------------------------------------------------------
def _load_mtx_trio(directory, prefix=None):
    matrix = _find_one(directory, ["*matrix.mtx.gz", "*matrix.mtx"]) if not prefix else _first_existing(
        [prefix + "matrix.mtx.gz", prefix + "matrix.mtx"]
    )
    barcodes = _find_one(directory, ["*barcodes.tsv.gz", "*barcodes.tsv"]) if not prefix else _first_existing(
        [prefix + "barcodes.tsv.gz", prefix + "barcodes.tsv"]
    )
    features = _find_one(directory, ["*features.tsv.gz", "*genes.tsv.gz", "*features.tsv", "*genes.tsv"]) if not prefix else _first_existing(
        [prefix + "features.tsv.gz", prefix + "genes.tsv.gz", prefix + "features.tsv", prefix + "genes.tsv"]
    )
    if not (matrix and barcodes and features):
        return None

    # NOTE 2026-09-09: previously only `matrix` went through gzip-magic
    # detection here; `barcodes`/`features` went straight through
    # pd.read_csv, which infers compression from the '.gz' extension alone
    # and crashes with BadGzipFile on files that are plain text wearing a
    # '.gz' name (confirmed real for GSE157398_Apex_RNA's barcodes file).
    # All three now go through the same _smart_open() magic-byte check.
    mat = sio.mmread(_smart_open(matrix)).tocsr()
    bc = pd.read_csv(_smart_open(barcodes), header=None, sep="\t")[0].values
    ft = pd.read_csv(_smart_open(features), header=None, sep="\t")
    gene_ids = ft[0].values if ft.shape[1] == 1 else ft[1].values

    if mat.shape[0] == len(gene_ids) and mat.shape[1] == len(bc):
        X = mat.T
    elif mat.shape[0] == len(bc) and mat.shape[1] == len(gene_ids):
        X = mat
    else:
        raise ValueError(
            f"mtx trio in {directory}: matrix shape {mat.shape} matches "
            f"neither (genes={len(gene_ids)}, cells={len(bc)}) orientation."
        )
    adata = ad.AnnData(X=sp.csr_matrix(X))
    adata.obs_names = [str(b) for b in bc]
    adata.var_names = [str(g) for g in gene_ids]
    _log(f"branch=mtx_trio matrix={os.path.basename(matrix)} shape={adata.shape}")
    return adata


# --------------------------------------------------------------------------
# Shared column-selection helpers
# --------------------------------------------------------------------------
def _coerce_numeric_columns(df, min_columns=2, deny_list=None, source_name=""):
    """
    Category A auto-fix: keep columns that parse as numeric with low NaN
    fraction, drop everything else. Refuses (raises) rather than silently
    returning a too-thin matrix — this is the universal safety net that
    turns any future unrecognized file shape into a loud failure instead
    of a silent 1-sample or metadata-as-sample matrix.
    """
    deny = set(deny_list or [])
    kept, dropped = [], []
    numeric_cols = {}
    for col in df.columns:
        if col in deny:
            dropped.append(col)
            continue
        coerced = pd.to_numeric(df[col], errors="coerce")
        if coerced.isna().mean() <= 0.05:
            numeric_cols[col] = coerced.fillna(0)
            kept.append(col)
        else:
            dropped.append(col)

    if len(kept) < min_columns:
        raise ValueError(
            f"[{source_name}] Only {len(kept)} numeric column(s) survived "
            f"coercion (kept={kept}, dropped={dropped}) — refusing to "
            f"build a matrix with fewer than {min_columns} samples. This "
            f"usually means the picked file isn't a real multi-sample "
            f"matrix (e.g. a DE-results/statistics table), or needs an "
            f"explicit `sample_column_regex` override naming the real "
            f"sample columns."
        )
    _log(f"[{source_name}] dtype-coercion kept={kept} dropped={dropped}")
    return pd.DataFrame(numeric_cols)


def _select_by_regex(df, regex, source_name=""):
    pat = re.compile(regex)
    cols = [c for c in df.columns if pat.match(str(c))]
    if not cols:
        raise ValueError(
            f"[{source_name}] sample_column_regex {regex!r} matched no "
            f"columns. Columns present: {list(df.columns)}"
        )
    dropped = [c for c in df.columns if c not in cols]
    _log(f"[{source_name}] sample_column_regex={regex!r} kept={cols} dropped={dropped}")
    return df[cols].apply(pd.to_numeric, errors="coerce").fillna(0)


def _select_columns(df, entry, source_name=""):
    """Single entry point used by every column-bearing branch: regex
    override takes priority, else deny-list + auto dtype-coercion."""
    regex = entry.get("sample_column_regex")
    if regex:
        return _select_by_regex(df, regex, source_name=source_name)
    deny_list = entry.get("exclude_columns", [])
    return _coerce_numeric_columns(df, deny_list=deny_list, source_name=source_name)


def _resolve_var_names(df, index_values, entry, source_name=""):
    """
    NOTE 2026-09-09: TF_activity.regulon() intersects RNA.index (i.e.
    whatever this loader sets as var_names) against trans_reg.index —
    confirmed real gene-SYMBOL space (e.g. '2700069I18Rik'), matching
    Module 6's cell_population_trans_regulatory.txt. If var_names ends up
    in a different identifier space (Ensembl gene_id, Agilent probe ID,
    etc.), that intersection is EMPTY and regulon() doesn't error — it
    just silently returns no TF activity for that dataset, which is what
    happened for every per_sample_merge/single_matrix/excel dataset whose
    `id_column`/index happened to be Ensembl rather than symbol.

    `symbol_column` (new config key): when the source file has a real gene-
    symbol column separate from whatever's being used as the row index/id
    (e.g. GSE203228's `geneSymbol` alongside `geneID`, GSE85519's `sym`
    alongside `geneid`), set var_names from THAT column instead of
    `index_values`. The id/index itself is left alone for row-identity/
    dedup purposes elsewhere in the caller — only the final var_names
    assignment changes.

    `df` must still have symbol_column as an actual column (i.e. called
    BEFORE any column-filtering step that would drop it).
    """
    symbol_column = entry.get("symbol_column") if entry else None
    if not symbol_column:
        return [str(g) for g in index_values]
    if symbol_column not in df.columns:
        raise ValueError(
            f"[{source_name}] symbol_column {symbol_column!r} not found — "
            f"columns present: {list(df.columns)}"
        )
    _log(f"[{source_name}] var_names set from symbol_column={symbol_column!r} (id/index kept for row-identity only)")
    return [str(g) for g in df[symbol_column].values]


# --------------------------------------------------------------------------
# Branch 2: single matrix file — Category A/A-variant, now with real column
# filtering instead of blindly trusting every column after index_col=0.
# --------------------------------------------------------------------------
def _load_single_matrix_file(directory, entry=None):
    hit = _find_one(
        directory,
        ["*counts*.txt*", "*counts*.csv*", "*counts*.tsv*", "*expression*.txt*",
         "*expression*.csv*", "*expression*.tsv*", "*.txt.gz", "*.txt", "*.csv.gz", "*.csv", "*.tsv.gz", "*.tsv"],
    )
    if hit is None:
        return None
    entry = entry or {}
    sep = "," if hit.endswith((".csv", ".csv.gz")) else "\t"
    df = pd.read_csv(hit, sep=sep, index_col=0)
    numeric_df = _select_columns(df, entry, source_name=os.path.basename(hit))
    adata = ad.AnnData(X=sp.csr_matrix(numeric_df.values.T))
    adata.obs_names = [str(c) for c in numeric_df.columns]
    adata.var_names = _resolve_var_names(df, df.index, entry, source_name=os.path.basename(hit))
    _log(f"branch=single_matrix file={os.path.basename(hit)} shape={adata.shape}")
    return adata


def _load_h5ad(directory):
    hit = _find_one(directory, ["*.h5ad"])
    if hit is None:
        return None
    adata = ad.read_h5ad(hit)
    _log(f"branch=h5ad file={os.path.basename(hit)} shape={adata.shape}")
    return adata


# --------------------------------------------------------------------------
# Branch 4: Excel — real .xlsx/.xls, or a disguised delimited-text file
# wearing an .xls extension (confirmed both exist in this dataset batch).
# --------------------------------------------------------------------------
def _load_excel(directory, entry=None):
    hit = _find_one(directory, ["*.xlsx", "*.xls.gz", "*.xls"])
    if hit is None:
        return None
    entry = entry or {}
    raw = gzip.open(hit, "rb").read() if hit.endswith(".gz") else None

    # sheet_name (config key): workbooks with multiple sheets previously
    # always parsed sheet_names[0], which silently grabbed a README/notes
    # tab instead of the real data tab for GSE154833
    # (['READ ME', 'All values', 'Means', 'Diff Exp'] — sheet 0 was the
    # README, confirmed via openpyxl inspection). An explicit override
    # names the real data sheet; if it doesn't exist in the workbook we
    # raise immediately rather than falling through to the delimited-text
    # fallback below, which would otherwise try to parse binary xlsx bytes
    # as text and produce a confusing, unrelated error.
    override_sheet = entry.get("sheet_name")

    try:
        buf = io.BytesIO(raw) if raw is not None else hit
        xl = pd.ExcelFile(buf)
        if override_sheet is not None and override_sheet not in xl.sheet_names:
            raise ValueError(
                f"[{os.path.basename(hit)}] sheet_name override {override_sheet!r} "
                f"not found — sheets present: {xl.sheet_names}"
            )
        sheet_name = override_sheet or xl.sheet_names[0]
        df = xl.parse(sheet_name, index_col=0)
        _log(
            f"[{os.path.basename(hit)}] read as real Excel, sheet={sheet_name!r}"
            + (" (explicit override)" if override_sheet else " (default: first sheet)")
        )
    except Exception as e:
        if override_sheet is not None:
            raise
        _log(f"[{os.path.basename(hit)}] read_excel failed ({e!r}); falling back to delimited-text read")
        buf = io.BytesIO(raw) if raw is not None else hit
        df = pd.read_csv(buf, sep="\t", index_col=0)

    numeric_df = _select_columns(df, entry, source_name=os.path.basename(hit))
    adata = ad.AnnData(X=sp.csr_matrix(numeric_df.values.T))
    adata.obs_names = [str(c) for c in numeric_df.columns]
    adata.var_names = _resolve_var_names(df, df.index, entry, source_name=os.path.basename(hit))
    _log(f"branch=excel file={os.path.basename(hit)} shape={adata.shape}")
    return adata


# --------------------------------------------------------------------------
# format: per_sample_merge — N sibling per-sample files -> one matrix.
# Confirmed real shapes this covers: RSEM genes.results (GSE149254,
# GSE150002_P1SC_reprogramming, GSE268504), RSEM isoforms.results with
# gene-level rollup (GSE149257), and two simpler single-value-column
# per-sample formats (GSE139145, GSE312476).
# --------------------------------------------------------------------------
def _load_per_sample_merge(directory, entry):
    """
    id_column: name of the file's real, LABELED id column (e.g. RSEM's
    "gene_id", or GSE275243/GSE83599's confirmed named first column) — read
    normally, then df.set_index(id_column).

    id_column omitted/None: the file's first column has NO header name at
    all — confirmed real shape for GSE116702 (header row has one fewer
    field than data rows) and, combined with header=False below, for
    GSE198167/GSE294736 (no header row whatsoever, HTSeq-count-style).
    Parsed with index_col=0, which correctly picks up the first column as
    the id either way.

    header: True (default) for a real header row; False for a confirmed
    genuinely headerless file (GSE198167/GSE294736) — pandas would
    otherwise silently consume the first data row as column names, which
    is exactly what produced the misleading "kept=['260']" /
    "kept=['0']" errors these two datasets hit before their real shape
    was inspected directly.

    value_column: fixed column NAME shared across every sibling file (e.g.
    RSEM's "expected_count").

    value_column_index: 0-based column POSITION instead, for confirmed
    real cases where the value column's NAME varies per file (GSE139145,
    GSE312476: literally named after the sample; GSE116702/GSE198167/
    GSE294736: the only remaining column after the id, whatever pandas
    auto-names it). Position is counted AMONG COLUMNS AFTER id_column is
    removed (i.e. excludes the id column itself), consistently whether
    id_column was named or implicit.
    """
    glob_pattern = entry["glob_pattern"]
    value_column = entry.get("value_column")
    value_column_index = entry.get("value_column_index")
    if (value_column is None) == (value_column_index is None):
        raise ValueError(
            "per_sample_merge: set exactly one of value_column (fixed "
            "name) or value_column_index (0-based position among "
            "columns after the id column)."
        )
    id_column = entry.get("id_column")
    header = entry.get("header", True)
    rollup = entry.get("rollup")
    rollup_id_column = entry.get("rollup_id_column", id_column)
    sample_name_regex = entry.get("sample_name_regex")
    exclude_files_regex = entry.get("exclude_files_regex")
    symbol_column = entry.get("symbol_column")

    files = sorted(glob.glob(os.path.join(directory, glob_pattern)))
    if exclude_files_regex:
        # NOTE: glob_pattern only supports real glob syntax (*, ?, [...])
        # — Python's glob module has NO bash-style {a,b} brace expansion,
        # confirmed directly (glob.glob('*.{yaml,py}') matches nothing).
        # This is the supported way to exclude a subset of files that a
        # single positive glob_pattern can't cleanly express on its own
        # (e.g. GSE275243's confirmed-corrupt "-HOMO-" replicates, which
        # share the same GSM*_P3-*.txt.gz prefix as the real HET/WT files).
        pat = re.compile(exclude_files_regex)
        excluded = [f for f in files if pat.search(os.path.basename(f))]
        files = [f for f in files if not pat.search(os.path.basename(f))]
        if excluded:
            _log(f"per_sample_merge: exclude_files_regex={exclude_files_regex!r} dropped {len(excluded)} file(s): {[os.path.basename(f) for f in excluded]}")
    if not files:
        return None

    series_by_sample = {}
    id_to_symbol = {}
    for f in files:
        base = os.path.basename(f)
        if sample_name_regex:
            m = re.search(sample_name_regex, base)
            sample_name = m.group(1) if m else base
        else:
            sample_name = re.sub(r"\.gz$", "", base)

        if id_column is None:
            header_arg = 0 if header else None
            df = pd.read_csv(f, sep="\t", index_col=0, header=header_arg)
            this_value_column = df.columns[value_column_index] if value_column_index is not None else value_column
            s = df[this_value_column]
        else:
            # NOTE 2026-09-09: previously used a fixed header_arg (0 or
            # None), which breaks on files with a caption/blank line before
            # the real header (confirmed real for 10 of GSE312476's 11
            # sibling files — see _detect_header_row docstring). When a
            # named id_column is expected, scan for the real header row
            # per-file and skip up to it explicitly rather than assuming
            # every sibling file has an identical preamble.
            if header:
                header_row = _detect_header_row(f, id_column)
                df = pd.read_csv(f, sep="\t", skiprows=header_row, header=0)
            else:
                df = pd.read_csv(f, sep="\t", header=None)
            this_value_column = df.columns[value_column_index] if value_column_index is not None else value_column
            if rollup == "sum_by_gene":
                s = df.groupby(rollup_id_column)[this_value_column].sum()
            else:
                s = df.set_index(id_column)[this_value_column]
            # NOTE 2026-09-09: TF_activity.regulon() needs var_names in
            # gene-SYMBOL space (see _resolve_var_names docstring) — most
            # per_sample_merge sources are RSEM/Ensembl-keyed with NO
            # symbol column at all (those are deferred, not fixed here),
            # but confirmed real exceptions (GSE312476's gene_symbol,
            # GSE139145's GeneSymbol) do carry one alongside id_column.
            # Captured here, from `df` before set_index/groupby discards
            # it, keyed by whichever column ends up as the final merged
            # index (rollup_id_column when rolling up, else id_column).
            if symbol_column and symbol_column in df.columns:
                mapping_key_column = rollup_id_column if rollup == "sum_by_gene" else id_column
                for gene_id, sym in zip(df[mapping_key_column], df[symbol_column]):
                    id_to_symbol.setdefault(gene_id, sym)
        series_by_sample[sample_name] = s
        _log(f"per_sample_merge: {base} -> sample={sample_name!r} ({len(s)} genes, value_column={this_value_column!r}, id_column={id_column!r}, header={header})")

    mat = pd.DataFrame(series_by_sample).fillna(0)
    adata = ad.AnnData(X=sp.csr_matrix(mat.values.T))
    adata.obs_names = list(mat.columns)
    if symbol_column:
        missing = [gene_id for gene_id in mat.index if gene_id not in id_to_symbol]
        if missing:
            _log(f"per_sample_merge: symbol_column={symbol_column!r} — {len(missing)} of {len(mat.index)} ids had no symbol mapping (keeping original id for those)")
        adata.var_names = [str(id_to_symbol.get(gene_id, gene_id)) for gene_id in mat.index]
    else:
        adata.var_names = [str(g) for g in mat.index]
    _log(f"branch=per_sample_merge n_files={len(files)} symbol_column={symbol_column!r} shape={adata.shape}")
    return adata


# --------------------------------------------------------------------------
# format: multi_file_merge — N condition-level files, each contributing its
# own subset of real sample columns, merged column-wise. Confirmed real
# shapes: GSE242321 (Atf6-fpkm.txt.gz / WT-fpkm.txt.gz), GSE296324
# (KO_summary.txt.gz / wildtype_summary.txt.gz).
# --------------------------------------------------------------------------
def _load_multi_file_merge(directory, entry):
    # NOTE 2026-09-09: symbol_column is entry-level (applies uniformly to
    # every sibling file), not per-file-spec — confirmed real case
    # (GSE242321's Atf6-fpkm.txt.gz/WT-fpkm.txt.gz) shares one gene_symbol
    # column across both files. The concat below must stay indexed by the
    # file's own id column (join="inner" needs a shared key across files);
    # the id->symbol mapping is built alongside that from each file's own
    # symbol_column, then applied to the merged matrix's index at the end.
    files_spec = entry["files"]
    symbol_column = entry.get("symbol_column")
    frames = []
    id_to_symbol = {}
    for spec in files_spec:
        hit = _find_one(directory, [spec["pattern"]])
        if hit is None:
            raise ValueError(f"multi_file_merge: pattern {spec['pattern']!r} matched nothing in {directory}")
        sep = "," if hit.endswith((".csv", ".csv.gz")) else "\t"
        df = pd.read_csv(hit, sep=sep, index_col=0)
        if symbol_column:
            if symbol_column not in df.columns:
                raise ValueError(
                    f"multi_file_merge: symbol_column {symbol_column!r} not found in "
                    f"{os.path.basename(hit)}'s columns: {list(df.columns)}"
                )
            for gene_id, sym in zip(df.index, df[symbol_column]):
                id_to_symbol.setdefault(gene_id, sym)
        sub_entry = {k: v for k, v in spec.items() if k != "pattern"}
        sub = _select_columns(df, sub_entry, source_name=os.path.basename(hit))
        frames.append(sub)
    mat = pd.concat(frames, axis=1, join="inner")
    adata = ad.AnnData(X=sp.csr_matrix(mat.values.T))
    adata.obs_names = [str(c) for c in mat.columns]
    if symbol_column:
        missing = [gene_id for gene_id in mat.index if gene_id not in id_to_symbol]
        if missing:
            _log(f"multi_file_merge: symbol_column={symbol_column!r} — {len(missing)} of {len(mat.index)} ids had no symbol mapping (keeping original id for those)")
        adata.var_names = [str(id_to_symbol.get(gene_id, gene_id)) for gene_id in mat.index]
    else:
        adata.var_names = [str(g) for g in mat.index]
    _log(f"branch=multi_file_merge n_files={len(files_spec)} symbol_column={symbol_column!r} shape={adata.shape}")
    return adata


# --------------------------------------------------------------------------
# format: cellranger_h5 — confirmed real schema via h5py inspection
# (mm10plus2/{data,indices,indptr,shape,barcodes,gene_names,genes} for
# GSE120462, an older CellRanger v2-style genome-keyed export).
# --------------------------------------------------------------------------
def _load_cellranger_h5(directory, entry):
    import scanpy as sc

    hit = _find_one(directory, entry.get("glob_patterns", ["*_raw_gene_bc_matrices_h5.h5", "*.h5"]))
    if hit is None:
        return None
    genome = entry.get("h5_genome_key")
    adata = sc.read_10x_h5(hit, genome=genome) if genome else sc.read_10x_h5(hit)
    adata.var_names_make_unique()
    _log(f"branch=cellranger_h5 file={os.path.basename(hit)} genome_key={genome!r} raw_shape={adata.shape}")

    # NOTE 2026-09-09: the file GEO shipped for GSE120462 is CellRanger's
    # RAW (unfiltered) matrix — all ~737k possible barcodes, not just real
    # cells — because these are role="bulk" 10x runs, not single-cell
    # experiments. Every other bulk-format branch in this loader already
    # returns one column per sample; this branch was the one exception,
    # handing back per-barcode structure that TF_activity.regulon() then
    # tries to densify whole (confirmed: this is what raised the 76.9 GiB
    # MemoryError for GSE120462_GFP/Ikzf2). Collapse to a single pseudobulk
    # column here for role="bulk" entries, consistent with the rest of the
    # loader's bulk-format outputs.
    if entry.get("role") == "bulk":
        summed = np.asarray(adata.X.sum(axis=0)).ravel()
        sample_name = entry.get("sample_id") or entry.get("ref_id") or os.path.basename(hit)
        pseudobulk = ad.AnnData(X=sp.csr_matrix(summed.reshape(1, -1)))
        pseudobulk.obs_names = [str(sample_name)]
        pseudobulk.var_names = adata.var_names
        _log(
            f"branch=cellranger_h5 role=bulk collapsed {adata.shape[0]} barcodes -> "
            f"pseudobulk shape={pseudobulk.shape}"
        )
        return pseudobulk

    return adata


# --------------------------------------------------------------------------
# format: cellranger_h5_merge — N sibling per-sample CellRanger .h5 files,
# each collapsed to one pseudobulk column (same collapse logic as
# _load_cellranger_h5's role="bulk" path), merged column-wise into one
# multi-sample matrix. Confirmed real shape: GSE281207_Tmie_control (6
# GSM-per-sample .h5 files — P21-HET-A/B/C/F, P21-KO-B/E/F (4 HET + 3 KO,
# confirmed via actual load: n_files=7) — CellRanger
# 7.0.0, Single Cell 3' v3, single genome embedded directly in
# matrix/features/genome, confirmed via h5py inspection, so no
# h5_genome_key needed here unlike GSE120462's older v2-style export).
#
# This directory also contains filelist.txt (a GEO-generated manifest:
# Name/Time/Type/Size columns) and a raw .tar of the same h5 files —
# glob_pattern must be specific enough to exclude both, since
# _load_single_matrix_file's own glob list (tried earlier in
# load_rna_only's auto-detect chain for fmt=None) was previously matching
# filelist.txt via its bare "*.txt" pattern before ever reaching this
# branch. Auto-detection was NOT extended to include this branch — a
# directory holding N sibling single-cell .h5 files, mixed with the raw
# .tar/manifest, has no safe generic signature to auto-detect against, so
# `format: cellranger_h5_merge` must be set explicitly per entry.
# --------------------------------------------------------------------------
def _load_cellranger_h5_merge(directory, entry):
    import h5py
    import scanpy as sc

    pattern = entry.get("glob_pattern", "*.h5")
    hits = sorted(glob.glob(os.path.join(directory, pattern)))
    if not hits:
        return None

    # sample_name_source: "library_id" (default) reads each file's own
    # library_ids h5 attr — more robust than filename parsing since it's
    # baked into the CellRanger export itself (confirmed present:
    # library_ids=[b'P21-HET-A'] etc. for every GSE281207 file). Falls back
    # to "filename" for h5 files that don't carry a usable attr.
    name_source = entry.get("sample_name_source", "library_id")

    pseudobulk_cols = {}
    var_names_ref = None
    for hit in hits:
        adata = sc.read_10x_h5(hit)
        adata.var_names_make_unique()
        if var_names_ref is None:
            var_names_ref = list(adata.var_names)
        elif list(adata.var_names) != var_names_ref:
            raise ValueError(
                f"cellranger_h5_merge: {os.path.basename(hit)}'s gene set/order "
                f"differs from the first file's ({os.path.basename(hits[0])}) — "
                f"sibling h5 files must share identical features to merge as "
                f"pseudobulk columns of one matrix."
            )
        if name_source == "library_id":
            with h5py.File(hit, "r") as f:
                lib_ids = f.attrs.get("library_ids")
            sample_name = (
                lib_ids[0].decode() if lib_ids is not None and len(lib_ids)
                else os.path.splitext(os.path.basename(hit))[0]
            )
        else:
            sample_name = os.path.splitext(os.path.basename(hit))[0]
        summed = np.asarray(adata.X.sum(axis=0)).ravel()
        pseudobulk_cols[sample_name] = summed
        _log(
            f"cellranger_h5_merge: {os.path.basename(hit)} -> "
            f"sample={sample_name!r} ({adata.shape[0]} barcodes collapsed)"
        )

    mat = pd.DataFrame(pseudobulk_cols, index=var_names_ref)
    adata_out = ad.AnnData(X=sp.csr_matrix(mat.values.T))
    adata_out.obs_names = list(mat.columns)
    adata_out.var_names = [str(g) for g in mat.index]
    _log(f"branch=cellranger_h5_merge n_files={len(hits)} shape={adata_out.shape}")
    return adata_out


def load_rna_only(entry):
    """
    entry: one dict from extra_bulk_rna_inputs or model_construction_refs
           (already .format()-resolved by the Snakefile's _resolve()).
    Returns an AnnData with .var['gene_ids'] and .obs['barcode'] set, which
    is exactly what TF_activity.regulon() expects (see LingerGRN source).
    """
    directory = entry.get("path") or os.path.dirname(entry.get("rna_prefix", ""))
    prefix = entry.get("rna_prefix")
    fmt = entry.get("format")
    name = entry.get("sample_id") or entry.get("ref_id")

    adata = None
    if fmt == "per_sample_merge":
        adata = _load_per_sample_merge(directory, entry)
    elif fmt == "multi_file_merge":
        adata = _load_multi_file_merge(directory, entry)
    elif fmt == "cellranger_h5":
        adata = _load_cellranger_h5(directory, entry)
    elif fmt == "cellranger_h5_merge":
        adata = _load_cellranger_h5_merge(directory, entry)
    elif fmt == "excel":
        adata = _load_excel(directory, entry)
    else:
        if fmt == "mtx" or fmt is None:
            adata = _load_mtx_trio(directory, prefix=prefix)
        if adata is None and (fmt == "matrix" or fmt is None):
            adata = _load_single_matrix_file(directory, entry=entry)
        if adata is None and (fmt == "h5ad" or fmt is None):
            adata = _load_h5ad(directory)
        if adata is None and fmt is None:
            adata = _load_excel(directory, entry)

    if adata is None:
        raise ValueError(
            f"[{name}] Could not auto-detect RNA layout under {directory!r}. "
            f"Inspect the downloaded files and add an explicit `format:` "
            f"override (mtx | matrix | h5ad | excel | per_sample_merge | "
            f"multi_file_merge | cellranger_h5) to this entry in "
            f"config_eddie.yaml."
        )

    adata.var_names_make_unique()
    adata.var["gene_ids"] = adata.var_names
    adata.obs["barcode"] = adata.obs_names
    _log(f"FINAL: {name} shape={adata.shape}")
    return adata
