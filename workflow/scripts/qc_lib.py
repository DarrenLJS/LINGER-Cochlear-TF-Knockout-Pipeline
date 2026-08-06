"""
qc_lib.py — shared loader + QC/doublet-removal functions for Module 1.

Ported unchanged in logic from the original standalone 01_qc_preprocessing.py.
The only structural change vs. the original: there, these functions were
called from a `for sample in SAMPLES` loop inside one long-running script with
per-sample try/except so one crash didn't lose progress on the rest. Under
Snakemake, each sample is its own wildcard-based job — failure isolation and
"don't reprocess what already succeeded" are handled by the DAG + output-file
staleness tracking instead, so that loop/try-except scaffolding is gone; the
per-sample logic itself (loaders, QC filters, Scrublet handling) is untouched.
"""
import gzip

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.io as sio
import scipy.sparse as sp

sc.settings.verbosity = 1


# =============================================================================
# FORMAT-DETECTING FILE HELPERS
# Some GEO submitters name a file "*.gz" without ever actually running gzip on
# it (GSE157398's ATAC barcodes/matrix/peaks files are like this — real gzip
# always starts with bytes \x1f\x8b; these start with plain ASCII). Never trust
# the extension — check the actual magic bytes and open accordingly. BGZF
# (block-gzip, used for tabix-indexed files like fragments.tsv.gz) does start
# with \x1f\x8b and is gzip-compatible, so it's handled by the same gzip path.
# =============================================================================

def is_gzipped(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def smart_open_text(path):
    return gzip.open(path, "rt") if is_gzipped(path) else open(path, "rt")


def smart_mmread(path):
    if is_gzipped(path):
        with gzip.open(path, "rb") as f:
            return sio.mmread(f)
    with open(path, "rb") as f:
        return sio.mmread(f)


# =============================================================================
# LOADERS — one per file-layout type in the sample manifest
# =============================================================================

def load_h5_combined(sample):
    """10x .h5 with both Gene Expression + Peaks feature_types in one matrix."""
    adata = sc.read_10x_h5(sample["h5_path"], gex_only=False)
    adata.var_names_make_unique()

    ft_col = "feature_types" if "feature_types" in adata.var.columns else None
    if ft_col is None:
        raise ValueError(
            f"[{sample['sample_id']}] .h5 has no feature_types column — "
            "cannot split RNA vs Peaks. Inspect adata.var manually."
        )

    rna = adata[:, adata.var[ft_col] == "Gene Expression"].copy()
    atac = adata[:, adata.var[ft_col] == "Peaks"].copy()

    if atac.n_vars == 0:
        raise ValueError(
            f"[{sample['sample_id']}] No 'Peaks' features found in combined .h5 — "
            "this file may be RNA-only. Check features.tsv feature_types values."
        )
    return rna, atac


def load_mtx_combined(sample):
    """10x mtx trio (barcodes/features/matrix) with combined GEX+Peaks features."""
    prefix = sample["combined_prefix"]
    adata = sc.read_mtx(prefix + "matrix.mtx.gz").T
    barcodes = pd.read_csv(prefix + "barcodes.tsv.gz", header=None, sep="\t")[0].values
    features = pd.read_csv(prefix + "features.tsv.gz", header=None, sep="\t")

    adata.obs_names = barcodes
    adata.var_names = features[1].values if features.shape[1] > 1 else features[0].values
    adata.var_names_make_unique()

    if features.shape[1] >= 3:
        adata.var["feature_types"] = features[2].values
    else:
        raise ValueError(
            f"[{sample['sample_id']}] features.tsv has no 3rd (feature_types) column — "
            "cannot split RNA vs Peaks automatically."
        )

    rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    atac = adata[:, adata.var["feature_types"] == "Peaks"].copy()
    if atac.n_vars == 0:
        raise ValueError(f"[{sample['sample_id']}] No Peaks features found in combined mtx.")
    return rna, atac


def load_mtx_separate(sample):
    """Fully separate RNA mtx trio and ATAC mtx trio + peaks.bed.gz (GSE157398 layout).

    GSE157398's RNA files (GSM4764088/89) are genuinely gzipped. Its ATAC files
    (GSM4764090/91) are a mix: barcodes/matrix/peaks are plain text wearing a
    ".gz" name, while fragments.tsv.gz is real BGZF. smart_open_text/smart_mmread
    handle both cases per-file rather than assuming from the extension.
    """
    rna_prefix = sample["rna_prefix"]
    rna_mat = smart_mmread(rna_prefix + "matrix.mtx.gz").T.tocsr()
    rna = ad.AnnData(X=rna_mat)
    with smart_open_text(rna_prefix + "barcodes.tsv.gz") as f:
        rna.obs_names = [line.strip() for line in f if line.strip()]
    with smart_open_text(rna_prefix + "features.tsv.gz") as f:
        rna_feature_rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    rna.var_names = [row[1] if len(row) > 1 else row[0] for row in rna_feature_rows]
    rna.var_names_make_unique()

    atac_prefix = sample["atac_prefix"]
    atac_mat = smart_mmread(atac_prefix + "matrix.mtx.gz").T.tocsr()
    atac = ad.AnnData(X=atac_mat)
    with smart_open_text(atac_prefix + "barcodes.tsv.gz") as f:
        atac.obs_names = [line.strip() for line in f if line.strip()]

    with smart_open_text(sample["atac_peaks_bed"]) as f:
        peak_lines = [line.strip().split("\t")[:3] for line in f if line.strip()]
    peak_names = [f"{c}:{s}-{e}" for c, s, e in peak_lines]
    if len(peak_names) != atac.n_vars:
        raise ValueError(
            f"[{sample['sample_id']}] peaks.bed row count ({len(peak_names)}) != "
            f"ATAC matrix column count ({atac.n_vars}) — matrix may be peaks x barcodes "
            "instead of barcodes x peaks, or the bed file doesn't match this matrix."
        )
    atac.var_names = peak_names
    atac.var_names_make_unique()
    return rna, atac


LOADERS = {
    "h5_combined": load_h5_combined,
    "mtx_combined": load_mtx_combined,
    "mtx_separate": load_mtx_separate,
}


# =============================================================================
# QC PIPELINE
# =============================================================================

def qc_filter_rna(rna, params):
    rna.var["mt"] = rna.var_names.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(rna, min_genes=params["min_genes_per_cell"])
    sc.pp.filter_genes(rna, min_cells=params["min_cells_per_gene"])
    rna = rna[rna.obs["pct_counts_mt"] < params["max_pct_mito"]].copy()
    return rna


def run_scrublet(rna, params, sample=None):
    random_state = params.get("scrublet_random_state", 0)
    manual_threshold = (sample or {}).get("scrublet_manual_threshold")
    try:
        sc.pp.scrublet(
            rna,
            expected_doublet_rate=params["scrublet_expected_doublet_rate"],
            random_state=random_state,
        )

        if manual_threshold is not None:
            # Per-sample override — set in config_eddie.yaml for samples where
            # automatic thresholding has proven unreliable run-to-run (e.g.
            # GSE182202_P8: 38.6% doublets on one run, ~0% on another rerun of
            # the identical input, no code/param changes in between).
            rna.obs["predicted_doublet"] = rna.obs["doublet_score"] > manual_threshold
            n_doublets = int(rna.obs["predicted_doublet"].sum())
            threshold_source = "manual"
            print(f"    Using manual threshold={manual_threshold:.3f} from config — "
                  f"flagging {n_doublets}/{rna.n_obs} cells as doublets")
        else:
            n_doublets = int(rna.obs["predicted_doublet"].sum())
            threshold_source = "automatic"
            print(f"    Scrublet flagged {n_doublets}/{rna.n_obs} cells as doublets (random_state={random_state})")

            # Automatic threshold can land above the entire score distribution,
            # returning 0 doublets even when the histogram shows a real upper
            # tail. Fall back to a manual threshold at the valley between the
            # singlet peak and the upper tail when that happens.
            if n_doublets == 0:
                threshold = _find_doublet_threshold(rna.obs["doublet_score"].values)
                if threshold is not None:
                    rna.obs["predicted_doublet"] = rna.obs["doublet_score"] > threshold
                    n_doublets = int(rna.obs["predicted_doublet"].sum())
                    threshold_source = "fallback"
                    print(f"    Automatic threshold found 0 doublets — "
                          f"using fallback threshold={threshold:.3f}, now flagging {n_doublets}/{rna.n_obs}")
                else:
                    print("    WARNING: automatic threshold found 0 doublets and no fallback "
                          "threshold could be determined — leaving as 0, inspect doublet_scores.png manually")

        # Doublet removal below drops predicted_doublet==True cells, so
        # rna.obs["predicted_doublet"] is trivially all-False in the SAVED
        # file regardless of how many were actually removed. Stash the real
        # numbers in .uns BEFORE filtering so sanity_checks.py can report
        # what actually happened, not a tautology computed after the fact.
        n_pre = rna.n_obs
        rna.uns["n_doublets_removed"] = n_doublets
        rna.uns["n_cells_pre_doublet_removal"] = n_pre
        rna.uns["doublet_rate_observed"] = (n_doublets / n_pre) if n_pre > 0 else 0.0
        rna.uns["doublet_threshold_source"] = threshold_source

        rna = rna[~rna.obs["predicted_doublet"]].copy()
    except Exception as e:
        print(f"    WARNING: Scrublet failed ({e}) — skipping doublet removal for this sample")
    return rna


def _find_doublet_threshold(scores, min_threshold=0.08, max_threshold=0.5):
    """Valley-finding fallback: locate the local minimum between the main
    singlet peak and the upper doublet tail via a smoothed KDE, restricted
    to a plausible range so it can't land inside the main peak itself."""
    from scipy.stats import gaussian_kde
    if scores.max() < min_threshold:
        return None
    kde = gaussian_kde(scores)
    xs = np.linspace(min_threshold, min(max_threshold, scores.max()), 400)
    density = kde(xs)
    minima = [i for i in range(1, len(xs) - 1) if density[i] < density[i - 1] and density[i] < density[i + 1]]
    if not minima:
        return None
    return float(xs[minima[0]])


def qc_filter_atac(atac, params):
    atac.obs["n_peaks"] = (atac.X > 0).sum(axis=1).A1 if sp.issparse(atac.X) else (atac.X > 0).sum(axis=1)
    atac = atac[atac.obs["n_peaks"] >= params["min_peaks_per_cell"]].copy()
    peak_ncells = (atac.X > 0).sum(axis=0).A1 if sp.issparse(atac.X) else (atac.X > 0).sum(axis=0)
    atac = atac[:, peak_ncells >= params["min_cells_per_peak"]].copy()
    return atac
