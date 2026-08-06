"""
requantify_consensus.py — Module 2c, ONE sample per invocation.

Ported from the original 02c_requantify_consensus.py. Confirmed against
installed snapatac2==2.9.0 (snap.pp.import_data doesn't exist in this
version — the correct function is snap.pp.import_fragments). Uses the
`whitelist` parameter to import only barcodes that already passed Module 1's
QC, rather than importing every raw barcode and subsetting afterward.
"""
import anndata as ad

try:
    import snapatac2 as snap
except ImportError as e:
    raise ImportError(
        "snapatac2 not installed in this env — it's pinned in envs/linger_preproc.yaml"
    ) from e


def load_chrom_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


chrom_sizes = load_chrom_sizes(snakemake.input.chrom_sizes)
sample_id = snakemake.wildcards.sample

qc_barcodes = sorted(set(ad.read_h5ad(snakemake.input.atac_filtered).obs_names))
print(f"{sample_id}: re-quantifying {len(qc_barcodes)} QC-passed barcodes against consensus peaks")

snap_data = snap.pp.import_fragments(
    snakemake.input.fragments,
    chrom_sizes=chrom_sizes,
    is_paired=True,
    file=None,                 # in-memory AnnData; whitelist keeps this small
    min_num_fragments=0,       # QC already applied in Module 1 — don't re-filter here
    sorted_by_barcode=False,   # unknown/unverified sort order — safer to assume unsorted
    whitelist=qc_barcodes,     # only import cells that already passed Module 1 QC
)

peak_matrix = snap.pp.make_peak_matrix(snap_data, peak_file=snakemake.input.consensus_bed)
peak_matrix.write_h5ad(snakemake.output.atac_consensus)

print(f"  -> {peak_matrix.shape[0]} cells x {peak_matrix.shape[1]} consensus peaks "
      f"-> {snakemake.output.atac_consensus}")
