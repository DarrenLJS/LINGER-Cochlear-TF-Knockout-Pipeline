"""extract_peaks.py — writes one BED file for one sample's atac_filtered.h5ad.
Ported unchanged in logic from the original 02a_extract_peaks.py; loop over
SAMPLES replaced by Snakemake wildcards (one job per sample)."""
import anndata as ad


def peakname_to_bed_row(name):
    # expects "chr:start-end"
    chrom, rest = name.split(":")
    start, end = rest.split("-")
    return chrom, start, end


atac = ad.read_h5ad(snakemake.input.atac)
n_written, n_skipped = 0, 0
with open(snakemake.output.bed, "w") as f:
    for name in atac.var_names:
        try:
            chrom, start, end = peakname_to_bed_row(name)
            f.write(f"{chrom}\t{start}\t{end}\n")
            n_written += 1
        except ValueError:
            n_skipped += 1

print(f"{snakemake.wildcards.sample}: wrote {n_written} peaks "
      f"({n_skipped} unparseable names skipped)")
