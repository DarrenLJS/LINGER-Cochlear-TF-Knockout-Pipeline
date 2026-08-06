"""
sync_rna_atac_barcodes.py — Module 2d, ONE sample per invocation.

Fixes RNA/ATAC barcode desync after consensus re-quantification: 02c
re-quantifies ATAC against the consensus peak set using a whitelist of
Module 1's QC-passed barcodes, but cells with zero fragments over the new
(smaller) peak set get dropped from atac_consensus.h5ad. rna_filtered.h5ad
is never touched by 02c, so any sample that lost ATAC cells ends up with
RNA and ATAC on different cell sets.

Architectural change vs. the original 02d_sync_rna_atac_barcodes.py: that
script overwrote rna_filtered.h5ad in place (with a manual
.pre_sync_backup.h5ad copy made first, guarded by an idempotency check so
reruns wouldn't clobber the backup). Under Snakemake, a file must not be
both the output of one rule (Module 1's qc_preprocess_sample) and an
in-place target of a later rule — that breaks the DAG's staleness tracking
and is exactly the kind of thing the manual backup/idempotency dance was
working around by hand. This script instead writes a NEW file,
rna_synced.h5ad, leaving rna_filtered.h5ad as Module 1's untouched,
original output. No backup step needed: rna_filtered.h5ad already *is* the
backup, permanently, by construction.
"""
import sys
import scanpy as sc

sample_id = snakemake.wildcards.sample

rna = sc.read_h5ad(snakemake.input.rna_filtered)
atac = sc.read_h5ad(snakemake.input.atac_consensus)

rna_bc = set(rna.obs_names)
atac_bc = set(atac.obs_names)

if not atac_bc.issubset(rna_bc):
    extra = atac_bc - rna_bc
    sys.exit(
        f"ABORT {sample_id}: atac_consensus has {len(extra)} barcodes NOT in "
        "rna_filtered — this is not the expected attrition pattern (ATAC can only "
        "lose cells relative to the Module 1 whitelist, never gain new ones). "
        "Needs manual inspection before syncing."
    )

if list(rna.obs_names) == list(atac.obs_names):
    print(f"{sample_id}: already synced ({rna.n_obs} cells) — writing through unchanged")
    rna.write_h5ad(snakemake.output.rna_synced)
else:
    rna_synced = rna[atac.obs_names].copy()
    rna_synced.write_h5ad(snakemake.output.rna_synced)
    print(f"{sample_id}: synced {rna.n_obs} -> {rna_synced.n_obs} cells "
          "(now matches atac_consensus)")
