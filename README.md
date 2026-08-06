# LINGER Cochlear GRN Pipeline — Snakemake port (Modules 1-3)

Architecture mirrors the tRNA-seq Three-pass Alignment Snakemake Pipeline:
config-driven paths/resources, one `.smk` file per module, scripts under
`workflow/scripts/`, per-rule conda envs, an `sge_extra()` helper for Eddie/SGE.

## Status

**Built and dry-run verified** (DAG resolves cleanly, 21 jobs, all scripts
byte-compile):
- Module 1 — QC & preprocessing (ported from `01_qc_preprocessing.py`)
- Module 2 — consensus peak set, re-quantification, RNA/ATAC barcode sync,
  sanity checks (ported from `02a`-`02d` + `sanity_checks.py`)
- Module 3 — Harmony integration, Leiden clustering, UMAP, marker scoring
  (new; cell-type label *assignment* stays manual by design)

**Built 2026-07-24, verified against real Eddie output the same day** —
Module 4 (LINGER initialisation) and Module 6 (GRN inference), against the
real `LingerGRN==1.105` API actually installed in your `LINGER` conda env
(confirmed via source diff against the 1.110 sdist I first inspected — every
function used has an identical signature in both). All Module 4 prerequisites
confirmed on real Eddie output: `GRNdir` path, HOMER mm10 genome package
installed, motif filename (`all_motif_rmdup_Mammal`), and the pre-built
`LINGER` conda env reused by direct path rather than rebuilt. See SETUP.md
§10 for the full verification trail. Module 4 is ready to run.

Key findings baked into the rules/scripts, from inspecting the real source:

- LINGER has two training modes. The full atlas-scale-pretrained `'LINGER'`
  method is hardcoded to hg19/hg38 throughout (`hg19_hg38_pair.bed`,
  `RE_gene_corr_hg19.bed`, etc.) — **not usable for mouse data**. `method='scNN'`
  is the only mode with native mouse support (confirmed against
  `docs/scNN.md`, which explicitly uses `genome='mm10'`) — this is what
  every Module 4/6 rule uses.
- `LINGER_tr.get_TSS()` and `RE_TG_dis()` read/write hardcoded relative
  `./data/...` paths, ignoring their `outdir` argument. The scripts below
  `os.chdir()` into a shared per-run workdir to work around this rather than
  patching LINGER itself.
- `pseudo_bulk.pseudo_bulk()` requires cell-type labels **before**
  pseudobulking — Module 4 depends on `labeled.h5ad` (Module 3's manual
  output), not `integrated.h5ad`.
- `LL_net.TF_RE_binding(method='scNN')` needs a HOMER-generated
  `MotifTarget.bed` — an external tool step, not a LingerGRN function.

**Not yet built** — Module 5 (chromatin priors), Modules 7-10 (bulk TF
activity, perturbation, validation, benchmarking). Module 5 in particular
has an open design question — see below.

## Module 5 open question — Hi-C/CUT&RUN aren't native LINGER inputs

The HTML plan describes Module 5 as converting Hi-C loops and CUT&RUN peaks
into "structured constraints LINGER's GRN inference can use directly." That
isn't accurate to the installed package: `LingerGRN==1.110`'s public API
(`preprocess.py`, `LL_net.py`, `LINGER_tr.py`) has no ingestion path for
Hi-C contact maps or CUT&RUN peaks as inference-time priors — GRN inference
only consumes motif-based TF-RE binding (HOMER) and RE-TG genomic distance.

Three ways to actually build Module 5, none of which is "just call a
function":
1. **Benchmarking-only** (cheapest): treat Hi-C/CUT&RUN purely as Module 10
   ground truth, as the HTML's own Module 10 description already implies.
   Drop Module 5 as a distinct inference-time step.
2. **Post-hoc filtering**: use CUT&RUN peaks to filter/reweight Module 6's
   inferred TF-RE edges after the fact (a light post-processing script, not
   a LINGER-internal change).
3. **True manifold-regularization injection**: modify LINGER's own training
   loop (`LINGER_tr.py`'s `sc_nn_NN`/`training`) to accept Hi-C/CUT&RUN as an
   additional regularization term — a real research contribution, not a
   pipeline-engineering task, and probably out of scope for this project's
   timeline.

Not building Module 5 until this is decided — see SETUP.md.

## What changed vs. the original standalone scripts

1. **Per-sample loop → per-sample Snakemake job.** Each script's internal
   `for sample in SAMPLES: try/except` became one wildcard-based rule
   per sample. Failure isolation and resume-without-reprocessing now come
   from `--keep-going --rerun-incomplete` + Snakemake's own output-based
   staleness tracking, not manual `if os.path.exists(...)` checks.
2. **02d no longer overwrites `rna_filtered.h5ad` in place.** It writes a
   new file, `rna_synced.h5ad`. A file being both a rule's declared output
   and another rule's in-place mutation target breaks Snakemake's DAG
   assumptions — this was also the reason the original needed a manual
   `.pre_sync_backup.h5ad` + idempotency guard, which is no longer needed:
   `rna_filtered.h5ad` is permanently untouched, by construction.
3. **`EXPECTED_N_PEAKS` in sanity checks is now computed at run time** from
   `consensus_peaks.bed`'s line count, not hardcoded. The old constant had
   already gone stale once (235228 → 224719 after `GSE182202_P1` was
   excluded from `SAMPLES`) and needed a manual code edit every time
   `SAMPLES` membership changed.
4. **`download_datasets.sh` and `setup_linger.sh` stay outside the DAG**,
   run manually once — same convention as the reference pipeline's
   `00_reference_prep.smk` PREREQUISITES block (manual download, not an
   automated rule).

## Usage

```bash
# Dry run (check the DAG)
snakemake -n --cores 4

# Local run
snakemake --use-conda --cores 8 --rerun-incomplete --keep-going

# Eddie / SGE — profile ships in-repo, no copying required
snakemake --profile profiles/eddie
```

`--profile` accepts any directory containing a `config.yaml`, so the Eddie
profile lives at `profiles/eddie/config.yaml` inside the repo itself rather
than needing to be copied to `~/.config/snakemake/eddie/` (unlike the
reference tRNA-seq pipeline, where `eddie_profile_config.yaml` sits at the
repo root as a file you copy over manually). Point `--profile` at an
absolute path to `profiles/eddie` if you're invoking Snakemake from
somewhere other than the repo root.

### Module 3's manual step

`integrate_harmony` + `score_markers` are in `rule all` and run
automatically. `apply_cell_type_labels` is **not** — after reviewing
`module3_integration/cluster_umap.png` and `cluster_marker_scores.csv`:

1. Copy `module3_integration/cluster_template.tsv` to the path in
   `config["integration"]["cluster_annotation_tsv"]`.
2. Fill in the `cell_type` column for every `leiden` cluster.
3. Run: `snakemake --use-conda --cores 4 <scratch>/module3_integration/labeled.h5ad`

Before editing `config/config_eddie.yaml`'s paths for your own run, update
`multiome_root`, `chrom_root`, `scratch`, and `references.mm10_chrom_sizes`
to match your Eddie scratch layout.
