# SETUP.md — getting the pipeline running on Eddie

Everything below assumes your existing `setup_linger.sh` and
`download_datasets.sh` runs are the ones referenced (same paths, same env
names). This guide is what to do *in addition to* those two scripts to get
the new Snakemake pipeline (Modules 1-3) actually executing.

## 1. Directory layout

Three separate trees, matching what `setup_linger.sh`/`download_datasets.sh`
already assume — nothing to invent here, just confirming the target layout:

```
/exports/eddie/scratch/s2906787/
├── cochlear_datasets/                  # raw GEO downloads (download_datasets.sh)
│   ├── 01_multiome/
│   └── 02_chromatin_priors/
├── linger_pipeline/                    # setup_linger.sh's WORK_DIR
│   ├── conda_envs/, conda_pkgs/        # conda env storage (redirected off $HOME)
│   ├── references/mm10/                # mm10.fa.gz, GRCm38.102.gtf.gz, mm10.chrom.sizes
│   ├── LINGER_repo/
│   ├── code/                           # <-- THIS Snakemake repo goes here
│   └── preprocessed/                   # <-- Snakemake pipeline outputs land here
└── ...
```

`code/` and `preprocessed/` don't exist yet — you're creating them now.

## 2. Get the repo onto Eddie

From your local machine, unzip `linger_pipeline_snakemake.zip` and:

```bash
scp -r linger_pipeline s2906787@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s2906787/linger_pipeline/code
```

(Same `scp` pattern your existing `run_preprocessing_job.sh` docstring
already documents for copying code over.)

## 3. Two environments, not one

This pipeline needs two *separate* conda envs:

- **An orchestrator env** (e.g. `snakemake_eddie`) — runs `snakemake` itself.
  Modern Snakemake (8+) requires **Python ≥3.11**, so this can't be
  `linger_preproc` (which `setup_linger.sh` deliberately keeps on
  Python 3.10 to match the scanpy/snapatac2 stack it verified there).
- **A pipeline env built from `envs/linger_preproc.yaml`** — this is what
  every rule's `conda:` directive actually references. `--use-conda` builds
  it fresh the first time any rule needs it; you don't create this one by
  hand.

**The one thing that has to be true before any of this works:** the `conda`
binary active in your orchestrator env must be **≥24.7.1** — Snakemake
enforces this itself whenever it has to build a conda env from a yaml spec,
regardless of which pipeline is asking. If you're on Eddie's module-loaded
conda (commonly older than this), fix it once, scoped to your own env only:

```bash
conda activate snakemake_eddie      # or whatever you named it, python>=3.11
conda install -c conda-forge "conda>=24.7.1"
conda --version                     # confirm >=24.7.1
pip install snakemake snakemake-executor-plugin-sge
```

This only changes `snakemake_eddie`'s own private `conda` binary (via
`PATH` precedence while that env is active) — it doesn't touch the
system module, `linger_preproc`, or anything else.

**Build the pipeline's own conda env once, ahead of time**, so the actual
run doesn't need outbound network access mid-job (useful if compute nodes
are more network-restricted than the login node):

```bash
cd /exports/eddie/scratch/s2906787/linger_pipeline/code
conda activate snakemake_eddie
snakemake --use-conda --conda-frontend conda --conda-create-envs-only --cores 1
```

This builds one env (hashed from `envs/linger_preproc.yaml`'s contents)
under `.snakemake/conda/` inside `code/`, shared by every rule across
Modules 1-3. `--use-conda` is required on every subsequent invocation too
— without it, Snakemake ignores `conda:` directives entirely and runs
rules in whatever env launched `snakemake` (which doesn't have
scanpy/snapatac2 — that fails immediately on the first script rule).

## 4. Point the config at your real paths

Open `code/config/config_eddie.yaml` and check/edit:

| Key | Should match |
|---|---|
| `multiome_root` | `.../cochlear_datasets/01_multiome` (download_datasets.sh's `${MULTI}`) |
| `chrom_root` | `.../cochlear_datasets/02_chromatin_priors` |
| `scratch` | where you want Module 1-3 **outputs** written — recommend `.../linger_pipeline/preprocessed` |
| `references.mm10_chrom_sizes` | `.../linger_pipeline/references/mm10/mm10.chrom.sizes` (setup_linger.sh downloads this automatically in Step 11 — confirm it exists before running Module 2) |

The `samples:` / `extra_bulk_rna_inputs:` / `extra_chromatin_prior_inputs:`
blocks already mirror your `sample_config.py` exactly — no edits needed
there unless `SAMPLES` membership itself changes (e.g. a corrected GEO
deposit shows up).

## 5. Set up the Eddie profile (one-time, per-user)

Nothing to copy — `profiles/eddie/config.yaml` ships inside the repo and is
used directly by path. Just confirm your username isn't hardcoded anywhere
it shouldn't be (it isn't currently — all paths flow from `config_eddie.yaml`).

## 6. Dry run

Always do this before submitting anything:

```bash
cd /exports/eddie/scratch/s2906787/linger_pipeline/code
conda activate snakemake_eddie
snakemake -n --cores 4 --use-conda --conda-frontend conda
```

Check the job count and target list look right (21 jobs for 4 samples
through Module 3, as verified in the build). If a rule errors here, it's a
path/config problem — fix it before it ever touches the scheduler.

## 7. Run it

In a `tmux` session (this will run for hours):

```bash
tmux new-session -s linger_run
cd /exports/eddie/scratch/s2906787/linger_pipeline/code
conda activate snakemake_eddie    # orchestrator env — NOT linger_preproc, see section 3
snakemake --profile profiles/eddie --use-conda --conda-frontend conda --rerun-incomplete --keep-going \
    > /exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/snakemake_run.log 2>&1 &
echo "Snakemake PID: $!"
```

Detach (`Ctrl-b d`), monitor with `qstat` and:

```bash
tail -f /exports/eddie/scratch/s2906787/linger_pipeline/preprocessed/snakemake_run.log
```

Per-rule logs land under `<scratch>/logs/`, exactly as declared in each
`.smk` rule's `log:` block.

## 8. What to check when it finishes (Modules 1-2)

- `<scratch>/module1_summary.csv` — per-sample cell/gene/peak counts
- `<scratch>/qc_sanity_checks/sanity_summary.csv` — flags anything needing
  a look before Module 3 (barcode mismatches, peak count drift, doublet
  rate outliers)
- `<scratch>/qc_sanity_checks/*.png` — QC grids + doublet score overlay

If `sanity_summary.csv` flags nothing, Module 3 already ran too (it's in
`rule all`) — check `<scratch>/module3_integration/` next.

## 9. The one manual step: cell type labels

Module 3's clustering/UMAP/marker-scoring runs automatically, but the actual
IHC/OHC/supporting-cell **call** doesn't, by design:

1. Open `<scratch>/module3_integration/cluster_umap.png` and
   `cluster_marker_scores.csv` (per-cluster mean expression of
   `Myo7a`/`Pou4f3`/`Gfi1`/`Slc26a5`/`Sox2`/`Hes1`).
2. Copy `<scratch>/module3_integration/cluster_template.tsv` to the path
   given in `config["integration"]["cluster_annotation_tsv"]` (same
   directory, `cluster_annotation.tsv`, by default).
3. Fill in the blank `cell_type` column for every `leiden` row.
4. Run:
   ```bash
   snakemake --profile profiles/eddie --use-conda <scratch>/module3_integration/labeled.h5ad
   ```

That's the full Phase 1-3 pipeline through to a labeled, integrated
AnnData ready for Module 4 (LINGER initialization) once that's built.

## Pending manual actions carried over from setup_linger.sh

These aren't part of the Snakemake DAG and still need doing separately if
you haven't already:
1. `qsub <WORK_DIR>/homer_mm10_install.sh` — **now a real hard prerequisite**,
   not just pending housekeeping: Module 4's `linger_motif_scan` rule calls
   HOMER's `annotatePeaks.pl` directly. Confirm this job actually completed
   (check `perl configureHomer.pl -list` shows `mm10` as `installed`) before
   running Module 4.
2. `qsub <WORK_DIR>/star_index_job.sh` (only if aligning bulk RNA-seq FASTQs — Module 7)
3. Confirm `mm10.chrom.sizes` downloaded successfully (Step 11 of setup_linger.sh) — Module 2's `requantify_consensus` rule will fail immediately without it.

## 10. Module 4/6 — using your EXISTING LINGER setup (corrected 2026-07-24, verified 2026-07-24)

Earlier guidance here assumed `GRNdir` and the `LINGER` conda env still
needed building. Your `setup_linger.sh` shows both are already done, and
you've since confirmed the remaining open questions on real Eddie output —
**Module 4 is now ready to run end-to-end, nothing left to verify.**

**Fixed against your real setup:**

| What | Was (wrong) | Now (confirmed) |
|---|---|---|
| `config["linger"]["grn_dir"]` | `references/LINGER_provide_data/` (guessed, didn't exist) | `provide_data/provide_data` |
| Module 4/6 conda env | rebuilt from `envs/linger.yaml` every run | reuses existing `LINGER` env directly by path (see below) |
| HOMER | bundled as a conda package | called at its real standalone install path |
| HOMER mm10 genome package | assumed pending | **confirmed installed** (`configureHomer.pl -list` shows `+  mm10  v7.0`) |
| `linger_motif_scan`'s motif file | guessed `all_motif_rmdup.motif` (didn't exist) | **confirmed** `all_motif_rmdup_Mammal` — GRNdir splits by taxon (`_fly`/`_Mammal`/`_Plant`); mouse uses `_Mammal` |
| `linger_weights/` | unknown, flagged for you to check | **confirmed empty** (0 files) — not used by anything, safe to ignore |
| `references/jaspar/` | unclear if needed | **confirmed unused** — nothing in LINGER's `scNN` code path reads raw `.jaspar` files; `GRNdir` already has everything needed (`Match_TF_motif_Mus_musculus.txt`, `genome_map_homer.txt`, `TSS_mm10.txt`, `all_motif_rmdup_Mammal`) |

**Why the conda env fix matters:** your `LINGER` env has exact pins
(`scanpy==1.9.5`, `anndata==0.9.2`, `scipy==1.11.3`, `LingerGRN==1.105
--no-deps`, `rpy2`) that `setup_linger.sh`'s own comments say were arrived
at after a failed attempt with looser versions broke `import LingerGRN`.
Module 4/6 rules now point `conda:` at
`config["linger"]["conda_env_linger"]` — an **absolute path to the env
directory itself**, not a `.yaml` spec. Verified this is real, supported
Snakemake behavior by reading `deployment/conda.py` directly (not assumed):
a `conda:` value that's an existing local directory is parsed as
`CondaEnvDirSpec`, which sets `Env.is_externally_managed = True` — Snakemake
activates that directory directly via `conda activate <path>` and never
tries to solve or recreate it. `envs/linger.yaml` is kept in the repo as
documentation only; no rule uses it.

**Run Module 4:**
```bash
snakemake --profile profiles/eddie \
    <scratch>/module4_linger_init/tss_redist.done \
    <scratch>/module4_linger_init/MotifTarget.bed
```
Note: `--use-conda` alone is enough for Module 4/6 — you don't need
`--conda-frontend conda`/`--conda-create-envs-only` the way Modules 1-3
do, since the externally-managed env path means Snakemake activates your
existing env rather than building one.

Module 6 (`module6_all`) becomes reachable once Module 3's `labeled.h5ad`
exists (Section 9 above) *and* Module 4 has completed.
