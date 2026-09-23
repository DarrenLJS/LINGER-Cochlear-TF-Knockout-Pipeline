# Setup guide — running the pipeline on Eddie

This guide takes a new user from an empty Eddie scratch space to a working
pipeline: environments, reference data, dataset downloads, configuration, and
the first runs through Module 4. Once Module 4 has completed, continue with
the **Run order** table in [README.md](README.md#usage).

**Conventions used below**

- `<user>` — your Eddie username (e.g. `s1234567`).
- `<scratch>` — the pipeline output root set as `scratch` in
  `config/config_eddie.yaml` (recommended: `/exports/eddie/scratch/<user>/linger_pipeline/preprocessed`).
- `<WORK_DIR>` — `setup_linger.sh`'s working directory
  (`/exports/eddie/scratch/<user>/linger_pipeline`).

**Before you start:** allow 150–500 GB of scratch space for the raw
downloads, 12–48 hours for the downloads to finish, and several hours for
environment setup. Both setup scripts are idempotent and safe to re-run.

---

## Contents

1. [Directory layout](#1-directory-layout)
2. [Run the setup scripts](#2-run-the-setup-scripts)
3. [Finish the manual installs](#3-finish-the-manual-installs)
4. [Get the repository onto Eddie](#4-get-the-repository-onto-eddie)
5. [Create the orchestrator environment](#5-create-the-orchestrator-environment)
6. [Build the Modules 1–3 environment](#6-build-the-modules-13-environment)
7. [Point the config at your paths](#7-point-the-config-at-your-paths)
8. [Check the Eddie profile](#8-check-the-eddie-profile)
9. [Dry run](#9-dry-run)
10. [Run Modules 1–3](#10-run-modules-13)
11. [Check the Module 1–2 output](#11-check-the-module-12-output)
12. [Label cell types (manual step)](#12-label-cell-types-manual-step)
13. [LINGER modules: Module 4 and onward](#13-linger-modules-module-4-and-onward)

---

## 1. Directory layout

Three separate trees, matching what `setup_linger.sh` and
`download_datasets.sh` assume:

```
/exports/eddie/scratch/<user>/
├── cochlear_datasets/                  # raw GEO downloads (download_datasets.sh)
│   ├── 01_multiome/
│   └── 02_chromatin_priors/
├── linger_pipeline/                    # setup_linger.sh's WORK_DIR
│   ├── conda_envs/, conda_pkgs/        # conda env storage (redirected off $HOME)
│   ├── references/mm10/                # mm10.fa.gz, GRCm38.102.gtf.gz, mm10.chrom.sizes
│   ├── provide_data/provide_data/      # LINGER reference data (GRNdir)
│   ├── homer/                          # standalone HOMER install
│   ├── LINGER_repo/
│   ├── code/                           # <-- this Snakemake repository goes here
│   └── preprocessed/                   # <-- pipeline outputs land here (<scratch>)
└── ...
```

`code/` and `preprocessed/` are created in the steps below.

## 2. Run the setup scripts

Both scripts are in `setup_scripts/`. Copy them to your scratch space first
(or run them from a local copy of the repository). They are independent, so
run them in parallel in two `tmux` sessions.

**Session 1 — environments and references (`setup_linger.sh`).** Run on an
interactive node:

```bash
ssh <user>@eddie.ecdf.ed.ac.uk
tmux new-session -s setup
qlogin -l h_vmem=8G -l h_rt=04:00:00
module load anaconda
bash /exports/eddie/scratch/<user>/setup_linger.sh \
     /exports/eddie/scratch/<user>/linger_pipeline
```

This creates the `linger_preproc` and `LINGER` conda environments, downloads
the mm10 references (including `mm10.chrom.sizes`), LINGER `provide_data`,
the HOMER core install and the LINGER repository, and writes the job scripts
used in step 3. Set `INSTALL_R=1` before the command to also install R and
Seurat (about 60 minutes extra). Watch progress with:

```bash
tail -f /exports/eddie/scratch/<user>/linger_pipeline/setup.log
```

**Session 2 — dataset downloads (`download_datasets.sh`).** Runs on a login
node; no `qlogin` or `module load` needed (it only uses `wget`):

```bash
tmux new-session -s download
bash /exports/eddie/scratch/<user>/download_datasets.sh \
     /exports/eddie/scratch/<user>/cochlear_datasets
```

All `.tar` bundles are extracted automatically. Interrupted downloads resume
on re-run. Watch progress with:

```bash
tail -f /exports/eddie/scratch/<user>/cochlear_datasets/download.log
```

## 3. Finish the manual installs

These steps sit outside the Snakemake DAG.

1. **HOMER mm10 genome package — required before Module 4.** Module 4's
   `linger_motif_scan` rule calls HOMER's `annotatePeaks.pl` directly.
   ```bash
   qsub <WORK_DIR>/homer_mm10_install.sh
   # once the job has finished:
   perl <WORK_DIR>/homer/configureHomer.pl -list | grep mm10   # expect: +  mm10  v7.0
   ```
2. **Confirm `mm10.chrom.sizes` exists** (downloaded in Step 11 of
   `setup_linger.sh`). Module 2's `requantify_consensus` rule fails
   immediately without it:
   ```bash
   ls -l <WORK_DIR>/references/mm10/mm10.chrom.sizes
   ```
3. **Patch the installed LingerGRN package — required before Module 6.** Two
   bugs in LingerGRN's single-cell-type `scNN` code are fixed by the patch
   scripts in `setup_scripts/`. Run each once with the `LINGER` env's Python
   (both are idempotent and keep a backup of the original file):
   ```bash
   <path-to-LINGER-env>/bin/python setup_scripts/patch_LL_net_RE_ordering.py
   <path-to-LINGER-env>/bin/python setup_scripts/patch_LL_net_cis_reg_load.py
   ```
   Find `<path-to-LINGER-env>` with `conda env list | grep LINGER`.
4. **STAR index — optional.** `qsub <WORK_DIR>/star_index_job.sh` is only
   needed if you plan to align bulk RNA-seq FASTQs yourself; Module 7 uses
   pre-quantified matrices by default.

## 4. Get the repository onto Eddie

From your local machine, copy the repository into `code/`:

```bash
scp -r LINGER-Cochlear-TF-Knockout-Pipeline-master \
    <user>@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/<user>/linger_pipeline/code
```

## 5. Create the orchestrator environment

The pipeline uses three separate conda environments:

- **Orchestrator (`snakemake_eddie`)** — runs `snakemake` itself. You create
  this one by hand (below). Modern Snakemake (8+) requires **Python ≥ 3.11**,
  so it cannot be `linger_preproc`, which is kept on Python 3.10 to match its
  scanpy/snapatac2 stack.
- **Modules 1–3 (`linger_preproc`, from `envs/linger_preproc.yaml`)** — what
  the Modules 1–3 rules' `conda:` directives reference. Snakemake builds it
  (step 6); you do not create it by hand.
- **LINGER (`LINGER`)** — built by `setup_linger.sh` and used by Module 4
  onward (step 13).

**The conda binary in the orchestrator env must be ≥ 24.7.1.** Snakemake
enforces this whenever it builds a conda env from a YAML spec. Eddie's
module-loaded conda is often older, so install a newer one into this env
only:

```bash
conda create -n snakemake_eddie python=3.11
conda activate snakemake_eddie
conda install -c conda-forge "conda>=24.7.1"
conda --version                     # confirm >= 24.7.1
pip install snakemake snakemake-executor-plugin-sge
```

This only changes `snakemake_eddie`'s own `conda` binary (through `PATH`
precedence while the env is active); it does not touch the system module,
`linger_preproc` or anything else.

> **Do not install anything else into `snakemake_eddie`** (for example
> `scanpy`). It exists only to run the `snakemake` CLI, and a solver-driven
> install can silently shift versions Snakemake depends on. Use
> `linger_preproc` for ad-hoc analysis of Modules 1–3 output (reading
> `.h5ad` files, checking markers), since that is the env that produced it.

## 6. Build the Modules 1–3 environment

Build the pipeline's own conda env once, ahead of the real run, so jobs do
not need outbound network access mid-run (compute nodes can be more
network-restricted than login nodes):

```bash
cd /exports/eddie/scratch/<user>/linger_pipeline/code
conda activate snakemake_eddie
snakemake --use-conda --conda-frontend conda --conda-create-envs-only --cores 1
```

This builds one env (named by a hash of `envs/linger_preproc.yaml`) under
`.snakemake/conda/` inside `code/`, shared by every Modules 1–3 rule.

`--use-conda` is required on **every** later invocation too. Without it,
Snakemake ignores `conda:` directives and runs rules in whatever env
launched `snakemake`, which has no scanpy/snapatac2 and fails on the first
script rule.

## 7. Point the config at your paths

Open `code/config/config_eddie.yaml` and check or edit:

| Key | Should match |
|---|---|
| `multiome_root` | `.../cochlear_datasets/01_multiome` (`download_datasets.sh`'s `${MULTI}`) |
| `chrom_root` | `.../cochlear_datasets/02_chromatin_priors` |
| `scratch` | Where outputs are written — recommended `.../linger_pipeline/preprocessed` |
| `references.mm10_chrom_sizes` | `.../linger_pipeline/references/mm10/mm10.chrom.sizes` |
| `linger.grn_dir` | `provide_data/provide_data` |
| `linger.conda_env_linger` | Absolute path to the `LINGER` env directory |

The `samples`, `extra_bulk_rna_inputs` and `extra_chromatin_prior_inputs`
blocks already describe the full dataset collection; no edits are needed
unless the sample set changes (for example, a corrected GEO deposit).

> **Eddie conda path:** `module load anaconda` loads the group-managed
> install first, which can take precedence over a personal `envs_dirs`
> setting. The `LINGER` env may therefore end up under the group path
> (`/exports/csce/eddie/biology/groups/.../envs/<user>/LINGER`) rather than
> under `<WORK_DIR>/conda_envs/`. Run `conda env list | grep LINGER` and use
> the path it shows for `linger.conda_env_linger`.

## 8. Check the Eddie profile

Nothing to copy: `profiles/eddie/config.yaml` ships with the repository and
is used in place (`--profile profiles/eddie`). No usernames are hard-coded in
it; all paths come from `config_eddie.yaml`.

## 9. Dry run

Always dry-run before submitting anything:

```bash
cd /exports/eddie/scratch/<user>/linger_pipeline/code
conda activate snakemake_eddie
snakemake --profile profiles/eddie --use-conda --conda-frontend conda \
    --rerun-incomplete --keep-going --latency-wait 60 -n
```

The dry run uses the same profile and flags as the real run, so it checks
the same configuration.

Check that the job count and target list look right (about 21 jobs for four
samples through Module 3). An error at this stage is a path or config
problem; fix it before anything reaches the scheduler.

## 10. Run Modules 1–3

`rule all` covers Modules 1–3 and the QC sanity checks. Run it in `tmux`,
since it takes several hours:

```bash
tmux new-session -s linger_run
cd /exports/eddie/scratch/<user>/linger_pipeline/code
conda activate snakemake_eddie    # the orchestrator env, not linger_preproc
snakemake --profile profiles/eddie --use-conda --conda-frontend conda \
    --rerun-incomplete --keep-going --latency-wait 60 \
    > /exports/eddie/scratch/<user>/linger_pipeline/preprocessed/snakemake_run.log 2>&1 &
echo "Snakemake PID: $!"
```

Detach with `Ctrl-b d`. Monitor with `qstat` and:

```bash
tail -f /exports/eddie/scratch/<user>/linger_pipeline/preprocessed/snakemake_run.log
```

When `rule all` finishes, build the per-sample QC summary. It is not part of
`rule all`, so it needs its own command:

```bash
snakemake --profile profiles/eddie --use-conda --conda-frontend conda \
    --rerun-incomplete --keep-going --latency-wait 60 \
    /exports/eddie/scratch/<user>/linger_pipeline/preprocessed/module1_summary.csv
```

Per-rule logs are written to `<scratch>/logs/`, as declared in each rule's
`log:` block.

## 11. Check the Module 1–2 output

- `<scratch>/module1_summary.csv` — per-sample cell, gene and peak counts.
- `<scratch>/qc_sanity_checks/sanity_summary.csv` — flags anything to review
  before Module 3 (barcode mismatches, peak-count drift, doublet-rate
  outliers).
- `<scratch>/qc_sanity_checks/*.png` — QC grids and doublet-score overlay.

If nothing is flagged, check `<scratch>/module3_integration/` next; Module 3's
automatic steps run as part of `rule all`.

## 12. Label cell types (manual step)

Module 3's clustering, UMAP and marker scoring run automatically, but the
cell-type call (hair cell, supporting cell, etc.) is manual by design:

1. Open `<scratch>/module3_integration/cluster_umap.png` and
   `cluster_marker_scores.csv` (per-cluster mean expression of `Myo7a`,
   `Pou4f3`, `Gfi1`, `Slc26a5`, `Sox2` and `Hes1`).
2. Copy `<scratch>/module3_integration/cluster_template.tsv` to the path set
   in `config["integration"]["cluster_annotation_tsv"]`
   (`cluster_annotation.tsv` in the same directory by default).
3. Fill in the `cell_type` column for every `leiden` row. Avoid commas,
   slashes and parentheses in labels; they break SGE array-job wildcards.
4. Build the labelled AnnData:
   ```bash
   snakemake --profile profiles/eddie --use-conda --conda-frontend conda \
       --rerun-incomplete --keep-going --latency-wait 60 \
       <scratch>/module3_integration/labeled.h5ad
   ```

This produces a labelled, integrated AnnData ready for Module 4.

## 13. LINGER modules: Module 4 and onward

Modules 4 and 6–11 run in the existing `LINGER` env built by
`setup_linger.sh`, not in a Snakemake-built env.

**Expected configuration.** These values have been confirmed against a
working setup:

| Setting | Value | Notes |
|---|---|---|
| `linger.grn_dir` | `provide_data/provide_data` | Contains `Match_TF_motif_Mus_musculus.txt`, `genome_map_homer.txt`, `TSS_mm10.txt`, `all_motif_rmdup_Mammal` |
| Conda env for Modules 4 and 6–11 | Existing `LINGER` env, by absolute path | Set via `linger.conda_env_linger`; never rebuilt from `envs/linger.yaml` |
| HOMER | Standalone install under `<WORK_DIR>/homer/` | Called at its install path, not as a conda package |
| HOMER mm10 genome package | Installed | `configureHomer.pl -list` shows `+  mm10  v7.0` |
| Motif file for `linger_motif_scan` | `all_motif_rmdup_Mammal` | GRNdir splits motifs by taxon (`_fly`, `_Mammal`, `_Plant`); mouse uses `_Mammal` |
| `linger_weights/` | Empty | Not used by any rule; safe to ignore |
| `references/jaspar/` | Unused | The `scNN` code path reads no raw `.jaspar` files; GRNdir has everything needed |

**Why the env is referenced by path.** The `LINGER` env has exact pins
(`scanpy==1.9.5`, `anndata==0.9.2`, `scipy==1.11.3`, `LingerGRN --no-deps`,
`rpy2`) chosen after looser versions broke `import LingerGRN`. The rules
point `conda:` at `config["linger"]["conda_env_linger"]`, an absolute path to
the env directory rather than a YAML spec. Snakemake parses an existing
local directory as `CondaEnvDirSpec` (`Env.is_externally_managed = True`)
and activates it directly with `conda activate <path>`, without solving or
recreating it. `envs/linger.yaml` is kept for documentation only.

**Run Module 4:**

```bash
snakemake --profile profiles/eddie --use-conda --conda-frontend conda \
    --rerun-incomplete --keep-going --latency-wait 60 \
    <scratch>/module4_linger_init/tss_redist.done \
    <scratch>/module4_linger_init/MotifTarget.bed
```

Because the `LINGER` env is externally managed, `--conda-frontend conda` has
no effect from Module 4 onward; it is kept so every command uses the same
flags.

Module 6 (`module6_all`) becomes runnable once `labeled.h5ad` exists
(step 12), Module 4 has completed, and the LingerGRN patches are applied
(step 3). From there, follow the **Run order** and **Manual checkpoints** in
[README.md](README.md#usage), including the Module 8 and 8b sanity checks
that must be read before trusting any knockout output.
