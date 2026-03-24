# nf-multicaller-exomecnv

A Nextflow pipeline for comprehensive genomic analysis combining multiple variant callers with exome and copy number variation (CNV) detection. This workflow is based on a workflow by **Dr Phelelani T Mpangase** ([phelelani](https://github.com/phelelani)): [nf-exomecnv](https://github.com/phelelani/nf-exomecnv), which incorporated CANOES, CLAMMS, XHMM and InDelible, but stopped at the actual variant calling component, without including structural variant merging and *in silico* validation with a ML classifier. The ```modules-icav2-dragen.nf``` module is based off a script by **Regan Cannell** ([RWCannell](https://github.com/RWCannell)): [cram_input_dragen_ica_workflow](https://github.com/SBIMB/ica-elwazi/tree/main/nextflow_workflows/cram_input_dragen_ica_workflow). The ```modules-cnvkit.nf``` module follows the instructions laid out in the [CNVkit Copy Number Calling Pipeline](https://cnvkit.readthedocs.io/en/stable/pipeline.html). The ```modules-gatk-gcnv.nf``` module follows the instructions laid out in the [GATK Guide on How to Call Rare Germline Copy Number Variants](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants). The ```modules-truvari.nf``` module follows instructions laid out in the [Truvari Wiki](https://github.com/ACEnglish/truvari/wiki). The ```modules-survivor.nf``` module follows instructions laid out in the [SURVIVOR Wiki](https://github.com/fritzsedlazeck/SURVIVOR/wiki). 

## Tools

**CANOES**: [GitHub](https://github.com/ShenLab/CANOES); [DOI](https://doi.org/10.1093/nar/gku345)

**CLAMMS**: [GitHub](https://github.com/rgcgithub/clamms); [DOI](https://doi.org/10.1093/bioinformatics/btv547)

**XHMM**: [GitHub](https://github.com/RRafiee/XHMM); [DOI](https://doi.org/10.1002/0471142905.hg0723s81)

**InDelible**: [GitHub](https://github.com/HurlesGroupSanger/indelible); [DOI](https://doi.org/10.1016/j.ajhg.2021.09.010)

**ICAv2**: [GitHub](https://github.com/umccr/illumination/tree/v2); [ICAv2 Commands](https://help.ica.illumina.com/command-line-interface/cli-indexcommands)

**Illumina's DRAGEN Germline Enrichment**: [Description of DRAGEN Enrichment](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/dragen-enrichment.html); [DOI](https://doi.org/10.1038/s41587-024-02382-1)

**The Broad Institute's GATK-gCNV**: [GitHub](https://github.com/broadinstitute/gatk); [DOI](https://doi.org/10.1038/s41588-023-01449-0)

**CNVkit**: [GitHub](https://github.com/etal/cnvkit); [DOI](https://doi.org/10.1371/journal.pcbi.1004873)

**Truvari**: [GitHub](https://github.com/ACEnglish/truvari); [DOI](https://doi.org/10.1186/s13059-022-02840-6)

**SURVIVOR**: [GitHub](https://github.com/fritzsedlazeck/SURVIVOR); [DOI](https://doi.org/10.1038/ncomms14061)

## Overview

This pipeline implements a robust bioinformatics workflow for analyzing genomic data using:
- Multiple variant calling approaches
- Exome analysis
- Copy number variation (CNV) detection

Built with Nextflow for reproducible, scalable, and portable data analysis.

## Workflow Steps

The pipeline is modular. Each workflow can be run independently, and detailed
methodology/process documentation is now split into module-specific pages.

### Module documentation

| Workflow | Module docs |
|---|---|
| `canoes` | [docs/modules/canoes.md](docs/modules/canoes.md) |
| `xhmm` | [docs/modules/xhmm.md](docs/modules/xhmm.md) |
| `clamms` | [docs/modules/clamms.md](docs/modules/clamms.md) |
| `gcnv` | [docs/modules/gatk-gcnv.md](docs/modules/gatk-gcnv.md) |
| `cnvkit` | [docs/modules/cnvkit.md](docs/modules/cnvkit.md) |
| `dragen` | [docs/modules/dragen.md](docs/modules/dragen.md) |
| `indelible` | [docs/modules/indelible.md](docs/modules/indelible.md) |
| `survivor` | [docs/modules/survivor.md](docs/modules/survivor.md) |
| `truvari` | [docs/modules/truvari.md](docs/modules/truvari.md) |
| `feature_extraction` | [docs/modules/feature-extraction.md](docs/modules/feature-extraction.md) |
| `train` | [docs/modules/train.md](docs/modules/train.md) |
| `evaluate` | [docs/modules/evaluate.md](docs/modules/evaluate.md) |
| `normalise` | [docs/modules/normalise.md](docs/modules/normalise.md) |
| Shared caller post-processing | [docs/modules/common-postprocessing.md](docs/modules/common-postprocessing.md) |

## Features

- **Multi-caller approach**: Integrates multiple variant calling algorithms for robust variant detection
- **Exome analysis**: Specialized processing for whole exome sequencing (WES) data
- **CNV detection**: Comprehensive copy number variation analysis
- **Nextflow-based**: Ensures reproducibility and scalability across different computing environments
- **Containerized**: Includes Docker/Singularity support for dependency management

## Requirements

### Operating System
- **Linux** or **Unix** (macOS, HPC clusters running Linux)
  > Windows is **not** supported. Use WSL2 or a remote Linux server if you are on Windows.

### Software
- **Apptainer** (formerly Singularity): >= 1.0  
  All tool containers are pulled and executed via Apptainer. Install it following the [official guide](https://apptainer.org/docs/admin/main/installation.html). ([DOI](https://doi.org/10.1371/journal.pone.0177459))
- **Nextflow**: >= 20.04  ([DOI](https://doi.org/10.1038/nbt.3820)) 
- Install with:
  ```bash
  curl -s https://get.nextflow.io | bash
  ```
- **Java**: >= 8 (required by Nextflow)  
  Install OpenJDK:
  - Debian/Ubuntu: `sudo apt install default-jdk`
  - RHEL/CentOS/Rocky: `sudo yum install java-1.8.0-openjdk`
  - macOS (Homebrew): `brew install openjdk`

### Libraries and Dependencies
All bioinformatics tools (CANOES, CLAMMS, XHMM, GATK, CNVkit, InDelible, ICAv2 CLI, SURVIVOR, Truvari, bcftools, bedtools, Picard) are provided as pre-built container images and are downloaded automatically by Apptainer at runtime. No manual tool installation is required beyond Apptainer and Nextflow themselves.

### Python dependencies for tests
The pipeline itself runs in containers, but the local Python test suite (`tests/*.py`) requires these Python packages in a fresh environment:

- `pytest`
- `numpy`
- `pandas`
- `scikit-learn`
- `xgboost`
- `imbalanced-learn`
- `pysam`

Install them with:
```bash
python -m pip install pytest numpy pandas scikit-learn xgboost imbalanced-learn pysam
```

### Computing Resources
- Adequate CPU cores (parallel execution recommended)
- Sufficient disk space for intermediate and output files
- Memory requirements depend on dataset size

## Installation

1. Clone this repository:
```bash
git clone https://github.com/Yonatan-Ariel-Wolberg/nf-multicaller-exomecnv.git
cd nf-multicaller-exomecnv
```

2. Build the ICAv2 container image (Docker or Apptainer)
```bash
# Docker build (from the new Dockerfile)
docker build -f bin/Dockerfile.icav2 -t nf-multicaller-exomecnv/icav2-cli:2.43.0 .

# Optional: build Apptainer/Singularity image from the existing recipe
apptainer build icav2_cli_v2.43.0.sif bin/icav2.def
```

3. (Optional) Build the training container image
```bash
docker build -f bin/Dockerfile.train -t nf-multicaller-exomecnv/train:2.1.4 .

# Optional: build Apptainer/Singularity image from the existing recipe
apptainer build train_xgboost_deps.sif bin/train.def
```

## Usage

Each workflow is run by providing `--workflow <name>` along with a params file. Edit the corresponding `params/*.json` file to point to your actual input files before running.

### Individual CNV callers

The individual callers should be run first. Their VCF outputs feed into the consensus modules (SURVIVOR / Truvari).

```bash
# CANOES — exome CNV caller using read-depth, NNLS reference weighting and Negative Binomial HMM
nextflow run main.nf --workflow canoes -params-file params/general/params-canoes.json

# XHMM — exome CNV caller using GATK depth-of-coverage and PCA
nextflow run main.nf --workflow xhmm -params-file params/general/params-xhmm.json

# CLAMMS — exome CNV caller using nearest-neighbour normalisation
nextflow run main.nf --workflow clamms -params-file params/general/params-clamms.json

# GATK-gCNV — cohort-mode germline CNV caller from the Broad Institute
nextflow run main.nf --workflow gcnv -params-file params/general/params-gatk-gcnv.json

# CNVkit — coverage-based CNV caller with pooled reference normalisation
nextflow run main.nf --workflow cnvkit -params-file params/general/params-cnvkit.json

# DRAGEN — Illumina DRAGEN germline enrichment pipeline via ICAv2
nextflow run main.nf --workflow dragen -params-file params/general/params-icav2-dragen.json

# InDelible — small insertion/deletion caller from split reads
nextflow run main.nf --workflow indelible -params-file params/general/params-indelible.json
```

### Consensus modules (run after individual callers)

Point the `*_dir` paths in the params file to the VCF output directories produced by the callers above, then run:

```bash
# SURVIVOR — merge per-caller VCFs into a consensus call set
nextflow run main.nf --workflow survivor -params-file params/general/params-survivor.json

# Truvari — collapse per-caller VCFs into a consensus call set
nextflow run main.nf --workflow truvari -params-file params/general/params-truvari.json
```

`survivor` / `truvari` only produce the consensus VCF outputs.

### Full end-to-end workflow (BAM/CRAM → XGBoost model)

Use `--workflow full` to run available caller modules from BAM/CRAM inputs,
merge caller outputs, extract features, and train the XGBoost model in one run.

```bash
nextflow run main.nf --workflow full -params-file params/general/params-full.json
```

Notes:
- The full workflow automatically uses whichever caller inputs are configured in
  `params/general/params-full.json`.
- At least **two** caller inputs must be configured so consensus merging can run.
- Set `merger_mode` to `survivor` (default) or `truvari`.
- `truth_labels` is required for the final `train` step.

### Combined consensus + feature extraction workflows

```bash
nextflow run main.nf --workflow survivor_with_features -params-file params/general/params-survivor-with-features.json
nextflow run main.nf --workflow truvari_with_features -params-file params/general/params-truvari-with-features.json
```

`*_with_features` workflows run the same consensus merge step as `survivor`/`truvari`,
then immediately run `feature_extraction` on the merged outputs. They do **not** run `train`.

### Parameter files

All required and optional parameters for each workflow are documented in the corresponding `params/*.json` template:

| Params file | Workflow |
|---|---|
| `params/general/params-canoes.json` | `--workflow canoes` |
| `params/general/params-xhmm.json` | `--workflow xhmm` |
| `params/general/params-clamms.json` | `--workflow clamms` |
| `params/general/params-gatk-gcnv.json` | `--workflow gcnv` |
| `params/general/params-cnvkit.json` | `--workflow cnvkit` |
| `params/general/params-icav2-dragen.json` | `--workflow dragen` |
| `params/general/params-indelible.json` | `--workflow indelible` |
| `params/general/params-survivor.json` | `--workflow survivor` |
| `params/general/params-truvari.json` | `--workflow truvari` |
| `params/general/params-full.json` | `--workflow full` |

### Truth-label TSV requirements

For `--workflow train`, the `truth_labels` input TSV must include the following columns:

- `sample_id`
- `chrom`
- `start`
- `end`
- `cnv_type`
- `truth_label`

Example:
- `truth_label = 1` for true CNV calls.
- `truth_label = 0` for non-CNV calls.

### Sexinfo file requirements

For `--workflow clamms`, the `sexinfo` input file must be a headerless, tab-separated file with exactly two columns:

- Column 1: `sample_id` – used to look up sex for each sample in the samplesheet.
- Column 2: `sex` – biological sex of the sample; must be `M` (male) or `F` (female).

Example:

```
SAMPLE001	M
SAMPLE002	F
SAMPLE003	M
```

Every sample in the cohort must have an entry in this file; CLAMMS will abort with
an error if a sample ID is missing. The file does not need to contain exactly the
same set of IDs as `samplesheet.tsv` (extra IDs are ignored), but it must include
all sample IDs present in the samplesheet.

### Known unhandled use-case errors (current gaps)

The pipeline already validates many required parameters and empty-channel inputs.
The following use-case errors are **not yet explicitly pre-validated** and may
instead fail later inside a process/tool step:

- `--workflow gcnv`: `samples_path` glob matches files with unsupported
  extensions or missing companion index/metadata expected by downstream steps.
- `--workflow cnvkit`: `bams` are provided but one or more BAM index files are
  missing or stale.
- `--workflow train`: `truth_labels` file exists and has required columns but
  contains bad coordinate/value types or inconsistent labels.
- `--workflow feature_extraction`: merged VCFs pass basic VCF header checks but
  are still missing INFO/FORMAT annotations required for specific features.
- Any workflow: some workflow-specific secondary assets may exist but still be
  unreadable or invalid in ways not fully pre-validated.

## Running on the Wits UI Cluster

The workflow has been tested on the **University of the Witwatersrand (Wits) UI HPC cluster**, which runs **Linux** and uses **SLURM** as its job scheduler.  Use the `wits` profile to activate the cluster-specific settings.

### Quick start on the Wits cluster

```bash
# Load the required modules (adjust module names to match the cluster's module system)
module load nextflow
module load apptainer   # or: module load singularity (apptainer is already loaded by default on the Wits UI Cluster)

# Run a workflow — replace <workflow> and the params file as needed
nextflow run main.nf \
    -profile wits,medium \
    --workflow canoes \
    -params-file params/general/params-canoes.json
```

### What the `wits` profile does

| Setting | Value |
|---|---|
| `executor` | `slurm` |
| `process.queue` | `batch` |
| `process.cpus` | 4 |
| `process.memory` | 16 GB |
| `process.time` | 72 h |
| `executor.queueSize` | 500 |
| `executor.submitRateLimit` | 10/sec |

### Cluster-specific bind mounts

The global `singularity.runOptions` in `nextflow.config` binds the paths listed
in `params.bind_paths` into every container. For the Wits DDD datasets, use:

```
-B /dataB/aux
-B /home/ywolberg
-B /dataG/ddd
-B /dataG/ddd-2023
```

Relevant Wits input paths used by the `params/*-wits.json` templates:

- CANOES / CLAMMS / XHMM DDD-AFRICA samplesheet:
  `/home/ywolberg/DECIPHERING_DD_DATA/DDD_AFRICA_DATA/batch_3/samplesheet.tsv`
- CNVkit / GATK-gCNV DDD-AFRICA BAMs:
  `/home/ywolberg/DECIPHERING_DD_DATA/DDD_AFRICA_DATA/batch_3/organized_data/**/*.{bam,bam.bai}`
- INDELIBLE DDD-AFRICA family directories:
  `/home/ywolberg/DECIPHERING_DD_DATA/DDD_AFRICA_DATA/batch_3/organized_data/{Extended,Father,Mother,Proband}`
- DRAGEN upload roots:
  `/home/ywolberg/DECIPHERING_DD_DATA/{DDD_UK_DATA/bams,DDD_UK_DATA/crams,DDD_AFRICA_DATA/batch_3/organized_data/Proband}/**/*.{bam,bam.bai,cram,cram.crai}`

Region-specific examples are also provided for every Wits module template:

- DDD-AFRICA examples: `params/ddd-africa/params-*-wits-ddd-africa.json`
- DDD-UK examples: `params/ddd-uk/params-*-wits-ddd-uk.json`

The DDD files in `/home/ywolberg/DECIPHERING_DD_DATA/...` are symbolic links to
`/dataG/ddd` and `/dataG/ddd-2023`, so both `/dataG` locations must be
bind-mounted as well.

Wits templates (except ICAv2-DRAGEN reference handling) use the shared
`/dataG/ddd/data/resources` assets:

- Reference FASTA:
  `/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa`
- Reference FAI:
  `/dataG/ddd/data/resources/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai`
- Targets BED:
  `/dataG/ddd/data/resources/canoes/probes_sanger.bed`
- Targets interval list:
  `/dataG/ddd/data/resources/canoes/probes_sanger.interval_list`

Edit `params.runOptions` in `nextflow.config` (or pass `--runOptions` on the
command line) if your data lives elsewhere on the cluster.

### GATK / Picard constraint

GATK and Picard processes request nodes with AVX2 support (`--constraint=avx2`) and exclude node `n04` (`--exclude=n04`).  Adjust `process.clusterOptions` in `nextflow.config` if the Wits cluster topology changes.

### ICAv2 / DRAGEN on the Wits cluster

The DRAGEN workflow submits runs to the ICAv2 cloud platform from a local container.  Before running `--workflow dragen`:

1. Build the ICAv2 CLI container (if not already done):
   ```bash
   docker build -f bin/Dockerfile.icav2 -t nf-multicaller-exomecnv/icav2-cli:2.43.0 .
   # or:
   apptainer build icav2_cli_v2.43.0.sif bin/icav2.def
   ```
2. Update the container path in `nextflow.config` (`withLabel: 'icav2-dragen'`) to the location where you stored the `.sif` file.
3. Ensure your ICA credentials are present at `~/.icav2/` on the cluster head node so they can be bind-mounted into the container.

## Scalability: running 50, 300, or 2000 samples

The pipeline ships with three cohort-size profiles that pre-configure the
scheduler, concurrency, and per-caller batch sizes for small, medium, and large
runs.  Combine a cohort-size profile with the site profile using `-profile`:

```bash
# ≤ 50 samples
nextflow run main.nf -profile wits,small  -params-file params/general/params-canoes.json

# 50–300 samples (default behaviour if no cohort-size profile is specified)
nextflow run main.nf -profile wits,medium -params-file params/general/params-canoes.json

# 300–2000 samples
nextflow run main.nf -profile wits,large  -params-file params/general/params-canoes.json
```

### What each profile controls

| Setting | `small` (≤50) | `medium` (50–300) | `large` (300–2000) |
|---|---|---|---|
| `executor.queueSize` | 100 | 500 | 2000 |
| `executor.submitRateLimit` | 5/sec | 10/sec | 20/sec |
| `process.maxForks` | 25 | 100 | 200 |
| `canoes_batch_size` | 50 | 100 | 200 |
| `xhmm_batch_size` | 25 | 50 | 100 |
| GATK memory | 16 GB | 16 GB | 32 GB |
| `large_mem` label | 64 GB / 8 CPUs | 64 GB / 8 CPUs | 128 GB / 16 CPUs |
| `gpu_or_high_cpu` label | 32 GB / 16 CPUs | 32 GB / 16 CPUs | 64 GB / 32 CPUs |

### Caller-specific notes

#### CANOES (`--workflow canoes`)
CANOES uses `bedtools multicov` to count reads across all BAMs simultaneously.
Running a single `multicov` job over thousands of BAMs at once exhausts both
memory and file-descriptor limits.  The pipeline therefore splits the BAM list
into fixed-size batches (`canoes_batch_size`) and runs one `multicov` job per
batch per chromosome; the per-batch matrices are then merged into a single
cohort-wide matrix **before** the CANOES R script normalises read counts,
selects reference samples via NNLS, and calls CNVs using an HMM in which the
Negative Binomial distribution is the integrated emission model (the NB
log-probabilities are fed directly into the Viterbi decoder — they are a
single combined model, not separate steps).  The batch size therefore controls
**only parallelism and
memory use**, not normalisation accuracy — every sample contributes to the
final normalisation regardless of batch size.

- **`canoes_batch_size`** (default `100`; override in params file or via profile):
  - Up to 50 samples → set to `50` (1 batch; minimal overhead).
  - Up to 300 samples → set to `100` (3 batches per chromosome).
  - Up to 2000 samples → set to `200` (10 batches per chromosome).

#### XHMM (`--workflow xhmm`)
XHMM uses GATK `DepthOfCoverage`, which can require significant memory per job
when given a large BAM list.  The BAM list is split into fixed-size groups
(`xhmm_batch_size`) and one `GATK_DOC` job is submitted per group.  All groups
are merged by `xhmm --mergeGATKdepths` before PCA normalisation, so the batch
size controls **only job memory**, not the normalisation cohort.

- **`xhmm_batch_size`** (default `50`; override in params file or via profile):
  - Up to 50 samples → set to `25`–`50` (1–2 batches).
  - Up to 300 samples → set to `50` (up to 6 batches).
  - Up to 2000 samples → set to `100` (up to 20 batches; keeps memory per job manageable).

#### CLAMMS (`--workflow clamms`)
All CLAMMS steps (depth-of-coverage, normalisation, model training, CNV calling)
are fully per-sample and run in parallel automatically via Nextflow's data-flow
model.  No batch-size tuning is needed; `process.maxForks` limits the number of
concurrent per-sample jobs.  The custom reference panel is built from **all**
normalised coverage files in the cohort, ensuring accurate nearest-neighbour
normalisation regardless of cohort size.

#### GATK-gCNV (`--workflow gcnv`)
gCNV trains a cohort-level model and is naturally designed for large cohorts.
Interval parallelism is controlled by `scatter_count` in
`params/general/params-gatk-gcnv.json` (default `5000`).  For the cohort-level steps
(`FilterIntervals`, `DetermineGermlineContigPloidy`, `GermlineCNVCaller`) that
collect all samples at once, increase the memory/CPU limits using the `large`
profile when running hundreds or thousands of samples.

#### CNVkit (`--workflow cnvkit`)
CNVkit creates a pooled reference from **all** sample coverages before calling
CNVs, so every sample contributes to bias correction and GC normalisation.
Concurrency is governed by `process.maxForks`.  Use `test_size` (integer) or
`test_list` (comma-separated sample IDs) in the params file to restrict a run
to a subset of samples during development.

#### DRAGEN / ICAv2 (`--workflow dragen`)
DRAGEN runs are submitted to the ICAv2 platform.  `maxUploadForks` (default `4`)
limits simultaneous CRAM uploads.  Increase this if your network bandwidth and
ICA account limits allow it.

#### SURVIVOR / Truvari (`--workflow survivor` / `--workflow truvari`)
These consensus modules are fast and lightweight regardless of cohort size; no
additional tuning is required.

### Minimum cohort sizes for accurate normalisation

Each CNV caller's normalisation method requires a minimum number of samples to
produce statistically reliable results.  Running with fewer samples than the
recommended minimum will not cause the pipeline to fail, but CNV call quality
will be reduced.

| Caller | Normalisation method | Recommended minimum | Notes |
|--------|---------------------|--------------------:|-------|
| **CANOES** | Mean-coverage normalisation + NNLS reference weighting + Negative Binomial HMM | 30 | NNLS reference-panel quality degrades with very small cohorts; 30+ samples recommended for reliable reference weighting |
| **XHMM** | Mean-centering + PCA | 30 | XHMM's own filters may exclude all targets/samples with very small cohorts; 50+ strongly recommended |
| **CLAMMS** | Nearest-neighbour reference panel | 30 | Nearest-neighbour selection quality degrades when the reference pool is small |
| **GATK-gCNV** | Probabilistic cohort model | 30 | The Broad recommends ≥ 30 samples for the COHORT mode model to be well-calibrated |
| **CNVkit** | Pooled reference (bias correction) | 10 | A pooled reference with as few as 10 samples provides reasonable GC and mappability correction; larger cohorts improve accuracy |
| **DRAGEN** | Illumina platform normalisation | N/A | Normalisation is handled internally by the DRAGEN pipeline on the ICAv2 platform. ≥ 5 samples (inclusive) are required to create the 'in-run' PoN. ≥ 30 samples required for Targeted CNV Calling |

> **Note:** `canoes_batch_size` and `xhmm_batch_size` control how many BAMs are
> processed per parallel job, not how many samples participate in normalisation.
> All per-batch results are merged into a single cohort-wide matrix before any
> normalisation step runs, so these parameters have no effect on normalisation
> accuracy.

### Disk and I/O considerations
- Intermediate files (depth-of-coverage, read-count matrices, normalised
  coverages) grow proportionally with cohort size.  Ensure the working directory
  has sufficient free space: budget ~2–5 GB per sample.
- Setting `process.stageInMode = 'symlink'` (the default) avoids copying large
  BAM/CRAM files into each work directory.
