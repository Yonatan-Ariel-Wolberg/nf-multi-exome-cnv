# nf-multicaller-exomecnv

A Nextflow pipeline for comprehensive genomic analysis combining multiple variant callers with exome and copy number variation (CNV) detection. Based on a workflow by **Dr Phelelani T Mpangase** (phelelani): https://github.com/phelelani/nf-exomecnv. The modules-icav2-dragen.nf module is based off a script by **Regan Cannell** (RWCannell): https://github.com/SBIMB/ica-elwazi/tree/main/nextflow_workflows/cram_input_dragen_ica_workflow.

## Tools

**CANOES**: https://github.com/ShenLab/CANOES

**CLAMMS**: https://github.com/rgcgithub/clamms

**XHMM**: https://github.com/RRafiee/XHMM

**InDelible**: https://github.com/HurlesGroupSanger/indelible

**ICAv2**: https://help.ica.illumina.com/command-line-interface/cli-indexcommands

**Illumina's DRAGEN Germline Enrichment**: https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/dragen-enrichment.html

**The Broad Institute's GATK-gCNV**: https://github.com/broadinstitute/gatk

**CNVkit**: https://github.com/etal/cnvkit

**Truvari**: https://github.com/ACEnglish/truvari

**SURVIVOR**: https://github.com/fritzsedlazeck/SURVIVOR

## Overview

This pipeline implements a robust bioinformatics workflow for analyzing genomic data using:
- Multiple variant calling approaches
- Exome analysis
- Copy number variation (CNV) detection

Built with Nextflow for reproducible, scalable, and portable data analysis.

## Features

- **Multi-caller approach**: Integrates multiple variant calling algorithms for robust variant detection
- **Exome analysis**: Specialized processing for whole exome sequencing (WES) data
- **CNV detection**: Comprehensive copy number variation analysis
- **Nextflow-based**: Ensures reproducibility and scalability across different computing environments
- **Containerized**: Includes Docker/Singularity support for dependency management

## Requirements

### Software
- **Nextflow**: >= 20.04
- **Java**: >= 8

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

2. Build the icav2_cli_v2.43.0.sif
```bash
apptainer build icav2_cli_v2.43.0.sif icav2.def
```

## Usage

Each workflow is run by providing `--workflow <name>` along with a params file. Edit the corresponding `params/*.json` file to point to your actual input files before running.

### Individual CNV callers

The individual callers should be run first. Their VCF outputs feed into the consensus modules (SURVIVOR / Truvari).

```bash
# CANOES — exome CNV caller using read-depth and PCA normalisation
nextflow run main.nf --workflow canoes -params-file params/params-canoes.json

# XHMM — exome CNV caller using GATK depth-of-coverage and PCA
nextflow run main.nf --workflow xhmm -params-file params/params-xhmm.json

# CLAMMS — exome CNV caller using nearest-neighbour normalisation
nextflow run main.nf --workflow clamms -params-file params/params-clamms.json

# GATK-gCNV — cohort-mode germline CNV caller from the Broad Institute
nextflow run main.nf --workflow gcnv -params-file params/params-gatk-gcnv.json

# CNVkit — coverage-based CNV caller with pooled reference normalisation
nextflow run main.nf --workflow cnvkit -params-file params/params-cnvkit.json

# DRAGEN — Illumina DRAGEN germline enrichment pipeline via ICAv2
nextflow run main.nf --workflow dragen -params-file params/params-icav2-dragen.json

# InDelible — small insertion/deletion caller from split reads
nextflow run main.nf --workflow indelible -params-file params/params-indelible.json
```

### Consensus modules (run after individual callers)

Point the `*_dir` paths in the params file to the VCF output directories produced by the callers above, then run:

```bash
# SURVIVOR — merge per-caller VCFs into a consensus call set
nextflow run main.nf --workflow survivor -params-file params/params-survivor.json

# Truvari — collapse per-caller VCFs into a consensus call set
nextflow run main.nf --workflow truvari -params-file params/params-truvari.json
```

### Parameter files

All required and optional parameters for each workflow are documented in the corresponding `params/*.json` template:

| Params file | Workflow |
|---|---|
| `params/params-canoes.json` | `--workflow canoes` |
| `params/params-xhmm.json` | `--workflow xhmm` |
| `params/params-clamms.json` | `--workflow clamms` |
| `params/params-gatk-gcnv.json` | `--workflow gcnv` |
| `params/params-cnvkit.json` | `--workflow cnvkit` |
| `params/params-icav2-dragen.json` | `--workflow dragen` |
| `params/params-indelible.json` | `--workflow indelible` |
| `params/params-survivor.json` | `--workflow survivor` |
| `params/params-truvari.json` | `--workflow truvari` |

## Scalability: running 50, 300, or 2000 samples

The pipeline ships with three cohort-size profiles that pre-configure the
scheduler, concurrency, and per-caller batch sizes for small, medium, and large
runs.  Combine a cohort-size profile with the site profile using `-profile`:

```bash
# ≤ 50 samples
nextflow run main.nf -profile wits,small  -params-file params/params-canoes.json

# 50–300 samples (default behaviour if no cohort-size profile is specified)
nextflow run main.nf -profile wits,medium -params-file params/params-canoes.json

# 300–2000 samples
nextflow run main.nf -profile wits,large  -params-file params/params-canoes.json
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
cohort-wide matrix **before** the CANOES R script performs GC correction and
PCA normalisation.  The batch size therefore controls **only parallelism and
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
`params/params-gatk-gcnv.json` (default `5000`).  For the cohort-level steps
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
| **CANOES** | GC correction + PCA | 30 | PCA becomes unreliable below ~20 samples; 30+ recommended for typical exome panels |
| **XHMM** | Mean-centering + PCA | 30 | XHMM's own filters may exclude all targets/samples with very small cohorts; 50+ strongly recommended |
| **CLAMMS** | Nearest-neighbour reference panel | 30 | Nearest-neighbour selection quality degrades when the reference pool is small |
| **GATK-gCNV** | Probabilistic cohort model | 30 | The Broad recommends ≥ 30 samples for the COHORT mode model to be well-calibrated |
| **CNVkit** | Pooled reference (bias correction) | 10 | A pooled reference with as few as 10 samples provides reasonable GC and mappability correction; larger cohorts improve accuracy |
| **DRAGEN** | Illumina platform normalisation | N/A | Normalisation is handled internally by the DRAGEN pipeline on the ICAv2 platform |

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
