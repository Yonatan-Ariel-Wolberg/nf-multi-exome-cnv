# nf-multicaller-exomecnv

# Tools

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


# nf-multicaller-exomecnv

A Nextflow pipeline for comprehensive genomic analysis combining multiple variant callers with exome and copy number variation (CNV) detection. Based on a workflow by **Dr Phelelani T Mpangase** (phelelani): https://github.com/phelelani/nf-exomecnv. The modules-icav2-dragen.nf module is based off a script by **Regan Cannell** (RWCannell): https://github.com/SBIMB/ica-elwazi/tree/main/nextflow_workflows/cram_input_dragen_ica_workflow.

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
batch per chromosome; the per-batch matrices are then merged.

- **`canoes_batch_size`** (default `100`; override in params file or via profile):
  - Up to 50 samples → set to `50` (1 batch; minimal overhead).
  - Up to 300 samples → set to `100` (3 batches per chromosome).
  - Up to 2000 samples → set to `200` (10 batches per chromosome).

#### XHMM (`--workflow xhmm`)
XHMM uses GATK `DepthOfCoverage`, which can require significant memory per job
when given a large BAM list.  The BAM list is split into fixed-size groups
(`xhmm_batch_size`) and one `GATK_DOC` job is submitted per group.

- **`xhmm_batch_size`** (default `50`; override in params file or via profile):
  - Up to 50 samples → set to `25`–`50` (1–2 batches).
  - Up to 300 samples → set to `50` (up to 6 batches).
  - Up to 2000 samples → set to `100` (up to 20 batches; keeps memory per job manageable).

#### CLAMMS (`--workflow clamms`)
All CLAMMS steps (depth-of-coverage, normalisation, model training, CNV calling)
are fully per-sample and run in parallel automatically via Nextflow's data-flow
model.  No batch-size tuning is needed; `process.maxForks` limits the number of
concurrent per-sample jobs.

#### GATK-gCNV (`--workflow gcnv`)
gCNV trains a cohort-level model and is naturally designed for large cohorts.
Interval parallelism is controlled by `scatter_count` in
`params/params-gatk-gcnv.json` (default `5000`).  For the cohort-level steps
(`FilterIntervals`, `DetermineGermlineContigPloidy`, `GermlineCNVCaller`) that
collect all samples at once, increase the memory/CPU limits using the `large`
profile when running hundreds or thousands of samples.

#### CNVkit (`--workflow cnvkit`)
CNVkit processes samples independently; concurrency is governed by
`process.maxForks`.  Use `test_size` (integer) or `test_list` (comma-separated
sample IDs) in the params file to restrict a run to a subset of samples during
development.

#### DRAGEN / ICAv2 (`--workflow dragen`)
DRAGEN runs are submitted to the ICAv2 platform.  `maxUploadForks` (default `4`)
limits simultaneous CRAM uploads.  Increase this if your network bandwidth and
ICA account limits allow it.

#### SURVIVOR / Truvari (`--workflow survivor` / `--workflow truvari`)
These consensus modules are fast and lightweight regardless of cohort size; no
additional tuning is required.

### Disk and I/O considerations
- Intermediate files (depth-of-coverage, read-count matrices, normalised
  coverages) grow proportionally with cohort size.  Ensure the working directory
  has sufficient free space: budget ~2–5 GB per sample.
- Setting `process.stageInMode = 'symlink'` (the default) avoids copying large
  BAM/CRAM files into each work directory.
