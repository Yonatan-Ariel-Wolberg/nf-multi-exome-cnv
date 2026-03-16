# nf-multicaller-exomecnv

# Tools

**CANOES**: https://github.com/ShenLab/CANOES

**CLAMMS**: https://github.com/rgcgithub/clamms

**XHMM**: https://github.com/RRafiee/XHMM

**InDelible**: https://github.com/HurlesGroupSanger/indelible

**ICAv2**: https://help.ica.illumina.com/command-line-interface/cli-indexcommands

**DRAGEN Germline Enrichment**: https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/dragen-enrichment.html

**GATK-gCNV**: https://github.com/broadinstitute/gatk

**CNVkit**: https://github.com/etal/cnvkit

**Truvari**: https://github.com/ACEnglish/truvari

**SURVIVOR**: https://github.com/fritzsedlazeck/SURVIVOR


# nf-multicaller-exomecnv

A Nextflow pipeline for comprehensive genomic analysis combining multiple variant callers with exome and copy number variation (CNV) detection. Based on a workflow by **Dr Phelelani T Mpangase**: https://github.com/phelelani/nf-exomecnv. The modules-icav2-dragen.nf module is based off a script by **Reagan Cannell**: https://github.com/SBIMB/ica-elwazi/tree/main/nextflow_workflows/cram_input_dragen_ica_workflow.

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

## Building the icav2_cli_v2.43.0.sif
```bash
apptainer build icav2_cli_v2.43.0.sif icav2.def
```
