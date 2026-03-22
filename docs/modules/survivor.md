# SURVIVOR module (`--workflow survivor`)

## Methodology
SURVIVOR performs distance-based CNV merging across multiple callers per sample and produces union/intersection consensus callsets.

## Processes in this module
1. **RUN_SURVIVOR_MERGE (union)** – merge calls with minimum support of one caller.
2. **RUN_SURVIVOR_MERGE (intersection)** – merge calls with minimum support of two callers.

## Required parameters
- `workflow`: `survivor`
- `outdir`
- At least two caller VCF directories from:
  - `canoes_dir`
  - `clamms_dir`
  - `xhmm_dir`
  - `cnvkit_dir`
  - `gcnv_dir`
  - `dragen_dir`
  - `indelible_dir`

## Params template
- `params/params-survivor.json`
