#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <ddd-africa|ddd-uk> <5|10|all>" >&2
  exit 1
fi

DATASET="$1"
SCOPE="$2"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

DDD_DATA_HOME="${DDD_DATA_HOME:-/home/${USER}/DECIPHERING_DD_DATA}"

case "$DATASET" in
  ddd-africa)
    PARAM_DIR="params/ddd-africa"
    SAMPLE_SHEET="${DDD_AFRICA_SAMPLE_SHEET:-${DDD_DATA_HOME}/DDD_AFRICA_DATA/samplesheet_africa.tsv}"
    BAM_ROOT="${DDD_AFRICA_BAM_ROOT:-${DDD_DATA_HOME}/DDD_AFRICA_DATA/batch_3/organized_data}"
    INDELIBLE_ROOT="${DDD_AFRICA_INDELIBLE_ROOT:-${DDD_DATA_HOME}/DDD_AFRICA_DATA/batch_3/organized_data}"
    DRAGEN_ROOT="${DDD_AFRICA_DRAGEN_ROOT:-${DDD_DATA_HOME}/DDD_AFRICA_DATA/batch_3/organized_data/Proband}"
    ;;
  ddd-uk)
    PARAM_DIR="params/ddd-uk"
    SAMPLE_SHEET="${DDD_UK_SAMPLE_SHEET:-${DDD_DATA_HOME}/DDD_UK_DATA/samplesheet_uk.tsv}"
    BAM_ROOT="${DDD_UK_BAM_ROOT:-${DDD_DATA_HOME}/DDD_UK_DATA/bams}"
    INDELIBLE_ROOT="${DDD_UK_INDELIBLE_ROOT:-${DDD_DATA_HOME}/DDD_UK_DATA/crams}"
    DRAGEN_ROOT="${DDD_UK_DRAGEN_ROOT:-${DDD_DATA_HOME}/DDD_UK_DATA/crams}"
    ;;
  *)
    echo "Unsupported dataset: $DATASET (expected ddd-africa or ddd-uk)" >&2
    exit 1
    ;;
esac

case "$SCOPE" in
  5|10)
    SAMPLE_LIMIT="$SCOPE"
    COHORT_PROFILE="small"
    SCOPE_LABEL="${SCOPE}-samples"
    ;;
  all)
    SAMPLE_LIMIT=""
    COHORT_PROFILE="medium"
    SCOPE_LABEL="all-samples"
    ;;
  *)
    echo "Unsupported scope: $SCOPE (expected 5, 10, or all)" >&2
    exit 1
    ;;
esac

# For subset samplesheets, set SAMPLESHEET_HAS_HEADER=true/false to force behavior.
# Default "auto" detects a header when the first row contains both "sample" and "bam".
SAMPLESHEET_HAS_HEADER="${SAMPLESHEET_HAS_HEADER:-auto}"

RUN_ROOT="$REPO_DIR/results/sbatch/${DATASET}/${SCOPE_LABEL}"
TMP_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
TMP_ROOT="${TMP_BASE}/nf-multicaller-exomecnv/${DATASET}/${SCOPE_LABEL}/${SLURM_JOB_ID:-manual}"
mkdir -p "$RUN_ROOT" "$TMP_ROOT"

if command -v module >/dev/null 2>&1; then
  module load nextflow || true
  module load apptainer || true
fi

run_nextflow() {
  local workflow="$1"
  local params_file="$2"
  shift 2

  printf '\n=== Running workflow: %s (%s, %s) ===\n' "$workflow" "$DATASET" "$SCOPE_LABEL"
  nextflow run main.nf \
    -profile "wits,${COHORT_PROFILE}" \
    --workflow "$workflow" \
    -params-file "$params_file" \
    --outdir "$RUN_ROOT/$workflow" \
    -resume \
    "$@"
}

make_subset_samplesheet() {
  local destination="$1"
  local first_line
  local has_header=0
  case "$SAMPLESHEET_HAS_HEADER" in
    true|TRUE|1)
      has_header=1
      ;;
    false|FALSE|0)
      has_header=0
      ;;
    auto|AUTO)
      first_line="$(head -n 1 "$SAMPLE_SHEET" || true)"
      if [[ "$first_line" =~ [Ss]ample ]] && [[ "$first_line" =~ [Bb][Aa][Mm] ]]; then
        has_header=1
      fi
      ;;
    *)
      echo "Invalid SAMPLESHEET_HAS_HEADER value: $SAMPLESHEET_HAS_HEADER (use auto, true, or false)" >&2
      exit 1
      ;;
  esac
  head -n "$((SAMPLE_LIMIT + has_header))" "$SAMPLE_SHEET" > "$destination"
  if [[ ! -s "$destination" ]]; then
    echo "Subset samplesheet is empty: $destination" >&2
    exit 1
  fi
}

make_subset_bams_for_gcnv() {
  local subset_dir="$1"
  mkdir -p "$subset_dir"
  mapfile -d '' bam_files < <(find "$BAM_ROOT" -type f -name '*.bam' -print0 | sort -z)

  if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "No BAM files found under $BAM_ROOT" >&2
    exit 1
  fi

  local max_count="$SAMPLE_LIMIT"
  if (( ${#bam_files[@]} < max_count )); then
    max_count="${#bam_files[@]}"
  fi

  for ((i = 0; i < max_count; i++)); do
    local bam="${bam_files[$i]}"
    local bam_name
    bam_name="$(basename "$bam")"
    ln -sf "$bam" "$subset_dir/$bam_name"

    local bai="${bam}.bai"
    if [[ ! -f "$bai" ]]; then
      bai="${bam%.bam}.bai"
    fi
    if [[ -f "$bai" ]]; then
      ln -sf "$bai" "$subset_dir/${bam_name}.bai"
    fi
  done
}

make_subset_crams_for_indelible() {
  local subset_dir="$1"
  mkdir -p "$subset_dir"
  mapfile -d '' probands < <(find "$INDELIBLE_ROOT" -type f -name '*_01_1.cram' -print0 | sort -z)

  if [[ ${#probands[@]} -eq 0 ]]; then
    echo "No proband CRAM files (*_01_1.cram) found under $INDELIBLE_ROOT" >&2
    exit 1
  fi

  local max_count="$SAMPLE_LIMIT"
  if (( ${#probands[@]} < max_count )); then
    max_count="${#probands[@]}"
  fi

  for ((i = 0; i < max_count; i++)); do
    local proband="${probands[$i]}"
    local stem="${proband%_01_1.cram}"
    for rel in _01_1 _02_2 _03_3; do
      for ext in cram cram.crai; do
        local f="${stem}${rel}.${ext}"
        if [[ -f "$f" ]]; then
          ln -sf "$f" "$subset_dir/$(basename "$f")"
        fi
      done
    done
  done
}

make_subset_crams_for_dragen() {
  local subset_dir="$1"
  mkdir -p "$subset_dir"
  mapfile -d '' cram_files < <(find "$DRAGEN_ROOT" -type f -name '*.cram' -print0 | sort -z)

  if [[ ${#cram_files[@]} -eq 0 ]]; then
    echo "No CRAM files found under $DRAGEN_ROOT" >&2
    exit 1
  fi

  local max_count="$SAMPLE_LIMIT"
  if (( ${#cram_files[@]} < max_count )); then
    max_count="${#cram_files[@]}"
  fi

  for ((i = 0; i < max_count; i++)); do
    local cram="${cram_files[$i]}"
    local cram_name
    cram_name="$(basename "$cram")"
    ln -sf "$cram" "$subset_dir/$cram_name"

    local crai="${cram}.crai"
    if [[ -f "$crai" ]]; then
      ln -sf "$crai" "$subset_dir/${cram_name}.crai"
    fi
  done
}

if [[ -n "$SAMPLE_LIMIT" ]]; then
  SUBSET_SAMPLESHEET="$TMP_ROOT/samplesheet_${SAMPLE_LIMIT}.tsv"
  make_subset_samplesheet "$SUBSET_SAMPLESHEET"

  GCNV_SUBSET_DIR="$TMP_ROOT/gcnv_bams"
  make_subset_bams_for_gcnv "$GCNV_SUBSET_DIR"

  INDELIBLE_SUBSET_DIR="$TMP_ROOT/indelible_crams"
  make_subset_crams_for_indelible "$INDELIBLE_SUBSET_DIR"

  DRAGEN_SUBSET_DIR="$TMP_ROOT/dragen_crams"
  make_subset_crams_for_dragen "$DRAGEN_SUBSET_DIR"
fi

if [[ -n "$SAMPLE_LIMIT" ]]; then
  run_nextflow canoes    "$PARAM_DIR/params-canoes-wits-${DATASET}.json" --samplesheet_bams "$SUBSET_SAMPLESHEET"
  run_nextflow xhmm      "$PARAM_DIR/params-xhmm-wits-${DATASET}.json" --samplesheet_bams "$SUBSET_SAMPLESHEET"
  run_nextflow clamms    "$PARAM_DIR/params-clamms-wits-${DATASET}.json" --samplesheet_bams "$SUBSET_SAMPLESHEET"
  run_nextflow cnvkit    "$PARAM_DIR/params-cnvkit-wits-${DATASET}.json" --test_size "$SAMPLE_LIMIT" --test_list all
  run_nextflow gcnv      "$PARAM_DIR/params-gatk-gcnv-wits-${DATASET}.json" --samples_path "$GCNV_SUBSET_DIR/*.{bam,bam.bai}"
  run_nextflow indelible "$PARAM_DIR/params-indelible-wits-${DATASET}.json" --crams "$INDELIBLE_SUBSET_DIR"
  run_nextflow dragen    "$PARAM_DIR/params-icav2-dragen-wits-${DATASET}.json" --cramFilePairsUploadPath "$DRAGEN_SUBSET_DIR/*.{cram,cram.crai}"
else
  run_nextflow canoes    "$PARAM_DIR/params-canoes-wits-${DATASET}.json"
  run_nextflow xhmm      "$PARAM_DIR/params-xhmm-wits-${DATASET}.json"
  run_nextflow clamms    "$PARAM_DIR/params-clamms-wits-${DATASET}.json"
  run_nextflow cnvkit    "$PARAM_DIR/params-cnvkit-wits-${DATASET}.json"
  run_nextflow gcnv      "$PARAM_DIR/params-gatk-gcnv-wits-${DATASET}.json"
  run_nextflow indelible "$PARAM_DIR/params-indelible-wits-${DATASET}.json"
  run_nextflow dragen    "$PARAM_DIR/params-icav2-dragen-wits-${DATASET}.json"
fi

printf '\nCompleted all caller workflows for %s (%s). Results: %s\n' "$DATASET" "$SCOPE_LABEL" "$RUN_ROOT"
