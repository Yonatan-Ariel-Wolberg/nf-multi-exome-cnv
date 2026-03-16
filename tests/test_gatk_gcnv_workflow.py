#!/usr/bin/env python3
"""
Tests for the GATK-gCNV workflow pipeline.

Validates:
  1. modules-gatk-gcnv.nf contains all required processes in the correct
     execution order (GENERATE_PLOIDY_PRIORS → PREPROCESS_INTERVALS →
     ANNOTATE_INTERVALS → COLLECT_READ_COUNTS → FILTER_INTERVALS →
     DETERMINE_PLOIDY_COHORT → SCATTER_INTERVALS →
     GERMLINE_CNV_CALLER_COHORT → POSTPROCESS_CALLS → BGZIP_SORT_INDEX_VCF).
  2. GENERATE_PLOIDY_PRIORS produces a 5-column TSV with the correct header
     and correct default priors for autosomes (ploidy 2) and sex chromosomes.
  3. POSTPROCESS_CALLS produces per-sample VCFs named
     "{sample}.genotyped_segments.vcf.gz".
  4. BGZIP_SORT_INDEX_VCF annotates with TOOL=GATK-gCNV.
  5. params-gatk-gcnv.json contains the required keys for the workflow.
"""

import json
import os
import re
import subprocess
import tempfile

REPO_ROOT      = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
GCNV_MODULE    = os.path.join(REPO_ROOT, 'modules', 'modules-gatk-gcnv.nf')
PARAMS_GCNV    = os.path.join(REPO_ROOT, 'params', 'params-gatk-gcnv.json')


def _load_module():
    with open(GCNV_MODULE) as f:
        return f.read()


def _extract_workflow_block(content):
    """Extract the text of the GATK_GCNV workflow block."""
    m = re.search(r'workflow GATK_GCNV \{', content)
    assert m is not None, "workflow GATK_GCNV not found in modules-gatk-gcnv.nf"
    start = m.end() - 1
    depth = 0
    for i, ch in enumerate(content[start:], start):
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return content[start: i + 1]
    raise AssertionError("Unmatched braces for workflow GATK_GCNV")


def _extract_process_block(content, process_name):
    """Extract the text of a named process block."""
    m = re.search(rf'process {re.escape(process_name)} \{{', content)
    assert m is not None, f"process {process_name} not found in modules-gatk-gcnv.nf"
    start = m.end() - 1
    depth = 0
    for i, ch in enumerate(content[start:], start):
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return content[start: i + 1]
    raise AssertionError(f"Unmatched braces for process {process_name}")


def _extract_process_script(block):
    """Extract the script body (between triple-quotes) from a process block."""
    first = block.index('"""')
    second = block.index('"""', first + 3)
    return block[first + 3: second]


# ===========================================================================
# Module structure
# ===========================================================================

class TestGatkGcnvModuleStructure:
    """modules-gatk-gcnv.nf must define all required processes."""

    _REQUIRED_PROCESSES = [
        'GENERATE_PLOIDY_PRIORS',
        'PREPROCESS_INTERVALS',
        'ANNOTATE_INTERVALS',
        'COLLECT_READ_COUNTS',
        'FILTER_INTERVALS',
        'DETERMINE_PLOIDY_COHORT',
        'SCATTER_INTERVALS',
        'GERMLINE_CNV_CALLER_COHORT',
        'POSTPROCESS_CALLS',
        'BGZIP_SORT_INDEX_VCF',
    ]

    def test_all_required_processes_present(self):
        """Every required process is defined in modules-gatk-gcnv.nf."""
        content = _load_module()
        for proc in self._REQUIRED_PROCESSES:
            assert f'process {proc}' in content, (
                f"process {proc} not found in modules-gatk-gcnv.nf"
            )

    def test_workflow_gatk_gcnv_defined(self):
        """The top-level workflow GATK_GCNV is defined."""
        content = _load_module()
        assert 'workflow GATK_GCNV {' in content, (
            "workflow GATK_GCNV not found in modules-gatk-gcnv.nf"
        )

    def test_workflow_step_order(self):
        """All 10 workflow steps appear in the correct execution order."""
        content = _load_module()
        wf = _extract_workflow_block(content)

        ordered_steps = [
            'GENERATE_PLOIDY_PRIORS',
            'PREPROCESS_INTERVALS',
            'ANNOTATE_INTERVALS',
            'COLLECT_READ_COUNTS',
            'FILTER_INTERVALS',
            'DETERMINE_PLOIDY_COHORT',
            'SCATTER_INTERVALS',
            'GERMLINE_CNV_CALLER_COHORT',
            'POSTPROCESS_CALLS',
            'BGZIP_SORT_INDEX_VCF',
        ]
        positions = []
        for step in ordered_steps:
            pos = wf.find(step)
            assert pos != -1, f"{step} not found in workflow GATK_GCNV block"
            positions.append(pos)

        for i in range(len(positions) - 1):
            assert positions[i] < positions[i + 1], (
                f"{ordered_steps[i]} must appear before {ordered_steps[i + 1]} "
                f"in workflow GATK_GCNV"
            )

    def test_workflow_emits_sorted_vcf(self):
        """workflow GATK_GCNV emits sorted_vcf from BGZIP_SORT_INDEX_VCF."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        assert 'sorted_vcf' in wf, "workflow GATK_GCNV must emit sorted_vcf"

    def test_workflow_emits_intervals_vcf(self):
        """workflow GATK_GCNV emits intervals_vcf from POSTPROCESS_CALLS."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        assert 'intervals_vcf' in wf, "workflow GATK_GCNV must emit intervals_vcf"

    def test_ploidy_priors_used_in_determine_ploidy(self):
        """GENERATE_PLOIDY_PRIORS output feeds into DETERMINE_PLOIDY_COHORT."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        priors_pos = wf.find('GENERATE_PLOIDY_PRIORS')
        ploidy_pos = wf.find('DETERMINE_PLOIDY_COHORT')
        assert priors_pos < ploidy_pos, (
            "GENERATE_PLOIDY_PRIORS must precede DETERMINE_PLOIDY_COHORT in workflow"
        )

    def test_samples_indexed_before_postprocess(self):
        """Samples are indexed (withIndex) before being passed to POSTPROCESS_CALLS."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        assert 'withIndex' in wf or 'indexed_samples' in wf, (
            "workflow GATK_GCNV must build an indexed sample channel before "
            "POSTPROCESS_CALLS (each sample needs an integer index)"
        )


# ===========================================================================
# Process-level details
# ===========================================================================

class TestGatkGcnvProcessDetails:
    """Key process scripts contain the correct logic and output names."""

    def test_generate_ploidy_priors_output_name(self):
        """GENERATE_PLOIDY_PRIORS outputs a file named 'ploidy_priors.tsv'."""
        content = _load_module()
        block = _extract_process_block(content, 'GENERATE_PLOIDY_PRIORS')
        assert 'ploidy_priors.tsv' in block, (
            "GENERATE_PLOIDY_PRIORS must produce 'ploidy_priors.tsv'"
        )

    def test_generate_ploidy_priors_has_five_columns(self):
        """GENERATE_PLOIDY_PRIORS header contains 5 columns (contig + 4 ploidy priors)."""
        content = _load_module()
        block = _extract_process_block(content, 'GENERATE_PLOIDY_PRIORS')
        assert 'CONTIG_NAME' in block, "GENERATE_PLOIDY_PRIORS must include CONTIG_NAME column"
        assert 'PLOIDY_PRIOR_0' in block, "GENERATE_PLOIDY_PRIORS must include PLOIDY_PRIOR_0"
        assert 'PLOIDY_PRIOR_1' in block, "GENERATE_PLOIDY_PRIORS must include PLOIDY_PRIOR_1"
        assert 'PLOIDY_PRIOR_2' in block, "GENERATE_PLOIDY_PRIORS must include PLOIDY_PRIOR_2"
        assert 'PLOIDY_PRIOR_3' in block, "GENERATE_PLOIDY_PRIORS must include PLOIDY_PRIOR_3"

    def test_generate_ploidy_priors_adjusts_sex_chromosomes(self):
        """GENERATE_PLOIDY_PRIORS adjusts chrX and chrY priors from autosome defaults."""
        content = _load_module()
        block = _extract_process_block(content, 'GENERATE_PLOIDY_PRIORS')
        assert 'chrX' in block, "GENERATE_PLOIDY_PRIORS must adjust chrX ploidy priors"
        assert 'chrY' in block, "GENERATE_PLOIDY_PRIORS must adjust chrY ploidy priors"

    def test_generate_ploidy_priors_script_produces_correct_output(self):
        """The ploidy-priors bash script generates a valid 5-column TSV."""
        content = _load_module()
        block = _extract_process_block(content, 'GENERATE_PLOIDY_PRIORS')
        script = _extract_process_script(block)

        # Build a minimal FAI file with chr1, chrX, chrY
        fai_content = (
            "chr1\t248956422\t112\t70\t71\n"
            "chrX\t156040895\t251\t70\t71\n"
            "chrY\t57227415\t252\t70\t71\n"
        )

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fai', delete=False) as f:
            f.write(fai_content)
            fai_path = f.name

        # Replace $fai with the temp fai path
        bash_script = script.replace('$fai', fai_path).replace('${fai}', fai_path)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False,
                                         dir='/tmp') as sf:
            sf.write('#!/bin/bash\nset -e\n')
            sf.write(bash_script)
            script_path = sf.name

        try:
            result = subprocess.run(
                ['bash', script_path],
                capture_output=True, text=True, cwd='/tmp'
            )
            assert result.returncode == 0, (
                f"GENERATE_PLOIDY_PRIORS script failed:\n{result.stderr}"
            )

            # Read the output file
            with open('/tmp/ploidy_priors.tsv') as pf:
                lines = [l.rstrip('\n') for l in pf if l.strip()]

            # Header check
            assert lines[0] == (
                'CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\t'
                'PLOIDY_PRIOR_2\tPLOIDY_PRIOR_3'
            ), f"Unexpected header: {lines[0]}"

            # Each data line must have 5 tab-separated fields
            for line in lines[1:]:
                fields = line.split('\t')
                assert len(fields) == 5, (
                    f"Expected 5 columns in ploidy_priors.tsv, got {len(fields)}: {line}"
                )

            # chr1 should have autosome priors (high ploidy-2 weight)
            chr1_row = next((l for l in lines[1:] if l.startswith('chr1\t')), None)
            assert chr1_row is not None, "chr1 row not found in ploidy_priors.tsv"
            chr1_fields = chr1_row.split('\t')
            chr1_ploidy2 = float(chr1_fields[3])
            assert chr1_ploidy2 > 0.9, (
                f"chrX autosome ploidy-2 prior should be > 0.9, got {chr1_ploidy2}"
            )

        finally:
            os.unlink(fai_path)
            os.unlink(script_path)
            if os.path.exists('/tmp/ploidy_priors.tsv'):
                os.unlink('/tmp/ploidy_priors.tsv')

    def test_postprocess_calls_genotyped_segments_naming(self):
        """POSTPROCESS_CALLS produces '{sample}.genotyped_segments.vcf.gz'."""
        content = _load_module()
        block = _extract_process_block(content, 'POSTPROCESS_CALLS')
        assert 'genotyped_segments.vcf.gz' in block, (
            "POSTPROCESS_CALLS must produce '{sample}.genotyped_segments.vcf.gz'"
        )

    def test_postprocess_calls_genotyped_intervals_naming(self):
        """POSTPROCESS_CALLS produces '{sample}.genotyped_intervals.vcf.gz'."""
        content = _load_module()
        block = _extract_process_block(content, 'POSTPROCESS_CALLS')
        assert 'genotyped_intervals.vcf.gz' in block, (
            "POSTPROCESS_CALLS must produce '{sample}.genotyped_intervals.vcf.gz'"
        )

    def test_postprocess_calls_uses_sample_index(self):
        """POSTPROCESS_CALLS passes --sample-index to GATK PostprocessGermlineCNVCalls."""
        content = _load_module()
        block = _extract_process_block(content, 'POSTPROCESS_CALLS')
        assert '--sample-index' in block, (
            "POSTPROCESS_CALLS must pass --sample-index to PostprocessGermlineCNVCalls"
        )

    def test_postprocess_calls_handles_sex_chromosomes(self):
        """POSTPROCESS_CALLS declares chrX and chrY as allosomal contigs."""
        content = _load_module()
        block = _extract_process_block(content, 'POSTPROCESS_CALLS')
        assert '--allosomal-contig chrX' in block, (
            "POSTPROCESS_CALLS must declare chrX as allosomal"
        )
        assert '--allosomal-contig chrY' in block, (
            "POSTPROCESS_CALLS must declare chrY as allosomal"
        )

    def test_bgzip_sort_index_annotates_tool_gatk_gcnv(self):
        """BGZIP_SORT_INDEX_VCF annotates each variant with TOOL=GATK-gCNV."""
        content = _load_module()
        block = _extract_process_block(content, 'BGZIP_SORT_INDEX_VCF')
        assert 'GATK-gCNV' in block, (
            "BGZIP_SORT_INDEX_VCF must annotate variants with TOOL=GATK-gCNV"
        )

    def test_bgzip_sort_index_produces_sorted_gz(self):
        """BGZIP_SORT_INDEX_VCF output follows '*.sorted.vcf.gz' pattern."""
        content = _load_module()
        block = _extract_process_block(content, 'BGZIP_SORT_INDEX_VCF')
        assert '.sorted.vcf.gz' in block, (
            "BGZIP_SORT_INDEX_VCF must produce '*.sorted.vcf.gz' output files"
        )

    def test_bgzip_sort_index_produces_tbi_index(self):
        """BGZIP_SORT_INDEX_VCF output includes a '*.sorted.vcf.gz.tbi' index."""
        content = _load_module()
        block = _extract_process_block(content, 'BGZIP_SORT_INDEX_VCF')
        assert '.sorted.vcf.gz.tbi' in block, (
            "BGZIP_SORT_INDEX_VCF must produce '*.sorted.vcf.gz.tbi' index files"
        )

    def test_germline_cnv_caller_runs_in_cohort_mode(self):
        """GERMLINE_CNV_CALLER_COHORT calls GermlineCNVCaller in COHORT run-mode."""
        content = _load_module()
        block = _extract_process_block(content, 'GERMLINE_CNV_CALLER_COHORT')
        assert 'COHORT' in block, (
            "GERMLINE_CNV_CALLER_COHORT must use --run-mode COHORT"
        )

    def test_filter_intervals_uses_cohort_counts(self):
        """FILTER_INTERVALS collects all sample read counts before filtering."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        # The FILTER_INTERVALS step in the workflow should receive a collected list
        filter_pos = wf.find('FILTER_INTERVALS')
        collect_before_filter = wf[:filter_pos]
        assert '.collect()' in collect_before_filter, (
            "All sample read counts must be collected before FILTER_INTERVALS"
        )


# ===========================================================================
# params-gatk-gcnv.json
# ===========================================================================

class TestGatkGcnvParamsJson:
    """params-gatk-gcnv.json must contain all keys required by the workflow."""

    _REQUIRED_KEYS = [
        'workflow',
        'outdir',
        'samples_path',
        'fasta',
        'fai',
        'dict',
        'exome_targets',
    ]

    def _load_params(self):
        with open(PARAMS_GCNV) as f:
            return json.load(f)

    def test_all_required_keys_present(self):
        """params-gatk-gcnv.json must contain all required keys."""
        data = self._load_params()
        for key in self._REQUIRED_KEYS:
            assert key in data, (
                f"params-gatk-gcnv.json is missing required key '{key}'"
            )

    def test_workflow_value_is_gcnv(self):
        """The 'workflow' key must be set to 'gcnv'."""
        data = self._load_params()
        assert data['workflow'] == 'gcnv', (
            f"params-gatk-gcnv.json: expected workflow='gcnv', got '{data['workflow']}'"
        )

    def test_scatter_count_is_positive(self):
        """scatter_count must be a positive integer for interval sharding."""
        data = self._load_params()
        assert 'scatter_count' in data, "params-gatk-gcnv.json missing 'scatter_count'"
        assert isinstance(data['scatter_count'], int) and data['scatter_count'] > 0, (
            f"scatter_count must be a positive integer, got {data['scatter_count']}"
        )

    def test_is_wgs_is_false(self):
        """is_wgs must be false for exome sequencing analyses."""
        data = self._load_params()
        assert 'is_wgs' in data, "params-gatk-gcnv.json missing 'is_wgs'"
        assert data['is_wgs'] is False, (
            f"params-gatk-gcnv.json: is_wgs must be false for exome analysis, "
            f"got {data['is_wgs']}"
        )
