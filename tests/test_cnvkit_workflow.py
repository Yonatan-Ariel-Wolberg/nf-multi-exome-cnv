#!/usr/bin/env python3
"""
Tests for the CNVkit workflow pipeline.

Validates:
  1. modules-cnvkit.nf contains all required processes in the correct
     execution order (GENERATE_ACCESS → AUTOBIN → COVERAGE →
     CREATE_POOLED_REFERENCE → CALL_CNV → EXPORT_RESULTS →
     BGZIP_SORT_INDEX_VCF).
  2. EXPORT_RESULTS produces per-sample VCFs named "{sample_id}_CNVKIT_output.vcf".
  3. BGZIP_SORT_INDEX_VCF annotates with TOOL=CNVkit.
  4. The is_large_run() helper function is defined and referenced for memory
     scaling.
  5. params-cnvkit.json contains the required keys for the workflow.
"""

import json
import os
import re

REPO_ROOT     = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
CNVKIT_MODULE = os.path.join(REPO_ROOT, 'modules', 'modules-cnvkit.nf')
PARAMS_CNVKIT = os.path.join(REPO_ROOT, 'params', 'params-cnvkit.json')


def _load_module():
    with open(CNVKIT_MODULE) as f:
        return f.read()


def _extract_workflow_block(content):
    """Extract the text of the CNVKIT workflow block."""
    m = re.search(r'workflow CNVKIT \{', content)
    assert m is not None, "workflow CNVKIT not found in modules-cnvkit.nf"
    start = m.end() - 1
    depth = 0
    for i, ch in enumerate(content[start:], start):
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return content[start: i + 1]
    raise AssertionError("Unmatched braces for workflow CNVKIT")


def _extract_process_block(content, process_name):
    """Extract the text of a named process block."""
    m = re.search(rf'process {re.escape(process_name)} \{{', content)
    assert m is not None, f"process {process_name} not found in modules-cnvkit.nf"
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


# ===========================================================================
# Module structure
# ===========================================================================

class TestCnvkitModuleStructure:
    """modules-cnvkit.nf must define all required processes."""

    _REQUIRED_PROCESSES = [
        'GENERATE_ACCESS',
        'AUTOBIN',
        'COVERAGE',
        'CREATE_POOLED_REFERENCE',
        'CALL_CNV',
        'EXPORT_RESULTS',
        'BGZIP_SORT_INDEX_VCF',
    ]

    def test_all_required_processes_present(self):
        """Every required process is defined in modules-cnvkit.nf."""
        content = _load_module()
        for proc in self._REQUIRED_PROCESSES:
            assert f'process {proc}' in content, (
                f"process {proc} not found in modules-cnvkit.nf"
            )

    def test_workflow_cnvkit_defined(self):
        """The top-level workflow CNVKIT is defined."""
        content = _load_module()
        assert 'workflow CNVKIT {' in content, "workflow CNVKIT not found in modules-cnvkit.nf"

    def test_workflow_step_order(self):
        """All 7 workflow steps appear in the correct execution order."""
        content = _load_module()
        wf = _extract_workflow_block(content)

        ordered_steps = [
            'GENERATE_ACCESS',
            'AUTOBIN',
            'COVERAGE',
            'CREATE_POOLED_REFERENCE',
            'CALL_CNV',
            'EXPORT_RESULTS',
            'BGZIP_SORT_INDEX_VCF',
        ]
        positions = []
        for step in ordered_steps:
            pos = wf.find(step)
            assert pos != -1, f"{step} not found in workflow CNVKIT block"
            positions.append(pos)

        for i in range(len(positions) - 1):
            assert positions[i] < positions[i + 1], (
                f"{ordered_steps[i]} must appear before {ordered_steps[i + 1]} "
                f"in workflow CNVKIT"
            )

    def test_workflow_emits_sorted_vcf(self):
        """workflow CNVKIT emits sorted_vcf from BGZIP_SORT_INDEX_VCF."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        assert 'sorted_vcf' in wf, "workflow CNVKIT must emit sorted_vcf"

    def test_workflow_emits_bed(self):
        """workflow CNVKIT emits bed from EXPORT_RESULTS."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        assert 'bed' in wf, "workflow CNVKIT must emit bed results"

    def test_pooled_reference_uses_all_coverages(self):
        """CREATE_POOLED_REFERENCE step collects all coverage files before building."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        # The workflow should collect coverage from all samples before creating the reference
        assert '.collect()' in wf, (
            "workflow CNVKIT must call .collect() to gather all sample coverage "
            "files before CREATE_POOLED_REFERENCE"
        )

    def test_coverage_joined_before_call_cnv(self):
        """Target and antitarget coverage channels are joined before CALL_CNV."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        join_pos = wf.find('.join(')
        call_pos = wf.find('CALL_CNV')
        assert join_pos != -1, "Coverage channels must be joined before CALL_CNV"
        assert join_pos < call_pos, (
            "Coverage join must appear before CALL_CNV in workflow CNVKIT"
        )


# ===========================================================================
# Process-level details
# ===========================================================================

class TestCnvkitProcessDetails:
    """Key process scripts contain the correct logic and output names."""

    def test_export_results_vcf_naming(self):
        """EXPORT_RESULTS produces VCFs named '{sample_id}_CNVKIT_output.vcf'."""
        content = _load_module()
        block = _extract_process_block(content, 'EXPORT_RESULTS')
        assert '_CNVKIT_output.vcf' in block, (
            "EXPORT_RESULTS must produce VCFs named '{sample_id}_CNVKIT_output.vcf'"
        )

    def test_export_results_uses_cnvkit_export_vcf(self):
        """EXPORT_RESULTS calls 'cnvkit.py export vcf' to produce VCF output."""
        content = _load_module()
        block = _extract_process_block(content, 'EXPORT_RESULTS')
        assert 'cnvkit.py export vcf' in block, (
            "EXPORT_RESULTS must call 'cnvkit.py export vcf'"
        )

    def test_export_results_produces_bed(self):
        """EXPORT_RESULTS also calls 'cnvkit.py export bed' to produce BED output."""
        content = _load_module()
        block = _extract_process_block(content, 'EXPORT_RESULTS')
        assert 'cnvkit.py export bed' in block, (
            "EXPORT_RESULTS must call 'cnvkit.py export bed'"
        )

    def test_bgzip_sort_index_annotates_tool_cnvkit(self):
        """BGZIP_SORT_INDEX_VCF annotates each variant with TOOL=CNVkit."""
        content = _load_module()
        block = _extract_process_block(content, 'BGZIP_SORT_INDEX_VCF')
        assert 'CNVkit' in block, (
            "BGZIP_SORT_INDEX_VCF must annotate variants with TOOL=CNVkit"
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

    def test_generate_access_uses_cnvkit_access(self):
        """GENERATE_ACCESS calls 'cnvkit.py access' on the reference FASTA."""
        content = _load_module()
        block = _extract_process_block(content, 'GENERATE_ACCESS')
        assert 'cnvkit.py access' in block, (
            "GENERATE_ACCESS must call 'cnvkit.py access'"
        )

    def test_generate_access_output_is_access_bed(self):
        """GENERATE_ACCESS produces 'access.hg38.bed' as output."""
        content = _load_module()
        block = _extract_process_block(content, 'GENERATE_ACCESS')
        assert 'access.hg38.bed' in block, (
            "GENERATE_ACCESS must output 'access.hg38.bed'"
        )

    def test_autobin_uses_cnvkit_autobin(self):
        """AUTOBIN calls 'cnvkit.py autobin'."""
        content = _load_module()
        block = _extract_process_block(content, 'AUTOBIN')
        assert 'cnvkit.py autobin' in block, (
            "AUTOBIN must call 'cnvkit.py autobin'"
        )

    def test_coverage_uses_cnvkit_coverage(self):
        """COVERAGE calls 'cnvkit.py coverage' for target and antitarget."""
        content = _load_module()
        block = _extract_process_block(content, 'COVERAGE')
        assert 'cnvkit.py coverage' in block, (
            "COVERAGE must call 'cnvkit.py coverage'"
        )

    def test_call_cnv_uses_fix_and_segment(self):
        """CALL_CNV uses 'cnvkit.py fix' and 'cnvkit.py segment'."""
        content = _load_module()
        block = _extract_process_block(content, 'CALL_CNV')
        assert 'cnvkit.py fix' in block, "CALL_CNV must call 'cnvkit.py fix'"
        assert 'cnvkit.py segment' in block, "CALL_CNV must call 'cnvkit.py segment'"

    def test_call_cnv_produces_cnr_and_cns(self):
        """CALL_CNV produces '.cnr' (ratios) and '.cns' (segments) output files."""
        content = _load_module()
        block = _extract_process_block(content, 'CALL_CNV')
        assert '.cnr' in block, "CALL_CNV must produce a '.cnr' ratios file"
        assert '.cns' in block, "CALL_CNV must produce a '.cns' segments file"

    def test_create_pooled_reference_uses_cnvkit_reference(self):
        """CREATE_POOLED_REFERENCE calls 'cnvkit.py reference'."""
        content = _load_module()
        block = _extract_process_block(content, 'CREATE_POOLED_REFERENCE')
        assert 'cnvkit.py reference' in block, (
            "CREATE_POOLED_REFERENCE must call 'cnvkit.py reference'"
        )

    def test_create_pooled_reference_output_is_pooled_cnn(self):
        """CREATE_POOLED_REFERENCE produces 'pooled_reference.cnn' as output."""
        content = _load_module()
        block = _extract_process_block(content, 'CREATE_POOLED_REFERENCE')
        assert 'pooled_reference.cnn' in block, (
            "CREATE_POOLED_REFERENCE must output 'pooled_reference.cnn'"
        )


# ===========================================================================
# is_large_run helper
# ===========================================================================

class TestCnvkitIsLargeRun:
    """The is_large_run() function must be defined and used for memory scaling."""

    def test_is_large_run_defined(self):
        """The is_large_run() helper function must be defined in modules-cnvkit.nf."""
        content = _load_module()
        assert 'def is_large_run()' in content, (
            "modules-cnvkit.nf must define is_large_run() for dynamic memory allocation"
        )

    def test_is_large_run_used_in_coverage(self):
        """is_large_run() is used in the COVERAGE process for memory scaling."""
        content = _load_module()
        block = _extract_process_block(content, 'COVERAGE')
        assert 'is_large_run()' in block, (
            "COVERAGE must use is_large_run() to dynamically set memory"
        )

    def test_is_large_run_used_in_create_pooled_reference(self):
        """is_large_run() is used in CREATE_POOLED_REFERENCE for memory scaling."""
        content = _load_module()
        block = _extract_process_block(content, 'CREATE_POOLED_REFERENCE')
        assert 'is_large_run()' in block, (
            "CREATE_POOLED_REFERENCE must use is_large_run() to dynamically set memory"
        )

    def test_is_large_run_references_test_size(self):
        """is_large_run() inspects params.test_size to determine run scale."""
        content = _load_module()
        m = re.search(r'def is_large_run\(\).*?\}', content, re.DOTALL)
        assert m is not None, "is_large_run() function body not found"
        func_body = m.group(0)
        assert 'test_size' in func_body, (
            "is_large_run() must reference params.test_size"
        )


# ===========================================================================
# params-cnvkit.json
# ===========================================================================

class TestCnvkitParamsJson:
    """params-cnvkit.json must contain all keys required by the CNVkit workflow."""

    _REQUIRED_KEYS = [
        'workflow',
        'outdir',
        'bams',
        'fasta',
        'targets',
        'refflat',
    ]

    def _load_params(self):
        with open(PARAMS_CNVKIT) as f:
            return json.load(f)

    def test_all_required_keys_present(self):
        """params-cnvkit.json must contain all required keys."""
        data = self._load_params()
        for key in self._REQUIRED_KEYS:
            assert key in data, (
                f"params-cnvkit.json is missing required key '{key}'"
            )

    def test_workflow_value_is_cnvkit(self):
        """The 'workflow' key must be set to 'cnvkit'."""
        data = self._load_params()
        assert data['workflow'] == 'cnvkit', (
            f"params-cnvkit.json: expected workflow='cnvkit', got '{data['workflow']}'"
        )

    def test_test_size_key_present(self):
        """params-cnvkit.json must include 'test_size' for performance tuning."""
        data = self._load_params()
        assert 'test_size' in data, (
            "params-cnvkit.json is missing 'test_size' (used by is_large_run())"
        )
