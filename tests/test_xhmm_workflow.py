#!/usr/bin/env python3
"""
Tests for the XHMM workflow pipeline.

Validates:
  1. The modules-xhmm.nf module contains all required processes in the correct
     execution order (GROUP_BAMS → GATK_DOC → COMBINE_DOC → CALC_GC_XHMM →
     FILTER_SAMPLES → RUN_PCA → NORMALISE_PCA → FILTER_ZSCORE → FILTER_RD →
     DISCOVER_CNVS → GENOTYPE_CNVS → SPLIT_VCF → FILTER_XHMM_CNVS →
     BGZIP_SORT_INDEX_VCF).
  2. SPLIT_VCF produces per-sample VCFs named "{sample}_DATA.vcf".
  3. FILTER_XHMM_CNVS applies EQ, SQ, and NDQ quality filters.
  4. BGZIP_SORT_INDEX_VCF annotates with TOOL=XHMM.
  5. The bin/xcnv_to_cnv script correctly converts XHMM .xcnv format to PLINK
     .cnv format for DEL, DUP, and mixed CNV types.
  6. params-xhmm.json contains the required keys for the workflow.
"""

import json
import os
import re
import subprocess
import textwrap

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
XHMM_MODULE   = os.path.join(REPO_ROOT, 'modules', 'modules-xhmm.nf')
XCNV_TO_CNV   = os.path.join(REPO_ROOT, 'bin', 'xcnv_to_cnv')
PARAMS_XHMM   = os.path.join(REPO_ROOT, 'params', 'params-xhmm.json')


def _load_module():
    with open(XHMM_MODULE) as f:
        return f.read()


def _extract_workflow_block(content):
    """Extract the text of the XHMM workflow block."""
    m = re.search(r'workflow XHMM \{', content)
    assert m is not None, "workflow XHMM not found in modules-xhmm.nf"
    start = m.end() - 1
    depth = 0
    for i, ch in enumerate(content[start:], start):
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                return content[start: i + 1]
    raise AssertionError("Unmatched braces for workflow XHMM")


def _extract_process_block(content, process_name):
    """Extract the text of a named process block."""
    m = re.search(rf'process {re.escape(process_name)} \{{', content)
    assert m is not None, f"process {process_name} not found in modules-xhmm.nf"
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

class TestXhmmModuleStructure:
    """The modules-xhmm.nf file must define all required processes."""

    _REQUIRED_PROCESSES = [
        'GROUP_BAMS',
        'GATK_DOC',
        'COMBINE_DOC',
        'CALC_GC_XHMM',
        'FILTER_SAMPLES',
        'RUN_PCA',
        'NORMALISE_PCA',
        'FILTER_ZSCORE',
        'FILTER_RD',
        'DISCOVER_CNVS',
        'GENOTYPE_CNVS',
        'SPLIT_VCF',
        'FILTER_XHMM_CNVS',
        'BGZIP_SORT_INDEX_VCF',
    ]

    def test_all_required_processes_present(self):
        """Every required process is defined in modules-xhmm.nf."""
        content = _load_module()
        for proc in self._REQUIRED_PROCESSES:
            assert f'process {proc}' in content, (
                f"process {proc} not found in modules-xhmm.nf"
            )

    def test_workflow_xhmm_defined(self):
        """The top-level workflow XHMM is defined."""
        content = _load_module()
        assert 'workflow XHMM {' in content, "workflow XHMM not found in modules-xhmm.nf"

    def test_workflow_step_order(self):
        """All 14 workflow steps appear in the correct order inside workflow XHMM."""
        content = _load_module()
        wf = _extract_workflow_block(content)

        ordered_steps = [
            'GROUP_BAMS',
            'GATK_DOC',
            'COMBINE_DOC',
            'CALC_GC_XHMM',
            'FILTER_SAMPLES',
            'RUN_PCA',
            'NORMALISE_PCA',
            'FILTER_ZSCORE',
            'FILTER_RD',
            'DISCOVER_CNVS',
            'GENOTYPE_CNVS',
            'SPLIT_VCF',
            'FILTER_XHMM_CNVS',
            'BGZIP_SORT_INDEX_VCF',
        ]
        positions = []
        for step in ordered_steps:
            pos = wf.find(step)
            assert pos != -1, f"{step} not found in workflow XHMM block"
            positions.append(pos)

        for i in range(len(positions) - 1):
            assert positions[i] < positions[i + 1], (
                f"{ordered_steps[i]} must appear before {ordered_steps[i + 1]} "
                f"in workflow XHMM"
            )

    def test_workflow_emits_sorted_vcf(self):
        """workflow XHMM emits sorted_vcf from BGZIP_SORT_INDEX_VCF."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        assert 'sorted_vcf' in wf, "workflow XHMM must emit sorted_vcf"
        assert 'BGZIP_SORT_INDEX_VCF' in wf, \
            "workflow XHMM must reference BGZIP_SORT_INDEX_VCF in emit"

    def test_combine_doc_precedes_filter_samples(self):
        """COMBINE_DOC must appear before FILTER_SAMPLES (normalisation start)."""
        content = _load_module()
        wf = _extract_workflow_block(content)
        combine_pos = wf.find('COMBINE_DOC')
        filter_pos = wf.find('FILTER_SAMPLES')
        assert combine_pos < filter_pos, (
            "COMBINE_DOC must precede FILTER_SAMPLES in workflow XHMM"
        )


# ===========================================================================
# Process-level details
# ===========================================================================

class TestXhmmProcessDetails:
    """Key process scripts contain the correct logic."""

    def test_split_vcf_output_pattern(self):
        """SPLIT_VCF output files follow the '{sample}_DATA.vcf' naming pattern."""
        content = _load_module()
        block = _extract_process_block(content, 'SPLIT_VCF')
        assert '_DATA.vcf' in block, (
            "SPLIT_VCF must produce per-sample files named '{sample}_DATA.vcf'"
        )

    def test_filter_xhmm_cnvs_uses_eq_sq_ndq(self):
        """FILTER_XHMM_CNVS must filter on EQ, SQ, and NDQ quality fields."""
        content = _load_module()
        block = _extract_process_block(content, 'FILTER_XHMM_CNVS')
        for field in ('EQ', 'SQ', 'NDQ'):
            assert field in block, (
                f"FILTER_XHMM_CNVS must apply a quality filter on FORMAT/{field}"
            )

    def test_filter_xhmm_cnvs_threshold_is_60(self):
        """FILTER_XHMM_CNVS uses a quality threshold of 60 for EQ/SQ/NDQ."""
        content = _load_module()
        block = _extract_process_block(content, 'FILTER_XHMM_CNVS')
        assert '60' in block, (
            "FILTER_XHMM_CNVS must use a quality score threshold of 60"
        )

    def test_bgzip_sort_index_annotates_tool_xhmm(self):
        """BGZIP_SORT_INDEX_VCF annotates each variant with TOOL=XHMM."""
        content = _load_module()
        block = _extract_process_block(content, 'BGZIP_SORT_INDEX_VCF')
        assert 'TOOL=XHMM' in block or '"XHMM"' in block or "'XHMM'" in block, (
            "BGZIP_SORT_INDEX_VCF must annotate variants with TOOL=XHMM"
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

    def test_group_bams_uses_batch_size_param(self):
        """GROUP_BAMS must use params.xhmm_batch_size to control batch size."""
        content = _load_module()
        block = _extract_process_block(content, 'GROUP_BAMS')
        assert 'xhmm_batch_size' in block, (
            "GROUP_BAMS must use params.xhmm_batch_size for the split batch size"
        )

    def test_combine_doc_produces_rd_txt(self):
        """COMBINE_DOC must produce DATA.RD.txt as output."""
        content = _load_module()
        block = _extract_process_block(content, 'COMBINE_DOC')
        assert 'DATA.RD.txt' in block, (
            "COMBINE_DOC must produce DATA.RD.txt (the combined read-depth matrix)"
        )


# ===========================================================================
# xcnv_to_cnv script
# ===========================================================================

class TestXcnvToCnvScript:
    """Unit tests for bin/xcnv_to_cnv – the XCNV-to-PLINK-CNV converter."""

    _XCNV_HEADER = "SAMPLE\tCNV\tINTERVAL\tKB\tCHR\tMID_BP\tTARGETS\tNUM_TARG\tQ_EXACT\tQ_SOME\tQ_NON_DIPLOID\tQ_START\tQ_STOP\tALL_Q\n"

    def _run_script(self, xcnv_content):
        """Write xcnv_content to a temp file, run xcnv_to_cnv, and return stdout."""
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xcnv', delete=False) as f:
            f.write(xcnv_content)
            path = f.name
        try:
            result = subprocess.run(
                ['bash', XCNV_TO_CNV, path],
                capture_output=True, text=True
            )
            return result.returncode, result.stdout, result.stderr
        finally:
            os.unlink(path)

    def _make_xcnv_row(self, sample, cnv_type, interval, num_targ=5, q_some=75.0):
        """Return a single xcnv data row (tab-separated)."""
        chr_part, coords = interval.split(':')
        bp1, bp2 = coords.split('-')
        kb = str(round((int(bp2) - int(bp1)) / 1000.0, 2))
        chr_num = chr_part.replace('chr', '')
        mid_bp = str((int(bp1) + int(bp2)) // 2)
        return (
            f"{sample}\t{cnv_type}\t{interval}\t{kb}\t{chr_num}\t{mid_bp}\t"
            f"target1;target2\t{num_targ}\t60\t{q_some}\t50\t55\t55\t70\n"
        )

    def test_script_exits_without_args(self):
        """xcnv_to_cnv prints usage and exits with status 1 when called without arguments."""
        result = subprocess.run(['bash', XCNV_TO_CNV], capture_output=True, text=True)
        assert result.returncode == 1, "Should exit with status 1 when no args given"
        assert 'xcnv' in result.stdout.lower() or 'usage' in result.stdout.lower(), (
            "Should print usage information when no args given"
        )

    def test_output_has_plink_header(self):
        """Output starts with the PLINK .cnv header (FID IID CHR BP1 BP2 TYPE SCORE SITE)."""
        content = self._XCNV_HEADER + self._make_xcnv_row(
            'SAMPLE1', 'DEL', 'chr1:1000-2000'
        )
        rc, stdout, _ = self._run_script(content)
        assert rc == 0, f"xcnv_to_cnv returned exit code {rc}"
        first_line = stdout.splitlines()[0]
        assert 'FID' in first_line and 'CHR' in first_line and 'TYPE' in first_line, (
            f"Output header does not match PLINK .cnv format: {first_line}"
        )

    def test_del_maps_to_type_1(self):
        """DEL CNVs must be encoded as TYPE=1 in the PLINK .cnv output."""
        content = self._XCNV_HEADER + self._make_xcnv_row(
            'SAMPLE1', 'DEL', 'chr1:1000-2000'
        )
        rc, stdout, _ = self._run_script(content)
        assert rc == 0
        data_lines = [l for l in stdout.splitlines() if not l.startswith('FID')]
        assert data_lines, "No data lines in output"
        fields = data_lines[0].split('\t')
        assert fields[5] == '1', (
            f"DEL should map to TYPE=1, got TYPE={fields[5]}"
        )

    def test_dup_maps_to_type_3(self):
        """DUP CNVs must be encoded as TYPE=3 in the PLINK .cnv output."""
        content = self._XCNV_HEADER + self._make_xcnv_row(
            'SAMPLE2', 'DUP', 'chr2:5000-8000'
        )
        rc, stdout, _ = self._run_script(content)
        assert rc == 0
        data_lines = [l for l in stdout.splitlines() if not l.startswith('FID')]
        assert data_lines, "No data lines in output"
        fields = data_lines[0].split('\t')
        assert fields[5] == '3', (
            f"DUP should map to TYPE=3, got TYPE={fields[5]}"
        )

    def test_sample_id_preserved(self):
        """The sample ID from the XCNV file is preserved in the FID column."""
        content = self._XCNV_HEADER + self._make_xcnv_row(
            'CHILD_01_1', 'DEL', 'chr1:1000-2000'
        )
        rc, stdout, _ = self._run_script(content)
        assert rc == 0
        data_lines = [l for l in stdout.splitlines() if not l.startswith('FID')]
        assert data_lines, "No data lines in output"
        fields = data_lines[0].split('\t')
        assert fields[0] == 'CHILD_01_1', (
            f"Sample ID not preserved in FID column: {fields[0]}"
        )

    def test_chromosome_prefix_stripped(self):
        """The 'chr' prefix is stripped from chromosome names in the output."""
        content = self._XCNV_HEADER + self._make_xcnv_row(
            'SAMPLE1', 'DEL', 'chr1:1000-2000'
        )
        rc, stdout, _ = self._run_script(content)
        assert rc == 0
        data_lines = [l for l in stdout.splitlines() if not l.startswith('FID')]
        assert data_lines, "No data lines in output"
        fields = data_lines[0].split('\t')
        assert fields[2] == '1', (
            f"'chr' prefix should be stripped; CHR column shows '{fields[2]}'"
        )

    def test_bp1_and_bp2_extracted_correctly(self):
        """BP1 and BP2 are correctly extracted from the INTERVAL column."""
        content = self._XCNV_HEADER + self._make_xcnv_row(
            'SAMPLE1', 'DUP', 'chr3:10000-20000'
        )
        rc, stdout, _ = self._run_script(content)
        assert rc == 0
        data_lines = [l for l in stdout.splitlines() if not l.startswith('FID')]
        assert data_lines, "No data lines in output"
        fields = data_lines[0].split('\t')
        assert fields[3] == '10000', f"BP1 incorrect: {fields[3]}"
        assert fields[4] == '20000', f"BP2 incorrect: {fields[4]}"

    def test_multiple_rows_all_converted(self):
        """All data rows in the XCNV file are present in the output."""
        content = (
            self._XCNV_HEADER
            + self._make_xcnv_row('SAMPLE1', 'DEL', 'chr1:1000-2000')
            + self._make_xcnv_row('SAMPLE2', 'DUP', 'chr2:5000-8000')
            + self._make_xcnv_row('SAMPLE3', 'DEL', 'chrX:100000-200000')
        )
        rc, stdout, _ = self._run_script(content)
        assert rc == 0
        data_lines = [l for l in stdout.splitlines() if not l.startswith('FID')]
        assert len(data_lines) == 3, (
            f"Expected 3 data rows in output, got {len(data_lines)}"
        )


# ===========================================================================
# params-xhmm.json
# ===========================================================================

class TestXhmmParamsJson:
    """params-xhmm.json must contain all keys required by the XHMM workflow."""

    _REQUIRED_KEYS = [
        'workflow',
        'outdir',
        'samplesheet_bams',
        'ref',
        'probes',
        'xhmm_conf',
        'xhmm_batch_size',
    ]

    def _load_params(self):
        with open(PARAMS_XHMM) as f:
            return json.load(f)

    def test_all_required_keys_present(self):
        """params-xhmm.json must contain all required keys."""
        data = self._load_params()
        for key in self._REQUIRED_KEYS:
            assert key in data, f"params-xhmm.json is missing required key '{key}'"

    def test_workflow_value_is_xhmm(self):
        """The 'workflow' key must be set to 'xhmm'."""
        data = self._load_params()
        assert data['workflow'] == 'xhmm', (
            f"params-xhmm.json: expected workflow='xhmm', got '{data['workflow']}'"
        )

    def test_xhmm_batch_size_is_positive(self):
        """The xhmm_batch_size must be a positive integer."""
        data = self._load_params()
        assert isinstance(data['xhmm_batch_size'], int) and data['xhmm_batch_size'] > 0, (
            f"params-xhmm.json: xhmm_batch_size must be a positive integer, "
            f"got {data['xhmm_batch_size']}"
        )
