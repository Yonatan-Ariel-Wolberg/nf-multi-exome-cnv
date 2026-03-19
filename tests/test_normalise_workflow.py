#!/usr/bin/env python3
"""
Tests for modules/modules-normalise.nf and its integration into main.nf.

Validates:
  1. The module file exists and has the correct DSL2 / process / workflow structure.
  2. NORMALISE_CNV_QUALITY_SCORES process inputs / outputs are correct.
  3. The script block calls normalise_cnv_caller_quality_scores.py with the
     correct arguments, bgzips, and tabix-indexes the output.
  4. main.nf includes the NORMALISE module, defines RUN_NORMALISE, and has a
     case['normalise'] block that reads --vcf_dir and --caller params.
  5. params/params-normalise.json exists and contains valid JSON with the
     correct "workflow" key.
"""

import json
import os
import re

import pytest

REPO_ROOT   = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
NF_PATH     = os.path.join(REPO_ROOT, 'modules', 'modules-normalise.nf')
MAIN_NF     = os.path.join(REPO_ROOT, 'main.nf')
PARAMS_FILE = os.path.join(REPO_ROOT, 'params', 'params-normalise.json')


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def nf_text():
    with open(NF_PATH) as fh:
        return fh.read()


@pytest.fixture(scope='module')
def main_text():
    with open(MAIN_NF) as fh:
        return fh.read()


# ===========================================================================
# 1. Module structure
# ===========================================================================

class TestModuleStructure:
    """modules-normalise.nf must be a well-formed DSL2 module."""

    def test_file_exists(self):
        assert os.path.isfile(NF_PATH), (
            'modules/modules-normalise.nf must exist'
        )

    def test_dsl2_enabled(self, nf_text):
        assert 'nextflow.enable.dsl=2' in nf_text

    def test_has_normalise_process(self, nf_text):
        assert 'process NORMALISE_CNV_QUALITY_SCORES' in nf_text, (
            'Module must define process NORMALISE_CNV_QUALITY_SCORES'
        )

    def test_has_normalise_workflow(self, nf_text):
        assert 'workflow NORMALISE' in nf_text, (
            'Module must define workflow NORMALISE'
        )

    def test_uses_pysam_label(self, nf_text):
        assert "label 'pysam'" in nf_text, (
            "Process must use label 'pysam' (same container as other normalise steps)"
        )

    def test_publishdir_present(self, nf_text):
        assert 'publishDir' in nf_text


# ===========================================================================
# 2. Process inputs / outputs
# ===========================================================================

class TestProcessIO:
    """NORMALISE_CNV_QUALITY_SCORES must accept and emit the correct channels."""

    def test_input_sample_name(self, nf_text):
        assert 'sample_name' in nf_text

    def test_input_vcf_path(self, nf_text):
        assert 'path(vcf)' in nf_text or 'path vcf' in nf_text

    def test_input_caller_val(self, nf_text):
        assert 'val(caller)' in nf_text or 'val caller' in nf_text

    def test_output_normalised_vcf(self, nf_text):
        assert 'normalised_vcf' in nf_text
        assert '.normalised.vcf.gz' in nf_text

    def test_output_normalised_vcf_index(self, nf_text):
        assert 'normalised_vcf_index' in nf_text
        assert '.normalised.vcf.gz.tbi' in nf_text

    def test_emit_normalised_vcf(self, nf_text):
        assert 'emit: normalised_vcf' in nf_text

    def test_emit_normalised_vcf_index(self, nf_text):
        assert 'emit: normalised_vcf_index' in nf_text


# ===========================================================================
# 3. Script block content
# ===========================================================================

class TestScriptBlock:
    """The script block must call the normalisation script, bgzip, and tabix."""

    def test_calls_normalise_script(self, nf_text):
        assert 'normalise_cnv_caller_quality_scores.py' in nf_text

    def test_passes_input_vcf(self, nf_text):
        assert '--input_vcf' in nf_text

    def test_passes_output_vcf(self, nf_text):
        assert '--output_vcf' in nf_text

    def test_passes_caller_arg(self, nf_text):
        assert '--caller' in nf_text

    def test_bgzips_output(self, nf_text):
        assert 'bgzip' in nf_text

    def test_tabix_indexes_output(self, nf_text):
        assert 'tabix' in nf_text


# ===========================================================================
# 4. Workflow NORMALISE take/emit
# ===========================================================================

class TestWorkflowNormalise:
    """The NORMALISE workflow must expose the expected take/main/emit structure."""

    def test_workflow_has_take(self, nf_text):
        workflow_block = re.search(
            r'workflow NORMALISE\s*\{(.+)',
            nf_text,
            re.DOTALL,
        )
        assert workflow_block, 'NORMALISE workflow block not found'
        block = workflow_block.group(1)
        assert 'take:' in block

    def test_workflow_has_main(self, nf_text):
        assert 'main:' in nf_text

    def test_workflow_calls_process(self, nf_text):
        workflow_block = re.search(
            r'workflow NORMALISE\s*\{(.+)',
            nf_text,
            re.DOTALL,
        )
        assert workflow_block, 'NORMALISE workflow block not found'
        block = workflow_block.group(1)
        assert 'NORMALISE_CNV_QUALITY_SCORES' in block

    def test_workflow_emits_normalised_vcf(self, nf_text):
        workflow_block = re.search(
            r'workflow NORMALISE\s*\{(.+)',
            nf_text,
            re.DOTALL,
        )
        assert workflow_block, 'NORMALISE workflow block not found'
        block = workflow_block.group(1)
        assert 'emit:' in block
        assert 'normalised_vcf' in block


# ===========================================================================
# 5. main.nf integration
# ===========================================================================

class TestMainNfIntegration:
    """main.nf must include and wire up the NORMALISE module correctly."""

    def test_includes_normalise_module(self, main_text):
        assert "include { NORMALISE }" in main_text or \
               "include { NORMALISE}" in main_text, (
            "main.nf must include { NORMALISE } from modules-normalise.nf"
        )
        assert 'modules-normalise.nf' in main_text, (
            "main.nf include statement must reference modules-normalise.nf"
        )

    def test_run_normalise_workflow_defined(self, main_text):
        assert 'workflow RUN_NORMALISE' in main_text, (
            "main.nf must define 'workflow RUN_NORMALISE'"
        )

    def test_normalise_case_in_switch(self, main_text):
        assert "case['normalise']" in main_text, (
            "main.nf switch block must contain case['normalise']"
        )

    def test_normalise_case_reads_vcf_dir(self, main_text):
        normalise_case = re.search(
            r"case\['normalise'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert normalise_case, "case['normalise'] block not found"
        block = normalise_case.group(1)
        assert 'vcf_dir' in block, (
            "case['normalise'] must read --vcf_dir param"
        )

    def test_normalise_case_reads_caller_param(self, main_text):
        normalise_case = re.search(
            r"case\['normalise'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert normalise_case, "case['normalise'] block not found"
        block = normalise_case.group(1)
        assert 'caller' in block, (
            "case['normalise'] must read --caller param"
        )

    def test_normalise_case_calls_run_normalise(self, main_text):
        normalise_case = re.search(
            r"case\['normalise'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert normalise_case, "case['normalise'] block not found"
        block = normalise_case.group(1)
        assert 'RUN_NORMALISE' in block, (
            "case['normalise'] must call RUN_NORMALISE"
        )

    def test_normalise_in_help_message(self, main_text):
        assert '--workflow normalise' in main_text, (
            "The default error message must list '--workflow normalise' as an option"
        )


# ===========================================================================
# 6. params/params-normalise.json
# ===========================================================================

class TestParamsFile:
    """params/params-normalise.json must exist and be valid JSON."""

    def test_params_file_exists(self):
        assert os.path.isfile(PARAMS_FILE), (
            'params/params-normalise.json must exist'
        )

    def test_params_file_valid_json(self):
        with open(PARAMS_FILE) as fh:
            data = json.load(fh)
        assert isinstance(data, dict), 'params-normalise.json must be a JSON object'

    def test_params_file_workflow_key(self):
        with open(PARAMS_FILE) as fh:
            data = json.load(fh)
        assert data.get('workflow') == 'normalise', (
            "params-normalise.json must contain \"workflow\": \"normalise\""
        )

    def test_params_file_vcf_dir_key(self):
        with open(PARAMS_FILE) as fh:
            data = json.load(fh)
        assert 'vcf_dir' in data, (
            'params-normalise.json must include a vcf_dir key'
        )

    def test_params_file_caller_key(self):
        with open(PARAMS_FILE) as fh:
            data = json.load(fh)
        assert 'caller' in data, (
            'params-normalise.json must include a caller key'
        )
