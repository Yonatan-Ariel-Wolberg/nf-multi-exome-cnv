#!/usr/bin/env python3
"""Tests for CodeQL workflow language/script analysis coverage."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CODEQL_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "codeql.yml"


def _read_workflow():
    return CODEQL_WORKFLOW.read_text(encoding="utf-8")


def test_codeql_matrix_covers_python_and_c():
    text = _read_workflow()
    assert "- language: python" in text
    assert "- language: c-cpp" in text


def test_codeql_has_script_analysis_job_for_nextflow_r_shell():
    text = _read_workflow()
    assert "analyze-scripts:" in text
    assert "Analyze Shell scripts" in text
    assert "shellcheck" in text
    assert "Analyze R scripts" in text
    assert "Rscript -e" in text
    assert "parse(file = commandArgs(trailingOnly = TRUE)[1])" in text
    assert "Analyze Nextflow scripts" in text
    assert "GroovyShell" in text
    assert "parse(new File(args[0]))" in text
