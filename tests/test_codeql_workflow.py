#!/usr/bin/env python3
"""Tests for CodeQL workflow language/script analysis coverage."""

from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent
CODEQL_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "codeql.yml"


def _load_workflow():
    return yaml.safe_load(CODEQL_WORKFLOW.read_text(encoding="utf-8"))


def test_codeql_matrix_covers_python_and_c():
    workflow = _load_workflow()
    matrix_entries = workflow["jobs"]["analyze"]["strategy"]["matrix"]["include"]
    languages = {entry["language"] for entry in matrix_entries}
    assert "python" in languages
    assert "c-cpp" in languages


def test_codeql_has_script_analysis_job_for_nextflow_r_shell():
    workflow = _load_workflow()
    steps = workflow["jobs"]["analyze-scripts"]["steps"]
    step_names = [step.get("name") for step in steps]
    assert "Install ShellCheck" in step_names
    assert "Analyze Shell scripts" in step_names
    assert "Analyze R scripts" in step_names
    assert "Analyze Nextflow scripts" in step_names

    shell_step = next((step for step in steps if step.get("name") == "Analyze Shell scripts"), None)
    assert shell_step is not None
    assert "shellcheck" in shell_step["run"]

    r_step = next((step for step in steps if step.get("name") == "Analyze R scripts"), None)
    assert r_step is not None
    assert "Rscript -e" in r_step["run"]
    assert "parse(file = commandArgs(trailingOnly = TRUE)[1])" in r_step["run"]
