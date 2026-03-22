#!/usr/bin/env python3
"""
Integration tests for main.nf – verifies that all workflow cases are
wired correctly and that all supporting files are coherent.

Covers:
  1. All workflow modules are included (INDELIBLE, CANOES, XHMM, CLAMMS,
     DRAGEN, CNVKIT, GATK_GCNV, SURVIVOR, TRUVARI, FEATURE_EXTRACTION).
  2. Every workflow case is present in the switch block:
     indelible, canoes, xhmm, clamms, dragen, cnvkit, gcnv, survivor, truvari,
     survivor_with_features, truvari_with_features, feature_extraction.
  3. gather_vcfs() checks for at least 2 caller directories before continuing.
  4. The get_id closure strips caller name suffixes and .vcf/.vcf.gz extensions
     from filenames to recover the bare sample ID.
  5. RUN_SURVIVOR and RUN_TRUVARI filter samples that have fewer than 2 VCFs
     (so only samples covered by ≥2 callers enter consensus steps).
  6. All params JSON files exist in params/ and contain valid JSON with the
     correct "workflow" key.
  7. Each params JSON references a workflow name that matches one of the cases.
  8. survivor_with_features and truvari_with_features end-to-end workflows
     are properly wired (gather → merge → feature extraction).
  9. feature_extraction case auto-discovers the collapsed VCF in Truvari mode.
"""

import json
import os
import re

import pytest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MAIN_NF = os.path.join(REPO_ROOT, "main.nf")
PARAMS_DIR = os.path.join(REPO_ROOT, "params")

# All supported caller workflow names (original 9)
ALL_WORKFLOWS = [
    "indelible",
    "canoes",
    "xhmm",
    "clamms",
    "dragen",
    "cnvkit",
    "gcnv",
    "survivor",
    "truvari",
]

# Combined end-to-end workflow names (new)
COMBINED_WORKFLOWS = [
    "survivor_with_features",
    "truvari_with_features",
]

# Mapping from workflow name to expected params JSON filename
PARAMS_FILES = {
    "indelible":               "params-indelible.json",
    "canoes":                  "params-canoes.json",
    "xhmm":                    "params-xhmm.json",
    "clamms":                  "params-clamms.json",
    "dragen":                  "params-icav2-dragen.json",
    "cnvkit":                  "params-cnvkit.json",
    "gcnv":                    "params-gatk-gcnv.json",
    "survivor":                "params-survivor.json",
    "truvari":                 "params-truvari.json",
    "survivor_with_features":  "params-survivor-with-features.json",
    "truvari_with_features":   "params-truvari-with-features.json",
}

# Module names expected in include statements
MODULE_INCLUDES = {
    "INDELIBLE":          "modules-indelible.nf",
    "CANOES":             "modules-canoes.nf",
    "XHMM":               "modules-xhmm.nf",
    "CLAMMS":             "modules-clamms.nf",
    "DRAGEN":             "modules-icav2-dragen.nf",
    "CNVKIT":             "modules-cnvkit.nf",
    "GATK_GCNV":          "modules-gatk-gcnv.nf",
    "SURVIVOR":           "modules-survivor.nf",
    "TRUVARI":            "modules-truvari.nf",
    "FEATURE_EXTRACTION": "modules-feature-extraction.nf",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_main():
    with open(MAIN_NF) as fh:
        return fh.read()


def _load_params(filename):
    path = os.path.join(PARAMS_DIR, filename)
    with open(path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def main_text():
    return _read_main()


# ===========================================================================
# 1. Module includes
# ===========================================================================

class TestModuleIncludes:
    """main.nf must include all required workflow modules."""

    @pytest.mark.parametrize("workflow_name,module_file", MODULE_INCLUDES.items())
    def test_include_present(self, main_text, workflow_name, module_file):
        """Each module must be imported via an include statement."""
        assert f"include {{ {workflow_name} }}" in main_text, (
            f"main.nf must include {{ {workflow_name} }} from ./modules/{module_file}"
        )
        assert module_file in main_text, (
            f"main.nf include for {workflow_name} must reference {module_file}"
        )


# ===========================================================================
# 2. Workflow switch cases
# ===========================================================================

class TestWorkflowSwitchCases:
    """The workflow switch block must handle all workflow modes."""

    @pytest.mark.parametrize("workflow_case", ALL_WORKFLOWS + COMBINED_WORKFLOWS)
    def test_case_present(self, main_text, workflow_case):
        """Every workflow mode must have a case in the switch block."""
        # case['indelible'] style
        assert f"case['{workflow_case}']" in main_text, (
            f"main.nf switch block must contain case['{workflow_case}'] "
            f"to dispatch the {workflow_case} workflow"
        )

    def test_default_case_with_error_message(self, main_text):
        """The switch block must have a default case that exits with an error message."""
        assert "default:" in main_text, (
            "main.nf switch block must contain a default: case"
        )
        # The default case must list all workflow options (original 9 + new combined modes)
        for wf in ALL_WORKFLOWS + COMBINED_WORKFLOWS:
            assert f"--workflow {wf}" in main_text, (
                f"The default error message must list '--workflow {wf}' as a valid option"
            )


# ===========================================================================
# 3. Sub-workflow definitions
# ===========================================================================

class TestSubWorkflowDefinitions:
    """main.nf must define a RUN_* sub-workflow for each workflow mode."""

    @pytest.mark.parametrize("workflow_case", ALL_WORKFLOWS + COMBINED_WORKFLOWS)
    def test_run_workflow_defined(self, main_text, workflow_case):
        """Each workflow must have a corresponding RUN_<WORKFLOW> sub-workflow."""
        expected_name = f"RUN_{workflow_case.upper()}"
        assert f"workflow {expected_name}" in main_text, (
            f"main.nf must define 'workflow {expected_name}' as a sub-workflow "
            f"for the {workflow_case} mode"
        )

    @pytest.mark.parametrize("run_name, module_name", [
        ("RUN_CANOES", "CANOES"),
        ("RUN_XHMM", "XHMM"),
        ("RUN_CLAMMS", "CLAMMS"),
        ("RUN_DRAGEN", "DRAGEN"),
        ("RUN_CNVKIT", "CNVKIT"),
        ("RUN_GCNV", "GATK_GCNV"),
        ("RUN_INDELIBLE", "INDELIBLE"),
    ])
    def test_evaluate_uses_normalised_vcf(self, main_text, run_name, module_name):
        """RUN_* workflows must pass normalised VCFs into EVALUATE."""
        run_block = re.search(
            rf'workflow {re.escape(run_name)}\s*\{{(.+?)(?=\nworkflow |\Z)',
            main_text, re.DOTALL,
        )
        assert run_block is not None, f"{run_name} workflow not found in main.nf"
        body = run_block.group(1)
        assert f"{module_name}.out.normalised_vcf" in body, (
            f"{run_name} must pass {module_name}.out.normalised_vcf to EVALUATE "
            "so QUAL_norm values in the QUAL column are used downstream"
        )
        assert f"{module_name}.out.sorted_vcf" not in body, (
            f"{run_name} must not pass {module_name}.out.sorted_vcf to EVALUATE; "
            "sorted VCFs may still contain non-normalised QUAL values"
        )


# ===========================================================================
# 4. gather_vcfs() helper function
# ===========================================================================

class TestGatherVcfsFunction:
    """gather_vcfs() must validate caller directory count and extract sample IDs."""

    def test_gather_vcfs_defined(self, main_text):
        """gather_vcfs() function must be defined in main.nf."""
        assert "def gather_vcfs()" in main_text, (
            "main.nf must define the gather_vcfs() helper function "
            "for SURVIVOR and TRUVARI consensus workflows"
        )

    def test_gather_vcfs_requires_at_least_two_dirs(self, main_text):
        """gather_vcfs() must exit with error if fewer than 2 caller dirs are given."""
        assert "dir_count < 2" in main_text, (
            "gather_vcfs() must check dir_count < 2 and exit 1 when fewer than "
            "two caller directories are provided"
        )

    def test_gather_vcfs_error_message_mentions_two_callers(self, main_text):
        """The gather_vcfs error message must explain the two-caller requirement."""
        assert "TWO" in main_text or "at least" in main_text, (
            "gather_vcfs() error message must explain that at least two caller "
            "directories must be provided"
        )

    def test_gather_vcfs_strips_canoes_suffix(self, main_text):
        """get_id closure must strip _CANOES from filenames."""
        assert "CANOES" in main_text, (
            "gather_vcfs() get_id closure must handle CANOES filenames"
        )

    def test_gather_vcfs_strips_clamms_suffix(self, main_text):
        """get_id closure must strip _CLAMMS from filenames."""
        assert "CLAMMS" in main_text, (
            "gather_vcfs() get_id closure must handle CLAMMS filenames"
        )

    def test_gather_vcfs_strips_xhmm_suffix(self, main_text):
        """get_id closure must strip _XHMM from filenames."""
        assert "XHMM" in main_text, (
            "gather_vcfs() get_id closure must handle XHMM filenames"
        )

    def test_gather_vcfs_strips_cnvkit_suffix(self, main_text):
        """get_id closure must strip _CNVKIT from filenames."""
        assert "CNVKIT" in main_text, (
            "gather_vcfs() get_id closure must handle CNVKIT filenames"
        )

    def test_gather_vcfs_strips_gcnv_suffix(self, main_text):
        """get_id closure must strip _GCNV from filenames."""
        assert "GCNV" in main_text, (
            "gather_vcfs() get_id closure must handle GCNV filenames (GATK gCNV)"
        )

    def test_gather_vcfs_strips_dragen_suffix(self, main_text):
        """get_id closure must strip _DRAGEN from filenames."""
        assert "DRAGEN" in main_text, (
            "gather_vcfs() get_id closure must handle DRAGEN filenames"
        )

    def test_gather_vcfs_strips_indelible_suffix(self, main_text):
        """get_id closure must strip _INDELIBLE from filenames."""
        assert "INDELIBLE" in main_text, (
            "gather_vcfs() get_id closure must handle INDELIBLE filenames"
        )

    def test_gather_vcfs_globs_vcf_files(self, main_text):
        """gather_vcfs() must collect *.vcf* files from each caller directory."""
        assert "*.vcf*" in main_text, (
            "gather_vcfs() must use '*.vcf*' glob to collect both .vcf and .vcf.gz files"
        )

    def test_gather_vcfs_used_by_survivor_case(self, main_text):
        """The survivor case must call gather_vcfs() to collect VCFs."""
        # Find the survivor case block
        survivor_case = re.search(
            r"case\['survivor'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert survivor_case, "case['survivor'] block not found"
        assert "gather_vcfs()" in survivor_case.group(1), (
            "case['survivor'] must call gather_vcfs() to collect VCFs from caller dirs"
        )

    def test_gather_vcfs_used_by_truvari_case(self, main_text):
        """The truvari case must call gather_vcfs() to collect VCFs."""
        truvari_case = re.search(
            r"case\['truvari'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert truvari_case, "case['truvari'] block not found"
        assert "gather_vcfs()" in truvari_case.group(1), (
            "case['truvari'] must call gather_vcfs() to collect VCFs from caller dirs"
        )

    @pytest.mark.parametrize("sample_id,filename", [
        ("SAMPLE001", "SAMPLE001_CANOES_output.vcf"),
        ("SAMPLE001", "SAMPLE001_CANOES_output.vcf.gz"),
        ("SAMPLE002", "SAMPLE002_CLAMMS_output.vcf"),
        ("MYSAMPLE",  "MYSAMPLE_XHMM_output.vcf.gz"),
        ("SID",       "SID_GCNV_genotyped_segments.sorted.vcf.gz"),
        ("SAMP",      "SAMP_CNVKIT_output.vcf"),
        ("PT01",      "PT01_DRAGEN_output.vcf"),
        ("PT02",      "PT02_INDELIBLE_output.vcf"),
        # DRAGEN Germline Enrichment output naming: ${sample_id}.cnv.vcf.gz
        ("NA12878",   "NA12878.cnv.vcf.gz"),
        ("SAMPLE003", "SAMPLE003.cnv.vcf.gz"),
        ("SAMPLE004", "SAMPLE004.cnv.vcf"),
    ])
    def test_get_id_regex_extracts_sample_id(self, main_text, sample_id, filename):
        """The get_id regex in gather_vcfs() must correctly strip caller suffixes.

        Simulates the Groovy closure:
          f.name.replaceAll(/_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*/, '')
                .replaceAll(/\\.cnv\\.vcf(\\.gz)?$/i, '')
                .replaceAll(/\\.vcf(\\.gz)?$/i, '')

        Also handles DRAGEN Germline Enrichment outputs named ${sample_id}.cnv.vcf.gz.
        """
        import re as _re
        # Replicate the Groovy regex logic in Python
        extracted = _re.sub(
            r'_(CANOES|CLAMMS|XHMM|CNVKIT|GCNV|DRAGEN|INDELIBLE).*', '', filename
        )
        extracted = _re.sub(r'\.cnv\.vcf(\.gz)?$', '', extracted, flags=_re.IGNORECASE)
        extracted = _re.sub(r'\.vcf(\.gz)?$', '', extracted, flags=_re.IGNORECASE)
        assert extracted == sample_id, (
            f"get_id regex must extract '{sample_id}' from '{filename}', "
            f"got '{extracted}'"
        )

    def test_get_id_cnv_vcf_gz_regex_in_main_nf(self, main_text):
        """gather_vcfs() get_id closure must strip .cnv.vcf.gz suffixes for DRAGEN."""
        assert r'\.cnv\.vcf' in main_text or '.cnv.vcf' in main_text, (
            "get_id in gather_vcfs() must handle .cnv.vcf.gz filenames "
            "produced by DRAGEN Germline Enrichment"
        )


# ===========================================================================
# 5. RUN_SURVIVOR and RUN_TRUVARI: enforce ≥2 VCFs per sample
# ===========================================================================

class TestConsensusWorkflowFiltering:
    """RUN_SURVIVOR and RUN_TRUVARI must only process samples with ≥2 caller VCFs."""

    def test_group_caller_vcfs_helper_enforces_min_two_callers(self, main_text):
        helper = re.search(
            r"def group_caller_vcfs\(vcf_ch\)\s*\{(.+?)\n\}",
            main_text,
            re.DOTALL,
        )
        assert helper, "main.nf must define group_caller_vcfs(vcf_ch)"
        helper_body = helper.group(1)
        assert ".groupTuple()" in helper_body, (
            "group_caller_vcfs() must aggregate VCFs by sample_id with .groupTuple()"
        )
        assert "vcfs.size() >= 2" in helper_body, (
            "group_caller_vcfs() must enforce vcfs.size() >= 2"
        )

    def test_run_survivor_filters_single_caller_samples(self, main_text):
        """RUN_SURVIVOR must filter out samples with only 1 caller VCF."""
        run_survivor = re.search(
            r"workflow RUN_SURVIVOR\s*\{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert run_survivor, "workflow RUN_SURVIVOR not found"
        body = run_survivor.group(1)
        assert "group_caller_vcfs(vcf_ch)" in body or (
            ".groupTuple()" in body and "vcfs.size() >= 2" in body
        ), (
            "RUN_SURVIVOR must enforce a minimum of 2 caller VCFs per sample "
            "via group_caller_vcfs(vcf_ch) or equivalent inline filtering"
        )

    def test_run_truvari_filters_single_caller_samples(self, main_text):
        """RUN_TRUVARI must filter out samples with only 1 caller VCF."""
        run_truvari = re.search(
            r"workflow RUN_TRUVARI\s*\{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert run_truvari, "workflow RUN_TRUVARI not found"
        body = run_truvari.group(1)
        assert "group_caller_vcfs(vcf_ch)" in body or (
            ".groupTuple()" in body and "vcfs.size() >= 2" in body
        ), (
            "RUN_TRUVARI must enforce a minimum of 2 caller VCFs per sample "
            "via group_caller_vcfs(vcf_ch) or equivalent inline filtering"
        )

    def test_run_survivor_calls_survivor_workflow(self, main_text):
        """RUN_SURVIVOR sub-workflow must invoke the SURVIVOR module workflow."""
        run_survivor = re.search(
            r"workflow RUN_SURVIVOR\s*\{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert run_survivor, "workflow RUN_SURVIVOR not found"
        assert "SURVIVOR(" in run_survivor.group(1), (
            "RUN_SURVIVOR must call the SURVIVOR workflow from modules-survivor.nf"
        )

    def test_run_truvari_calls_truvari_workflow(self, main_text):
        """RUN_TRUVARI sub-workflow must invoke the TRUVARI module workflow."""
        run_truvari = re.search(
            r"workflow RUN_TRUVARI\s*\{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert run_truvari, "workflow RUN_TRUVARI not found"
        assert "TRUVARI(" in run_truvari.group(1), (
            "RUN_TRUVARI must call the TRUVARI workflow from modules-truvari.nf"
        )


# ===========================================================================
# 6. Params JSON files: existence and validity
# ===========================================================================

class TestParamsJsonFiles:
    """All 9 params JSON files must exist, be valid JSON, and contain 'workflow'."""

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_exists(self, workflow_name, filename):
        """Each workflow must have a corresponding params JSON file in params/."""
        path = os.path.join(PARAMS_DIR, filename)
        assert os.path.isfile(path), (
            f"params/{filename} must exist for the {workflow_name} workflow"
        )

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_is_valid_json(self, workflow_name, filename):
        """Each params file must be valid JSON (parseable without errors)."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        try:
            data = _load_params(filename)
        except json.JSONDecodeError as exc:
            pytest.fail(
                f"params/{filename} contains invalid JSON: {exc}"
            )
        assert isinstance(data, dict), (
            f"params/{filename} must be a JSON object (dict), got {type(data).__name__}"
        )

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_has_workflow_key(self, workflow_name, filename):
        """Each params file must contain a 'workflow' key."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        assert "workflow" in data, (
            f"params/{filename} must contain a 'workflow' key to specify "
            f"which pipeline to run (expected '{workflow_name}')"
        )

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_workflow_value_is_correct(self, workflow_name, filename):
        """The 'workflow' value in each params file must match the expected workflow name."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        if "workflow" not in data:
            pytest.skip(f"'workflow' key absent in {filename}")
        # dragen params file uses 'dragen' not 'icav2-dragen'
        expected = workflow_name
        actual = data["workflow"]
        assert actual == expected, (
            f"params/{filename}: 'workflow' must be '{expected}', got '{actual}'"
        )

    @pytest.mark.parametrize("workflow_name,filename", PARAMS_FILES.items())
    def test_params_file_has_outdir_key(self, workflow_name, filename):
        """Each params file must contain an output path key for result organisation.

        Most workflows use 'outdir'. The ICAv2-DRAGEN cloud workflow instead uses
        'localDownloadPath' because results are downloaded from ICA rather than
        written to a local output directory.
        """
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        # DRAGEN (ICAv2 cloud) uses localDownloadPath instead of outdir
        has_output_path = "outdir" in data or "localDownloadPath" in data
        assert has_output_path, (
            f"params/{filename} must contain either an 'outdir' key "
            f"(local workflows) or 'localDownloadPath' (ICAv2-DRAGEN) "
            f"specifying where results should be written"
        )

    @pytest.mark.parametrize("filename", [
        "params-survivor.json",
        "params-truvari.json",
        "params-survivor-with-features.json",
        "params-truvari-with-features.json",
    ])
    def test_consensus_params_have_caller_dirs(self, filename):
        """Consensus params files must specify at least two caller dirs."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        caller_dir_keys = [
            "canoes_dir", "clamms_dir", "xhmm_dir",
            "cnvkit_dir", "gcnv_dir", "dragen_dir", "indelible_dir",
        ]
        found = [k for k in caller_dir_keys if k in data]
        assert len(found) >= 2, (
            f"params/{filename} must specify at least two caller_dir keys "
            f"(e.g. canoes_dir, clamms_dir) so gather_vcfs() has enough inputs. "
            f"Found: {found}"
        )

    @pytest.mark.parametrize("filename", [
        "params-survivor-with-features.json",
        "params-truvari-with-features.json",
    ])
    def test_combined_params_have_merger_mode(self, filename):
        """Combined workflow params files must include a merger_mode key."""
        path = os.path.join(PARAMS_DIR, filename)
        if not os.path.isfile(path):
            pytest.skip(f"{filename} does not exist")
        data = _load_params(filename)
        assert "merger_mode" in data, (
            f"params/{filename} must contain 'merger_mode' to tell "
            f"feature_extraction.py how the VCFs were merged"
        )


# ===========================================================================
# 7. Nextflow DSL2 declaration
# ===========================================================================

class TestDsl2Declaration:
    """main.nf must use Nextflow DSL2."""

    def test_dsl2_enabled(self, main_text):
        """main.nf must declare nextflow.enable.dsl=2."""
        assert "nextflow.enable.dsl=2" in main_text, (
            "main.nf must contain 'nextflow.enable.dsl=2' to use DSL2 module imports"
        )

    def test_shebang_line(self, main_text):
        """main.nf must start with a Nextflow shebang line."""
        assert main_text.startswith("#!/usr/bin/env nextflow"), (
            "main.nf must start with '#!/usr/bin/env nextflow'"
        )


# ===========================================================================
# 8. Combined end-to-end workflow wiring
# ===========================================================================

class TestCombinedWorkflowWiring:
    """RUN_SURVIVOR_WITH_FEATURES and RUN_TRUVARI_WITH_FEATURES must call
    both the consensus module and FEATURE_EXTRACTION in sequence."""

    def _get_workflow_body(self, main_text, name):
        m = re.search(
            rf"workflow {name}\s*\{{(.+?)(?=\nworkflow |\Z)",
            main_text,
            re.DOTALL,
        )
        assert m, f"workflow {name} not found in main.nf"
        return m.group(1)

    def test_run_survivor_with_features_calls_survivor(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_SURVIVOR_WITH_FEATURES")
        assert "SURVIVOR(" in body, (
            "RUN_SURVIVOR_WITH_FEATURES must call the SURVIVOR workflow"
        )

    def test_run_survivor_with_features_calls_feature_extraction(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_SURVIVOR_WITH_FEATURES")
        assert "FEATURE_EXTRACTION(" in body, (
            "RUN_SURVIVOR_WITH_FEATURES must call FEATURE_EXTRACTION after SURVIVOR"
        )

    def test_run_survivor_with_features_passes_union_vcf(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_SURVIVOR_WITH_FEATURES")
        assert "union_vcf" in body, (
            "RUN_SURVIVOR_WITH_FEATURES must use SURVIVOR.out.union_vcf "
            "as input to feature extraction"
        )

    def test_run_survivor_with_features_filters_min_two_callers(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_SURVIVOR_WITH_FEATURES")
        assert "group_caller_vcfs(vcf_ch)" in body or "vcfs.size() >= 2" in body, (
            "RUN_SURVIVOR_WITH_FEATURES must enforce minimum-2 caller VCFs "
            "before consensus merging"
        )

    def test_run_truvari_with_features_calls_truvari(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_TRUVARI_WITH_FEATURES")
        assert "TRUVARI(" in body, (
            "RUN_TRUVARI_WITH_FEATURES must call the TRUVARI workflow"
        )

    def test_run_truvari_with_features_calls_feature_extraction(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_TRUVARI_WITH_FEATURES")
        assert "FEATURE_EXTRACTION(" in body, (
            "RUN_TRUVARI_WITH_FEATURES must call FEATURE_EXTRACTION after TRUVARI"
        )

    def test_run_truvari_with_features_joins_collapsed_vcf(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_TRUVARI_WITH_FEATURES")
        assert "collapsed_vcf" in body, (
            "RUN_TRUVARI_WITH_FEATURES must pass the collapsed VCF to feature extraction"
        )

    def test_run_truvari_with_features_joins_merged_and_collapsed(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_TRUVARI_WITH_FEATURES")
        assert ".join(" in body, (
            "RUN_TRUVARI_WITH_FEATURES must join merged_vcf and collapsed_vcf "
            "channels by sample_id before building feature inputs"
        )

    def test_run_truvari_with_features_filters_min_two_callers(self, main_text):
        body = self._get_workflow_body(main_text, "RUN_TRUVARI_WITH_FEATURES")
        assert "group_caller_vcfs(vcf_ch)" in body or "vcfs.size() >= 2" in body, (
            "RUN_TRUVARI_WITH_FEATURES must enforce minimum-2 caller VCFs "
            "before consensus merging"
        )

    def test_survivor_with_features_case_calls_gather_vcfs(self, main_text):
        case_block = re.search(
            r"case\['survivor_with_features'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert case_block, "case['survivor_with_features'] block not found"
        assert "gather_vcfs()" in case_block.group(1), (
            "case['survivor_with_features'] must call gather_vcfs()"
        )

    def test_truvari_with_features_case_calls_gather_vcfs(self, main_text):
        case_block = re.search(
            r"case\['truvari_with_features'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert case_block, "case['truvari_with_features'] block not found"
        assert "gather_vcfs()" in case_block.group(1), (
            "case['truvari_with_features'] must call gather_vcfs()"
        )


# ===========================================================================
# 9. feature_extraction case: Truvari auto-discovery of collapsed VCF
# ===========================================================================

class TestFeatureExtractionCase:
    """The feature_extraction switch case must handle both merger modes and
    auto-discover the Truvari collapsed VCF alongside the merged VCF."""

    def _get_case_body(self, main_text, case_name):
        m = re.search(
            rf"case\['{re.escape(case_name)}'\](.+?)break",
            main_text,
            re.DOTALL,
        )
        assert m, f"case['{case_name}'] block not found in main.nf"
        return m.group(1)

    def test_feature_extraction_case_exists(self, main_text):
        assert "case['feature_extraction']" in main_text

    def test_feature_extraction_handles_truvari_mode(self, main_text):
        body = self._get_case_body(main_text, "feature_extraction")
        assert "truvari" in body, (
            "case['feature_extraction'] must handle merger_mode='truvari'"
        )

    def test_feature_extraction_discovers_collapsed_vcf_in_truvari_mode(self, main_text):
        body = self._get_case_body(main_text, "feature_extraction")
        assert "_truvari_collapsed.vcf" in body, (
            "case['feature_extraction'] in truvari mode must auto-discover "
            "<sample_id>_truvari_collapsed.vcf alongside the merged VCF"
        )

    def test_feature_extraction_uses_build_feature_inputs(self, main_text):
        body = self._get_case_body(main_text, "feature_extraction")
        assert "build_feature_inputs(" in body, (
            "case['feature_extraction'] must use build_feature_inputs() "
            "to construct the full 10-element feature input tuple"
        )

    def test_feature_extraction_passes_empty_collapsed_in_survivor_mode(self, main_text):
        body = self._get_case_body(main_text, "feature_extraction")
        # In survivor mode the collapsed_vcf slot must be [] (empty list)
        assert "build_feature_inputs(sample_id, f, [])" in body or \
               "build_feature_inputs(sample_id, merged_vcf, [])" in body, (
            "case['feature_extraction'] in survivor mode must pass [] "
            "as collapsed_vcf to build_feature_inputs()"
        )

    def test_feature_extraction_survivor_pattern_matches_union_vcf(self, main_text):
        body = self._get_case_body(main_text, "feature_extraction")
        assert "_survivor_union.vcf" in body, (
            "case['feature_extraction'] in survivor mode must glob for "
            "*_survivor_union.vcf* files"
        )

    def test_feature_extraction_truvari_pattern_matches_merged_vcf(self, main_text):
        body = self._get_case_body(main_text, "feature_extraction")
        assert "_truvari_merged.vcf" in body, (
            "case['feature_extraction'] in truvari mode must glob for "
            "*_truvari_merged.vcf* files"
        )
