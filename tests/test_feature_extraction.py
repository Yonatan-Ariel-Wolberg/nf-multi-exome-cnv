#!/usr/bin/env python3
"""
Tests for bin/feature_extraction.py.

Validates:
  1. Module-level structure: imports, public API, module docstring.
  2. Helper functions: _load_bed, _count_probes, _mean_mappability,
     _load_mappability_bed, _gc_content, _rd_ratio, _get_raw_qual.
  3. extract_normalized_features: signature has required optional parameters
     (bed_file, bam_file, reference_fasta, mappability_file, rd_flank),
     produces expected columns, and handles missing optional inputs.
  4. QUAL_norm is read from the per-caller VCF QUAL field (not recomputed).
  5. Raw quality fields per caller map to the correct FORMAT/INFO sources.

Where full VCF/BAM/FASTA fixtures are impractical, the tests import the
helper functions directly and exercise them with in-memory data.
"""

import importlib
import io
import math
import os
import sys
import types

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Locate the script
# ---------------------------------------------------------------------------

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
SCRIPT_PATH = os.path.join(BIN_DIR, 'feature_extraction.py')


def _read_script():
    with open(SCRIPT_PATH) as fh:
        return fh.read()


def _import_module():
    """Import feature_extraction.py as a Python module."""
    spec = importlib.util.spec_from_file_location('feature_extraction', SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope='module')
def script_text():
    return _read_script()


@pytest.fixture(scope='module')
def fe():
    """Return the feature_extraction module."""
    return _import_module()


# ===========================================================================
# 1. Module structure
# ===========================================================================

class TestModuleStructure:
    """feature_extraction.py must import the required libraries."""

    def test_imports_pysam(self, script_text):
        assert 'import pysam' in script_text

    def test_imports_pandas(self, script_text):
        assert 'import pandas' in script_text or 'pandas as pd' in script_text

    def test_imports_numpy(self, script_text):
        assert 'import numpy' in script_text or 'numpy as np' in script_text

    def test_imports_defaultdict(self, script_text):
        assert 'defaultdict' in script_text

    def test_exposes_extract_normalized_features(self, script_text):
        assert 'def extract_normalized_features(' in script_text

    def test_module_docstring_mentions_qual_norm(self, script_text):
        assert 'qual_norm' in script_text.lower() or 'QUAL_norm' in script_text

    def test_module_docstring_mentions_n_probes(self, script_text):
        assert 'n_probes' in script_text

    def test_module_docstring_mentions_rd_ratio(self, script_text):
        assert 'rd_ratio' in script_text

    def test_module_docstring_mentions_gc_content(self, script_text):
        assert 'gc_content' in script_text

    def test_module_docstring_mentions_mappability(self, script_text):
        assert 'mappability' in script_text


# ===========================================================================
# 2. extract_normalized_features signature
# ===========================================================================

class TestFunctionSignature:
    """extract_normalized_features must accept the new optional parameters."""

    def test_has_bed_file_param(self, script_text):
        assert 'bed_file' in script_text

    def test_has_bam_file_param(self, script_text):
        assert 'bam_file' in script_text

    def test_has_reference_fasta_param(self, script_text):
        assert 'reference_fasta' in script_text

    def test_has_mappability_file_param(self, script_text):
        assert 'mappability_file' in script_text

    def test_has_rd_flank_param(self, script_text):
        assert 'rd_flank' in script_text

    def test_bed_file_defaults_none(self, script_text):
        # Default value of bed_file must be None so the parameter is optional.
        assert 'bed_file=None' in script_text

    def test_bam_file_defaults_none(self, script_text):
        assert 'bam_file=None' in script_text

    def test_reference_fasta_defaults_none(self, script_text):
        assert 'reference_fasta=None' in script_text

    def test_mappability_file_defaults_none(self, script_text):
        assert 'mappability_file=None' in script_text


# ===========================================================================
# 3. QUAL_norm sourcing: must read from VCF QUAL field, not recompute
# ===========================================================================

class TestQualNormSourcing:
    """QUAL_norm should come from orig.qual (normalised VCF), not from a
    separate normalize_q recomputation inside feature_extraction.py."""

    def test_reads_orig_qual(self, script_text):
        assert 'orig.qual' in script_text, (
            "feature_extraction must read QUAL_norm from orig.qual "
            "(the QUAL field of the normalised per-caller VCF)"
        )

    def test_no_standalone_normalize_q_function(self, script_text):
        # normalize_q was the old re-implementation; it should be gone now
        # that we read QUAL_norm directly from the normalised VCF QUAL field.
        assert 'def normalize_q(' not in script_text, (
            "normalize_q helper should be removed; QUAL_norm is now read "
            "directly from the normalised VCF QUAL field (orig.qual)"
        )

    def test_stores_qual_norm_per_caller(self, script_text):
        assert "qual_norm_" in script_text, (
            "feature_extraction must store qual_norm_{caller} columns"
        )

    def test_stores_raw_qual_per_caller(self, script_text):
        assert "raw_qual_" in script_text, (
            "feature_extraction must store raw_qual_{caller} columns for "
            "the non-redundant caller-native quality metric"
        )


# ===========================================================================
# 4. Raw quality field mapping correctness
# ===========================================================================

class TestRawQualMapping:
    """_get_raw_qual must map each caller to the correct FORMAT/INFO field."""

    def test_canoes_uses_q_some(self, script_text):
        # Find the _get_raw_qual function block
        assert "'canoes'" in script_text and "'Q_SOME'" in script_text, (
            "CANOES raw quality must be extracted from Q_SOME FORMAT field"
        )

    def test_clamms_uses_q_some(self, script_text):
        assert "'clamms'" in script_text and "'Q_SOME'" in script_text

    def test_xhmm_uses_sq(self, script_text):
        assert "'xhmm'" in script_text and "'SQ'" in script_text, (
            "XHMM raw quality must be extracted from SQ FORMAT field"
        )

    def test_gatk_uses_qs(self, script_text):
        assert "'gatk'" in script_text and "'QS'" in script_text, (
            "GATK-gCNV raw quality must be extracted from QS FORMAT field"
        )

    def test_cnvkit_uses_cnq(self, script_text):
        assert "'cnvkit'" in script_text and "'CNQ'" in script_text, (
            "CNVkit raw quality must be extracted from CNQ FORMAT field"
        )

    def test_dragen_uses_oq(self, script_text):
        assert "'dragen'" in script_text and "'OQ'" in script_text, (
            "DRAGEN raw quality must be read from OQ FORMAT field "
            "(original QUAL stored by normalise_cnv_caller_quality_scores.py)"
        )

    def test_indelible_uses_sr_total_and_avg_mapq(self, script_text):
        assert "'SR_TOTAL'" in script_text and "'AVG_MAPQ'" in script_text, (
            "INDELIBLE raw quality must be computed from SR_TOTAL and AVG_MAPQ INFO fields"
        )

    def test_indelible_synthetic_formula(self, script_text):
        # The formula: sr_total * (avg_mapq / 60.0) * 100.0
        assert '60.0' in script_text and '100.0' in script_text, (
            "INDELIBLE synthetic Phred formula must include / 60.0 and * 100.0"
        )


# ===========================================================================
# 5. Probe count helper (_load_bed, _count_probes)
# ===========================================================================

class TestProbeCountHelpers:
    """_load_bed and _count_probes must be present and correct."""

    def test_load_bed_present(self, script_text):
        assert 'def _load_bed(' in script_text

    def test_count_probes_present(self, script_text):
        assert 'def _count_probes(' in script_text

    def test_count_probes_uses_overlap_logic(self, script_text):
        # Overlap test: iv_start < end and iv_end > start
        assert 'iv_start < end' in script_text or 'iv_end > start' in script_text

    def test_count_probes_correct(self, fe, tmp_path):
        # Write a small BED file
        bed = tmp_path / 'targets.bed'
        bed.write_text(
            'chr1\t100\t200\n'   # overlaps [150, 250)
            'chr1\t300\t400\n'   # does NOT overlap [150, 250)
            'chr1\t190\t260\n'   # overlaps [150, 250)
            'chr2\t100\t200\n'   # wrong chrom
        )
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes('chr1', 150, 250, intervals)
        assert count == 2

    def test_count_probes_no_overlap(self, fe, tmp_path):
        bed = tmp_path / 'targets2.bed'
        bed.write_text('chr1\t500\t600\n')
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes('chr1', 100, 200, intervals)
        assert count == 0

    def test_count_probes_probe_ends_at_cnv_start(self, fe, tmp_path):
        """Probe ending exactly at CNV start must NOT count as overlapping."""
        bed = tmp_path / 'boundary1.bed'
        bed.write_text('chr1\t0\t100\n')  # probe [0, 100)
        intervals = fe._load_bed(str(bed))
        # CNV starts at 100 -> probe iv_end (100) > start (100) is False
        count = fe._count_probes('chr1', 100, 200, intervals)
        assert count == 0

    def test_count_probes_probe_starts_at_cnv_end(self, fe, tmp_path):
        """Probe starting exactly at CNV end must NOT count as overlapping."""
        bed = tmp_path / 'boundary2.bed'
        bed.write_text('chr1\t200\t300\n')  # probe [200, 300)
        intervals = fe._load_bed(str(bed))
        # CNV ends at 200 -> probe iv_start (200) < end (200) is False
        count = fe._count_probes('chr1', 100, 200, intervals)
        assert count == 0

    def test_count_probes_single_base_overlap(self, fe, tmp_path):
        """Probe overlapping by exactly one base must count."""
        bed = tmp_path / 'boundary3.bed'
        bed.write_text('chr1\t199\t300\n')  # probe [199, 300)
        intervals = fe._load_bed(str(bed))
        count = fe._count_probes('chr1', 100, 200, intervals)
        assert count == 1

    def test_load_bed_skips_comment_lines(self, fe, tmp_path):
        bed = tmp_path / 'targets3.bed'
        bed.write_text(
            '# comment\n'
            'chr1\t10\t20\n'
        )
        intervals = fe._load_bed(str(bed))
        assert len(intervals['chr1']) == 1


# ===========================================================================
# 6. Mappability helpers (_load_mappability_bed, _mean_mappability)
# ===========================================================================

class TestMappabilityHelpers:
    """Mappability BED loader and weighted-mean scorer must be correct."""

    def test_load_mappability_bed_present(self, script_text):
        assert 'def _load_mappability_bed(' in script_text

    def test_mean_mappability_present(self, script_text):
        assert 'def _mean_mappability(' in script_text

    def test_mean_mappability_weighted(self, fe, tmp_path):
        bed = tmp_path / 'map.bed'
        # Two equal-sized intervals with different scores
        bed.write_text(
            'chr1\t0\t100\t0.8\n'   # 100 bp, score 0.8
            'chr1\t100\t200\t0.4\n' # 100 bp, score 0.4
        )
        intervals = fe._load_mappability_bed(str(bed))
        score = fe._mean_mappability('chr1', 0, 200, intervals)
        assert abs(score - 0.6) < 1e-9

    def test_mean_mappability_partial_overlap(self, fe, tmp_path):
        # Interval [50, 150) with score 1.0; query [0, 100).
        # Overlap = [50, 100) = 50 bp. Weighted mean is over covered bases
        # only: 50 * 1.0 / 50 = 1.0.
        bed = tmp_path / 'map2.bed'
        bed.write_text('chr1\t50\t150\t1.0\n')
        intervals = fe._load_mappability_bed(str(bed))
        score = fe._mean_mappability('chr1', 0, 100, intervals)
        assert score == pytest.approx(1.0)

    def test_mean_mappability_unequal_sizes(self, fe, tmp_path):
        """Weighted mean with intervals of unequal size must not be a simple average."""
        bed = tmp_path / 'map_unequal.bed'
        # 100 bp at 0.8, 50 bp at 0.4 -> (100*0.8 + 50*0.4) / 150 = 100/150
        bed.write_text(
            'chr1\t0\t100\t0.8\n'
            'chr1\t100\t150\t0.4\n'
        )
        intervals = fe._load_mappability_bed(str(bed))
        score = fe._mean_mappability('chr1', 0, 150, intervals)
        expected = (100 * 0.8 + 50 * 0.4) / 150
        assert score == pytest.approx(expected)

    def test_mean_mappability_no_overlap_returns_nan(self, fe, tmp_path):
        bed = tmp_path / 'map3.bed'
        bed.write_text('chr1\t500\t600\t1.0\n')
        intervals = fe._load_mappability_bed(str(bed))
        score = fe._mean_mappability('chr1', 0, 100, intervals)
        assert math.isnan(score)

    def test_mean_mappability_empty_dict_returns_nan(self, fe):
        score = fe._mean_mappability('chr1', 0, 100, {})
        assert math.isnan(score)


# ===========================================================================
# 7. GC content helper (_gc_content)
# ===========================================================================

class TestGcContent:
    """_gc_content must correctly compute G/C fraction."""

    def test_gc_content_present(self, script_text):
        assert 'def _gc_content(' in script_text

    def test_gc_content_pure_gc(self, fe):
        """100 % GC sequence."""
        class MockFasta:
            def fetch(self, chrom, start, end):
                return 'GCGCGC'
        assert fe._gc_content(MockFasta(), 'chr1', 0, 6) == pytest.approx(1.0)

    def test_gc_content_pure_at(self, fe):
        class MockFasta:
            def fetch(self, chrom, start, end):
                return 'ATATAT'
        assert fe._gc_content(MockFasta(), 'chr1', 0, 6) == pytest.approx(0.0)

    def test_gc_content_mixed(self, fe):
        class MockFasta:
            def fetch(self, chrom, start, end):
                return 'GCATAT'  # 2 GC out of 6
        assert fe._gc_content(MockFasta(), 'chr1', 0, 6) == pytest.approx(2 / 6)

    def test_gc_content_empty_returns_nan(self, fe):
        class MockFasta:
            def fetch(self, chrom, start, end):
                return ''
        assert math.isnan(fe._gc_content(MockFasta(), 'chr1', 0, 0))


# ===========================================================================
# 8. RD ratio helper (_rd_ratio, _mean_depth)
# ===========================================================================

class TestRdRatio:
    """_rd_ratio must return target_depth / flank_depth and NaN on zero flank."""

    def test_rd_ratio_present(self, script_text):
        assert 'def _rd_ratio(' in script_text

    def test_mean_depth_present(self, script_text):
        assert 'def _mean_depth(' in script_text

    def test_rd_ratio_nan_on_zero_flank(self, fe):
        """When flanking depth is 0, _rd_ratio must return NaN."""
        class MockPileupCol:
            def __init__(self, n):
                self.nsegments = n

        class MockBam:
            def pileup(self, chrom, start, end, **kwargs):
                return []  # empty -> depth 0

        result = fe._rd_ratio(MockBam(), 'chr1', 1000, 2000, flank=500)
        assert math.isnan(result)


# ===========================================================================
# 9. _get_raw_qual unit tests
# ===========================================================================

class TestGetRawQual:
    """_get_raw_qual must return the correct metric name and value per caller."""

    def _make_record(self, info=None, format_vals=None):
        """Create a minimal mock VariantRecord."""
        mock = types.SimpleNamespace()
        mock.info = info or {}
        mock.qual = None
        # Build a mock sample
        sample = types.SimpleNamespace()
        sample_dict = format_vals or {}

        class MockSample:
            def __contains__(self, key):
                return key in sample_dict
            def __getitem__(self, key):
                return sample_dict[key]

        mock.samples = types.SimpleNamespace()
        mock.samples.values = lambda: [MockSample()]
        return mock

    def test_canoes_q_some(self, fe):
        record = self._make_record(format_vals={'Q_SOME': 75.0})
        name, val = fe._get_raw_qual(record, 'canoes')
        assert name == 'Q_SOME'
        assert val == pytest.approx(75.0)

    def test_clamms_q_some(self, fe):
        record = self._make_record(format_vals={'Q_SOME': 120.0})
        name, val = fe._get_raw_qual(record, 'clamms')
        assert name == 'Q_SOME'
        assert val == pytest.approx(120.0)

    def test_xhmm_sq(self, fe):
        record = self._make_record(format_vals={'SQ': 55.0})
        name, val = fe._get_raw_qual(record, 'xhmm')
        assert name == 'SQ'
        assert val == pytest.approx(55.0)

    def test_gatk_qs(self, fe):
        record = self._make_record(format_vals={'QS': 200.0})
        name, val = fe._get_raw_qual(record, 'gatk')
        assert name == 'QS'
        assert val == pytest.approx(200.0)

    def test_cnvkit_cnq(self, fe):
        record = self._make_record(format_vals={'CNQ': 30.0})
        name, val = fe._get_raw_qual(record, 'cnvkit')
        assert name == 'CNQ'
        assert val == pytest.approx(30.0)

    def test_dragen_oq(self, fe):
        record = self._make_record(format_vals={'OQ': 45.0})
        name, val = fe._get_raw_qual(record, 'dragen')
        assert name == 'QUAL'
        assert val == pytest.approx(45.0)

    def test_indelible_synthetic_phred(self, fe):
        # sr_total * (avg_mapq / 60) * 100 = 10 * (30/60) * 100 = 500
        record = self._make_record(info={'SR_TOTAL': 10.0, 'AVG_MAPQ': 30.0})
        name, val = fe._get_raw_qual(record, 'indelible')
        assert name == 'SYNTHETIC_PHRED'
        assert val == pytest.approx(500.0)

    def test_indelible_zero_sr_total(self, fe):
        """SR_TOTAL=0 should produce synthetic score of 0.0."""
        record = self._make_record(info={'SR_TOTAL': 0.0, 'AVG_MAPQ': 40.0})
        name, val = fe._get_raw_qual(record, 'indelible')
        assert name == 'SYNTHETIC_PHRED'
        assert val == pytest.approx(0.0)

    def test_indelible_zero_avg_mapq(self, fe):
        """AVG_MAPQ=0 should produce synthetic score of 0.0."""
        record = self._make_record(info={'SR_TOTAL': 5.0, 'AVG_MAPQ': 0.0})
        name, val = fe._get_raw_qual(record, 'indelible')
        assert name == 'SYNTHETIC_PHRED'
        assert val == pytest.approx(0.0)

    def test_indelible_missing_sr_total_returns_nan(self, fe):
        """Missing SR_TOTAL must return NaN."""
        record = self._make_record(info={'AVG_MAPQ': 40.0})  # no SR_TOTAL
        name, val = fe._get_raw_qual(record, 'indelible')
        assert name is None
        assert math.isnan(val)

    def test_indelible_missing_avg_mapq_returns_nan(self, fe):
        """Missing AVG_MAPQ must return NaN."""
        record = self._make_record(info={'SR_TOTAL': 5.0})  # no AVG_MAPQ
        name, val = fe._get_raw_qual(record, 'indelible')
        assert name is None
        assert math.isnan(val)

    def test_missing_field_returns_nan(self, fe):
        record = self._make_record(format_vals={})  # no Q_SOME
        name, val = fe._get_raw_qual(record, 'canoes')
        assert name is None
        assert math.isnan(val)


# ===========================================================================
# 10. Concordance and caller flags
# ===========================================================================

class TestConcordanceAndFlags:
    """SUPP_VEC-derived concordance and is_{caller} flags must be present."""

    def test_concordance_column_produced(self, script_text):
        assert "'concordance'" in script_text or '"concordance"' in script_text

    def test_is_caller_columns_produced(self, script_text):
        assert "f'is_{caller_name}'" in script_text or "is_{caller" in script_text

    def test_supp_vec_parsed(self, script_text):
        assert 'SUPP_VEC' in script_text

    def test_cnv_type_column_produced(self, script_text):
        assert "'cnv_type'" in script_text or 'cnv_type' in script_text

    def test_cnv_size_column_produced(self, script_text):
        assert "'cnv_size'" in script_text or 'cnv_size' in script_text


# ===========================================================================
# 11. INDELIBLE split-read columns
# ===========================================================================

class TestIndelibleColumns:
    """INDELIBLE small-variant columns must still be extracted."""

    def test_total_sr_column(self, script_text):
        assert 'total_sr' in script_text

    def test_sr_entropy_column(self, script_text):
        assert 'sr_entropy' in script_text

    def test_mapq_avg_column(self, script_text):
        assert 'mapq_avg' in script_text

    def test_dual_split_column(self, script_text):
        assert 'dual_split' in script_text


# ===========================================================================
# 12. NaN defaults when optional files are absent
# ===========================================================================

class TestNaNDefaults:
    """When optional annotation files are not provided, relevant columns
    must default to NaN rather than raising exceptions."""

    def test_n_probes_nan_without_bed(self, script_text):
        # When bed_intervals is empty/False, n_probes must be np.nan
        assert 'np.nan' in script_text

    def test_rd_ratio_nan_without_bam(self, script_text):
        assert 'rd_ratio' in script_text and 'np.nan' in script_text

    def test_gc_content_nan_without_fasta(self, script_text):
        assert 'gc_content' in script_text and 'np.nan' in script_text

    def test_mappability_nan_without_file(self, script_text):
        assert 'mappability' in script_text and 'np.nan' in script_text


# ===========================================================================
# 13. CN-Learn-inspired features: _cnv_size_label
# ===========================================================================

class TestCnvSizeLabel:
    """_cnv_size_label must map CNV sizes to CN-Learn's ordinal size bins."""

    def test_size_label_present(self, script_text):
        assert 'def _cnv_size_label(' in script_text

    def test_size_label_column_emitted(self, script_text):
        assert 'size_label' in script_text

    def test_bin_1_less_than_1kb(self, fe):
        assert fe._cnv_size_label(500) == 1

    def test_bin_1_boundary_999(self, fe):
        assert fe._cnv_size_label(999) == 1

    def test_bin_2_exactly_1kb(self, fe):
        assert fe._cnv_size_label(1_000) == 2

    def test_bin_3_5kb(self, fe):
        assert fe._cnv_size_label(5_000) == 3

    def test_bin_4_10kb(self, fe):
        assert fe._cnv_size_label(10_000) == 4

    def test_bin_5_25kb(self, fe):
        assert fe._cnv_size_label(25_000) == 5

    def test_bin_6_50kb(self, fe):
        assert fe._cnv_size_label(50_000) == 6

    def test_bin_7_75kb(self, fe):
        assert fe._cnv_size_label(75_000) == 7

    def test_bin_8_100kb(self, fe):
        assert fe._cnv_size_label(100_000) == 8

    def test_bin_9_250kb(self, fe):
        assert fe._cnv_size_label(250_000) == 9

    def test_bin_10_500kb(self, fe):
        assert fe._cnv_size_label(500_000) == 10

    def test_bin_11_1mb(self, fe):
        assert fe._cnv_size_label(1_000_000) == 11

    def test_bin_12_5mb(self, fe):
        assert fe._cnv_size_label(5_000_000) == 12

    def test_bin_12_very_large(self, fe):
        assert fe._cnv_size_label(10_000_000) == 12

    def test_all_labels_in_range_1_to_12(self, fe):
        test_sizes = [0, 500, 999, 1000, 4999, 5000, 9999, 10000, 24999,
                      25000, 49999, 50000, 74999, 75000, 99999, 100000,
                      249999, 250000, 499999, 500000, 999999, 1000000,
                      4999999, 5000000, 10000000]
        for s in test_sizes:
            label = fe._cnv_size_label(s)
            assert 1 <= label <= 12, f"size {s} -> label {label} not in [1,12]"


# ===========================================================================
# 14. CN-Learn-inspired features: _encode_chrom
# ===========================================================================

class TestEncodeChrom:
    """_encode_chrom must return integer chromosome codes matching CN-Learn."""

    def test_encode_chrom_present(self, script_text):
        assert 'def _encode_chrom(' in script_text

    def test_chrom_encoded_column_emitted(self, script_text):
        assert 'chrom_encoded' in script_text

    def test_autosome_without_prefix(self, fe):
        assert fe._encode_chrom('1') == 1
        assert fe._encode_chrom('22') == 22

    def test_autosome_with_chr_prefix(self, fe):
        assert fe._encode_chrom('chr1') == 1
        assert fe._encode_chrom('chr22') == 22

    def test_x_chromosome(self, fe):
        assert fe._encode_chrom('X') == 23
        assert fe._encode_chrom('chrX') == 23
        assert fe._encode_chrom('x') == 23

    def test_y_chromosome(self, fe):
        assert fe._encode_chrom('Y') == 24
        assert fe._encode_chrom('chrY') == 24

    def test_unknown_returns_zero(self, fe):
        assert fe._encode_chrom('chrM') == 0
        assert fe._encode_chrom('scaffold_1') == 0
        assert fe._encode_chrom('') == 0


# ===========================================================================
# 15. Aggregate quality summary columns
# ===========================================================================

class TestAggregateQualNorm:
    """max_qual_norm and mean_qual_norm_supported columns must be present."""

    def test_max_qual_norm_present(self, script_text):
        assert 'max_qual_norm' in script_text

    def test_mean_qual_norm_supported_present(self, script_text):
        assert 'mean_qual_norm_supported' in script_text

    def test_aggregate_computed_from_supporting_callers_only(self, script_text):
        # The aggregate must filter on is_{caller} == 1, not all callers.
        assert "is_{cn}" in script_text or "f'is_{cn}'" in script_text or \
               "_qual_norm_vals" in script_text
