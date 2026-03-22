#!/usr/bin/env python3
"""
Tests for bin/train_xgboost.py.

Validates:
  1. SUPPORTED_CALLERS lists all seven pipeline callers.
  2. MIN_CALLERS_FOR_TRAINING equals 2 (concordance requires >= 2 callers).
  3. validate_min_callers raises ValueError for 0 or 1 caller and passes for >= 2.
  4. prepare_training_data uses named is_* columns when present (feature-
     extraction output) and validates minimum callers.
  5. prepare_training_data falls back to the legacy supp_vec path and
     validates minimum callers there too.
"""

import importlib.util
import os
import sys

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Locate and import the script
# ---------------------------------------------------------------------------

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
SCRIPT_PATH = os.path.join(BIN_DIR, 'train_xgboost.py')


def _read_script():
    with open(SCRIPT_PATH) as fh:
        return fh.read()


def _import_module():
    spec = importlib.util.spec_from_file_location('train_xgboost', SCRIPT_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope='module')
def script_text():
    return _read_script()


@pytest.fixture(scope='module')
def tx():
    """Return the train_xgboost module."""
    return _import_module()


# ===========================================================================
# 1. SUPPORTED_CALLERS constant
# ===========================================================================

class TestSupportedCallers:

    def test_supported_callers_present(self, script_text):
        assert 'SUPPORTED_CALLERS' in script_text

    def test_supported_callers_has_seven_entries(self, tx):
        assert len(tx.SUPPORTED_CALLERS) == 7

    def test_supported_callers_contains_canoes(self, tx):
        assert 'canoes' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_clamms(self, tx):
        assert 'clamms' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_xhmm(self, tx):
        assert 'xhmm' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_gatk_gcnv(self, tx):
        assert 'gatk_gcnv' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_cnvkit(self, tx):
        assert 'cnvkit' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_dragen(self, tx):
        assert 'dragen' in tx.SUPPORTED_CALLERS

    def test_supported_callers_contains_indelible(self, tx):
        assert 'indelible' in tx.SUPPORTED_CALLERS


# ===========================================================================
# 2. MIN_CALLERS_FOR_TRAINING constant
# ===========================================================================

class TestMinCallersConstant:

    def test_min_callers_present(self, script_text):
        assert 'MIN_CALLERS_FOR_TRAINING' in script_text

    def test_min_callers_equals_two(self, tx):
        assert tx.MIN_CALLERS_FOR_TRAINING == 2

    def test_min_callers_less_than_supported_callers(self, tx):
        assert tx.MIN_CALLERS_FOR_TRAINING < len(tx.SUPPORTED_CALLERS)


# ===========================================================================
# 3. validate_min_callers function
# ===========================================================================

class TestValidateMinCallers:

    def test_validate_min_callers_present(self, script_text):
        assert 'def validate_min_callers(' in script_text

    def test_zero_callers_raises(self, tx):
        with pytest.raises(ValueError, match='At least 2'):
            tx.validate_min_callers([])

    def test_one_caller_raises(self, tx):
        with pytest.raises(ValueError, match='At least 2'):
            tx.validate_min_callers(['is_canoes'])

    def test_two_callers_passes(self, tx):
        tx.validate_min_callers(['is_canoes', 'is_clamms'])  # must not raise

    def test_seven_callers_passes(self, tx):
        cols = [f'is_{c}' for c in tx.SUPPORTED_CALLERS]
        tx.validate_min_callers(cols)  # must not raise

    def test_error_message_mentions_supported_callers(self, tx):
        with pytest.raises(ValueError, match='Supported callers'):
            tx.validate_min_callers([])

    def test_error_message_mentions_found_count(self, tx):
        with pytest.raises(ValueError, match='Found 1'):
            tx.validate_min_callers(['is_canoes'])


# ===========================================================================
# 4. prepare_training_data – named is_* columns (feature-extraction output)
# ===========================================================================

class TestPrepareTrainingDataNamedColumns:

    def _make_df_named(self, callers=None):
        """Build a minimal feature-extraction-style DataFrame."""
        if callers is None:
            callers = ['canoes', 'clamms']
        data = {
            'concordance': [1, 2],
            'cnv_size': [5000, 10000],
        }
        for c in callers:
            data[f'is_{c}'] = [1, 0]
            data[f'qual_norm_{c}'] = [80.0, np.nan]
        return pd.DataFrame(data)

    def test_returns_dataframe(self, tx):
        df = self._make_df_named()
        result = tx.prepare_training_data(df)
        assert isinstance(result, pd.DataFrame)

    def test_is_columns_preserved(self, tx):
        df = self._make_df_named(['canoes', 'clamms', 'xhmm'])
        result = tx.prepare_training_data(df)
        assert 'is_canoes' in result.columns
        assert 'is_clamms' in result.columns
        assert 'is_xhmm' in result.columns

    def test_one_caller_raises(self, tx):
        df = self._make_df_named(['canoes'])
        with pytest.raises(ValueError):
            tx.prepare_training_data(df)

    def test_no_supp_vec_column_needed(self, tx):
        """Feature-extraction output has no supp_vec; must not raise KeyError."""
        df = self._make_df_named(['canoes', 'clamms'])
        assert 'supp_vec' not in df.columns
        result = tx.prepare_training_data(df)  # must not raise
        assert result is not None


# ===========================================================================
# 5. prepare_training_data – legacy supp_vec path
# ===========================================================================

class TestPrepareTrainingDataLegacy:

    def _make_df_supp(self, n_callers=7):
        """Build a minimal DataFrame with a supp_vec column."""
        vec = '1' * n_callers
        return pd.DataFrame({
            'supp_vec': [vec, vec],
            'concordance': [n_callers, n_callers],
        })

    def test_legacy_path_seven_callers(self, tx):
        df = self._make_df_supp(7)
        result = tx.prepare_training_data(df, num_callers=7)
        assert 'supp_vec' not in result.columns
        for i in range(7):
            assert f'caller_{i}_flag' in result.columns

    def test_legacy_path_two_callers(self, tx):
        df = self._make_df_supp(2)
        result = tx.prepare_training_data(df, num_callers=2)
        assert 'caller_0_flag' in result.columns
        assert 'caller_1_flag' in result.columns

    def test_legacy_path_one_caller_raises(self, tx):
        df = self._make_df_supp(1)
        with pytest.raises(ValueError):
            tx.prepare_training_data(df, num_callers=1)

    def test_legacy_path_zero_callers_raises(self, tx):
        df = self._make_df_supp(0)
        with pytest.raises(ValueError):
            tx.prepare_training_data(df, num_callers=0)


# ===========================================================================
# 6. truth-label schema / join contract
# ===========================================================================

class TestTruthLabelContract:

    def test_truth_labels_required_columns_include_cnv_type(self, script_text):
        assert "{'sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label'}" in script_text

    def test_merge_uses_cnv_type_in_join_keys(self, script_text):
        assert "on=['sample_id', 'chrom', 'start', 'end', 'cnv_type']" in script_text

    def test_truth_labels_cli_help_mentions_cnv_type(self, script_text):
        assert 'sample_id, chrom, start, end, cnv_type, truth_label' in script_text
