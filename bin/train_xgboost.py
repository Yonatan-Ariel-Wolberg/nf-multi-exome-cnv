import argparse
import glob
import os
import sys

import pandas as pd
import xgboost as xgb
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
import numpy as np

# All CNV callers supported by the pipeline.  Running all seven callers
# produces fully-populated is_{caller} and qual_norm_{caller} columns in the
# feature matrix.  Callers not provided will have NaN / 0 for their columns.
SUPPORTED_CALLERS = ('canoes', 'clamms', 'xhmm', 'gatk_gcnv', 'cnvkit', 'dragen', 'indelible')

# Minimum number of distinct callers required for training the ML classifier.
# With fewer than two callers, the concordance feature is always 1 and
# provides no discriminative signal; at least two callers must contribute to
# the merged input VCF.
MIN_CALLERS_FOR_TRAINING = 2


def validate_min_callers(caller_columns):
    """Raise ValueError when fewer than MIN_CALLERS_FOR_TRAINING callers are present.

    Parameters
    ----------
    caller_columns : list[str]
        Names of per-caller flag columns found in the training data.
        These are either ``is_{caller}`` names (feature-extraction output)
        or ``caller_{i}_flag`` names (legacy supp_vec split).

    Raises
    ------
    ValueError
        If ``len(caller_columns) < MIN_CALLERS_FOR_TRAINING``.
    """
    if len(caller_columns) < MIN_CALLERS_FOR_TRAINING:
        raise ValueError(
            f"At least {MIN_CALLERS_FOR_TRAINING} CNV callers are required for "
            f"training the classifier (concordance features are not meaningful "
            f"with fewer callers). "
            f"Found {len(caller_columns)} caller column(s): {list(caller_columns)}. "
            f"Supported callers: {list(SUPPORTED_CALLERS)}."
        )


def prepare_training_data(raw_df, num_callers=7):
    # Detect named is_* caller columns produced by feature_extraction.py.
    named_caller_cols = [c for c in raw_df.columns if c.startswith('is_')]

    if named_caller_cols:
        # Feature-extraction output: is_{caller} columns already present.
        validate_min_callers(named_caller_cols)
        return raw_df

    # Legacy path: raw SURVIVOR/Truvari output with a supp_vec column.
    validate_min_callers([f'caller_{i}_flag' for i in range(num_callers)])
    caller_flags = raw_df['supp_vec'].apply(lambda x: pd.Series(list(str(x))))
    caller_flags.columns = [f'caller_{i}_flag' for i in range(num_callers)]
    caller_flags = caller_flags.astype(int)

    feature_df = pd.concat([raw_df.drop(columns=['supp_vec']), caller_flags], axis=1)
    return feature_df

def train_validation_model(X, y):
    smote = SMOTE(random_state=42)
    X_res, y_res = smote.fit_resample(X, y)
    
    model = xgb.XGBClassifier(
        n_estimators=500,
        max_depth=4,
        learning_rate=0.05,
        scale_pos_weight=len(y[y==0]) / len(y[y==1]),
        use_label_encoder=False,
        missing=np.nan
    )
    
    return model, X_res, y_res

def cross_validate_model(model, X, y, n_splits=5):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = cross_val_score(model, X, y, cv=skf, scoring='f1')
    
    print(f"Cross-validated F1 scores: {scores}")
    print(f"Mean F1 score: {scores.mean()}")

def main():
    """Entry point for command-line invocation from Nextflow or the shell.

    Usage
    -----
    python train_xgboost.py \\
        --features_dir  ./out_FEATURES \\
        --truth_labels  truth_labels.tsv \\
        [--output_model  cnv_model.json] \\
        [--output_report training_report.txt]

    The ``features_dir`` is searched recursively for files matching
    ``*_features.tsv`` (the output pattern of ``feature_extraction.py``).

    The ``truth_labels`` TSV must contain at minimum the columns:
        sample_id, chrom, start, end, cnv_type, truth_label
    where ``truth_label`` is 1 for a true CNV and 0 for a false positive.
    """
    parser = argparse.ArgumentParser(
        description='Train an XGBoost classifier on CNV feature matrices.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '--features_dir', required=True,
        help='Directory (searched recursively) containing *_features.tsv files '
             'produced by feature_extraction.py.',
    )
    parser.add_argument(
        '--truth_labels', required=True,
        help='TSV file with columns: sample_id, chrom, start, end, cnv_type, truth_label.',
    )
    parser.add_argument(
        '--output_model', default='cnv_model.json',
        help='Output path for the trained XGBoost model (XGBoost JSON format).',
    )
    parser.add_argument(
        '--output_report', default='training_report.txt',
        help='Output path for the training metrics report.',
    )
    args = parser.parse_args()

    # ── Load all feature TSV files ────────────────────────────────────────────
    tsv_files = sorted(
        glob.glob(os.path.join(args.features_dir, '**', '*_features.tsv'), recursive=True)
        + glob.glob(os.path.join(args.features_dir, '*_features.tsv'))
    )
    # De-duplicate (glob may return the same file twice when features_dir == '.')
    tsv_files = list(dict.fromkeys(tsv_files))

    if not tsv_files:
        sys.exit(
            f"ERROR: No *_features.tsv files found under '{args.features_dir}'. "
            "Run the feature_extraction workflow first."
        )

    features_df = pd.concat(
        [pd.read_csv(f, sep='\t') for f in tsv_files],
        ignore_index=True,
    )

    # ── Load truth labels ─────────────────────────────────────────────────────
    labels_df = pd.read_csv(args.truth_labels, sep='\t')
    required_cols = {'sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label'}
    missing_cols = required_cols - set(labels_df.columns)
    if missing_cols:
        sys.exit(
            f"ERROR: truth_labels file is missing required columns: {missing_cols}"
        )

    # ── Merge features with truth labels ─────────────────────────────────────
    merged = features_df.merge(
        labels_df[['sample_id', 'chrom', 'start', 'end', 'cnv_type', 'truth_label']],
        on=['sample_id', 'chrom', 'start', 'end', 'cnv_type'],
        how='inner',
    )

    if merged.empty:
        sys.exit(
            "ERROR: No matching records found after joining features with truth labels. "
            "Check that sample_id / chrom / start / end / cnv_type values align."
        )

    y = merged['truth_label']
    X_raw = merged.drop(columns=['truth_label'])

    # prepare_training_data validates caller columns and handles legacy supp_vec
    X_raw = prepare_training_data(X_raw)

    # Drop non-numeric / identifier columns before passing to the model
    id_cols = [c for c in ('sample_id', 'chrom', 'cnv_type') if c in X_raw.columns]
    X = X_raw.drop(columns=id_cols).select_dtypes(include=[np.number])

    # ── Train and cross-validate ──────────────────────────────────────────────
    model, X_res, y_res = train_validation_model(X, y)
    cross_validate_model(model, X_res, y_res)
    model.fit(X_res, y_res)

    # ── Persist the trained model ─────────────────────────────────────────────
    model.save_model(args.output_model)

    # ── Write training report ─────────────────────────────────────────────────
    n_pos = int(y.sum())
    n_neg = int((y == 0).sum())
    with open(args.output_report, 'w') as fh:
        fh.write(f"Feature TSV files:         {len(tsv_files)}\n")
        fh.write(f"Training samples (calls):  {len(X)}\n")
        fh.write(f"True CNVs  (label=1):      {n_pos}\n")
        fh.write(f"False CNVs (label=0):      {n_neg}\n")
        fh.write(f"Feature columns:           {list(X.columns)}\n")
        fh.write(f"Model saved to:            {args.output_model}\n")

    print(f"Model saved to {args.output_model}")
    print(f"Report written to {args.output_report}")


if __name__ == '__main__':
    main()
