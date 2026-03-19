#!/usr/bin/env python3
"""
Evaluate CNV caller performance (precision and sensitivity) at the probe level.

Truth set and call set BED files must contain 5 tab-separated columns in this
order:
    CHR  START  STOP  CNV_TYPE  SAMPLE_ID

The probes (capture target) BED file must contain at least 3 tab-separated
columns:
    CHR  START  STOP  [...]

Usage:
    evaluate_caller_performance.py \\
        --truth_bed   truthset.bed \\
        --callset_bed callset.bed \\
        --probes_bed  probes.bed \\
        --output      performance.txt
"""

import argparse
import sys

import pandas as pd


def load_bed_file(bed_file):
    """Load a capture-target BED file (CHR, START, STOP[, ...])."""
    return pd.read_csv(bed_file, sep='\t', header=None,
                       names=['chr', 'start', 'end'],
                       usecols=[0, 1, 2])


def load_cnv_file(cnv_file):
    """Load a 5-column CNV BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)."""
    return pd.read_csv(cnv_file, sep='\t', header=None,
                       names=['chr', 'start', 'end', 'cnv_type', 'sample'])


def find_overlaps(probes, cnvs):
    overlaps = []
    for index, probe in probes.iterrows():
        for _, cnv in cnvs.iterrows():
            if (probe['chr'] == cnv['chr'] and
                    probe['start'] < cnv['end'] and
                    probe['end'] > cnv['start']):
                overlaps.append((probe, cnv))
    return overlaps


def categorize_probes(probes, truth_cnv, callset_cnv):
    categorized = {
        'TP': [],
        'FN': [],
        'FP': [],
        'TN': [],
    }

    truth_overlaps = find_overlaps(probes, truth_cnv)
    for probe in probes.itertuples(index=False):
        is_in_truth = any(
            (probe.chr == overlap[0]['chr'] and
             probe.start < overlap[0]['end'] and
             probe.end > overlap[0]['start'])
            for overlap in truth_overlaps
        )

        is_called = any(
            (probe.chr == cnv.chr and
             probe.start < cnv.end and
             probe.end > cnv.start)
            for cnv in callset_cnv.itertuples(index=False)
        )

        if is_in_truth and is_called:
            categorized['TP'].append(probe)
        elif is_in_truth and not is_called:
            categorized['FN'].append(probe)
        elif not is_in_truth and is_called:
            categorized['FP'].append(probe)
        else:
            categorized['TN'].append(probe)

    return categorized


def compute_metrics(categorized):
    TP = len(categorized['TP'])
    FN = len(categorized['FN'])
    FP = len(categorized['FP'])
    TN = len(categorized['TN'])

    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0

    return sensitivity, precision, TN


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate CNV caller performance (precision and sensitivity).')
    parser.add_argument('--truth_bed',   required=True,
                        help='Truth set BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)')
    parser.add_argument('--callset_bed', required=True,
                        help='Call set BED file (CHR, START, STOP, CNV_TYPE, SAMPLE_ID)')
    parser.add_argument('--probes_bed',  required=True,
                        help='Capture target BED file (CHR, START, STOP[, ...])')
    parser.add_argument('--output',      default=None,
                        help='Output file for performance metrics (default: stdout)')
    args = parser.parse_args()

    probes      = load_bed_file(args.probes_bed)
    truth_cnv   = load_cnv_file(args.truth_bed)
    callset_cnv = load_cnv_file(args.callset_bed)

    categorized = categorize_probes(probes, truth_cnv, callset_cnv)
    sensitivity, precision, TN = compute_metrics(categorized)

    lines = [
        f"Sensitivity: {sensitivity}",
        f"Precision: {precision}",
        f"True Negatives: {TN}",
    ]

    if args.output:
        with open(args.output, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')
    else:
        print('\n'.join(lines))


if __name__ == '__main__':
    main()
