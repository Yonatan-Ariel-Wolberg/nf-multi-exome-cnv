"""
Feature extraction for CNV machine-learning models.

Extracts a rich feature matrix from a SURVIVOR- or Truvari-merged SV VCF
together with per-caller normalised VCFs (produced by
``normalise_cnv_caller_quality_scores.py``) and optional genomic annotation
files (capture BED, BAM/CRAM, reference FASTA, mappability BED).

Feature design is informed by CN-Learn (Pounraja et al. 2019,
https://github.com/girirajanlab/CN_Learn), adapted to this pipeline's
additional callers (XHMM, GATK-gCNV, CNVkit, DRAGEN Germline, INDELIBLE) and
the QUAL_norm quality-score normalisation introduced here.

Features extracted
------------------
Structural / location
    chrom, start, end
    chrom_encoded   -- chromosome as integer (1-22=autosomes, 23=X, 24=Y, 0=other)
                       Mirrors CN-Learn's factorised CHR predictor.
    cnv_size        -- CNV length in base-pairs (continuous)
    size_label      -- ordinal size-bin (1-12) matching CN-Learn's SIZE_LABEL
                       scheme; provides a non-linear size encoding that the RF
                       can exploit in addition to the raw bp length.
    cnv_type        -- 1 for DUP, 0 for DEL/other  (CN-Learn: TYPE_IND)

Concordance / overlap  (CN-Learn: NUM_OVERLAPS)
    concordance -- number of callers supporting this event

Per-caller flags  (CN-Learn: caller_list binary indicators)
    is_{caller} -- 1 if caller supports this event, else 0

Per-caller quality scores  (absent in CN-Learn; added because QUAL_norm makes
    cross-caller quality directly comparable)
    qual_norm_{caller}   -- QUAL_norm score from the normalised VCF QUAL field.
                           This is the *only* per-caller quality column: the
                           caller-native raw scores (Q_SOME, SQ, QS, CNQ, native
                           QUAL, synthetic INDELIBLE Phred) are already the direct
                           input to QUAL_norm and are therefore redundant.

Aggregate quality summaries  (complements per-caller columns)
    max_qual_norm           -- maximum QUAL_norm across all supporting callers
    mean_qual_norm_supported -- mean QUAL_norm across callers that support this
                                event (NaN when concordance == 0)

Genomic annotations  (CN-Learn: RD_PROP, GC, MAP, NUM_TARGETS)
    n_probes     -- number of capture-target probes overlapping the CNV
                    (NaN when bed_file is absent)
    rd_ratio     -- mean read depth in CNV / mean read depth in flanking regions
                    (NaN when bam_file is absent)
    gc_content   -- GC fraction of the CNV interval
                    (NaN when reference_fasta is absent)
    mappability  -- weighted-mean mappability score over the CNV interval
                    (NaN when mappability_file is absent)

Caller-specific secondary metrics  (NaN when not present in the VCF)
    xhmm_rd        -- XHMM RD Z-score (INFO field)
    cnvkit_weight  -- CNVkit bin weight (INFO field)
    cnvkit_log2    -- CNVkit log2 ratio (INFO field)
    dragen_sm      -- DRAGEN sample median (INFO field)
    dragen_sd      -- DRAGEN sample SD (INFO field)

INDELIBLE split-read metrics  (NaN when not matched in indelible_counts)
    total_sr, sr_entropy, mapq_avg, dual_split
"""

import pysam
import pandas as pd
import numpy as np
from collections import defaultdict


# ── Chromosome encoder ────────────────────────────────────────────────────────

# Mapping mirrors CN-Learn's factorised CHR predictor (girirajanlab/CN_Learn).
_CHROM_MAP = {
    **{str(i): i for i in range(1, 23)},
    **{f'chr{i}': i for i in range(1, 23)},
    'x': 23, 'chrx': 23,
    'y': 24, 'chry': 24,
}


def _encode_chrom(chrom):
    """Return an integer encoding for a chromosome name.

    Autosomes 1-22 map to 1-22; X maps to 23; Y maps to 24.
    All other values (e.g. mitochondrial, scaffolds) map to 0.
    """
    return _CHROM_MAP.get(str(chrom).lower(), 0)


# ── CNV size label ────────────────────────────────────────────────────────────

# Ordinal size bins taken directly from CN-Learn
# (girirajanlab/CN_Learn/scripts/cn_learn.py).
_SIZE_BINS = [
    (0,        1_000,   1),   # A) < 1 KB
    (1_000,    5_000,   2),   # B) 1 KB – 5 KB
    (5_000,    10_000,  3),   # C) 5 KB – 10 KB
    (10_000,   25_000,  4),   # D) 10 KB – 25 KB
    (25_000,   50_000,  5),   # E) 25 KB – 50 KB
    (50_000,   75_000,  6),   # F) 50 KB – 75 KB
    (75_000,   100_000, 7),   # G) 75 KB – 100 KB
    (100_000,  250_000, 8),   # H) 100 KB – 250 KB
    (250_000,  500_000, 9),   # I) 250 KB – 500 KB
    (500_000,  1_000_000, 10), # J) 500 KB – 1 MB
    (1_000_000, 5_000_000, 11), # K) 1 MB – 5 MB
    (5_000_000, float('inf'), 12), # L) > 5 MB
]


def _cnv_size_label(size_bp):
    """Return an ordinal size-bin label (1-12) for a CNV of *size_bp* bases.

    Labels match CN-Learn's SIZE_LABEL scheme so that features are directly
    comparable when re-using CN-Learn-trained models or combining datasets.

    Parameters
    ----------
    size_bp : int or float
        CNV length in base-pairs (must be >= 0).

    Returns
    -------
    int
        A value from 1 (< 1 KB) to 12 (> 5 MB).
    """
    for lo, hi, label in _SIZE_BINS:
        if lo <= size_bp < hi:
            return label
    return 12  # fallback for very large values


# ── BED / probe helpers ───────────────────────────────────────────────────────

def _load_bed(bed_file):
    """Return dict chrom -> list[(start, end)] loaded from a BED file."""
    intervals = defaultdict(list)
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals[chrom].append((start, end))
    return intervals


def _count_probes(chrom, start, end, bed_intervals):
    """Count BED target intervals that overlap the half-open interval [start, end)."""
    count = 0
    for iv_start, iv_end in bed_intervals.get(chrom, []):
        if iv_start < end and iv_end > start:
            count += 1
    return count


# ── Read-depth helpers ────────────────────────────────────────────────────────

def _mean_depth(bam, chrom, start, end):
    """Compute mean read depth over [start, end) using pysam pileup."""
    if start >= end:
        return 0.0
    depths = [
        col.nsegments
        for col in bam.pileup(chrom, start, end, truncate=True, min_base_quality=0)
    ]
    return float(np.mean(depths)) if depths else 0.0


def _rd_ratio(bam, chrom, start, end, flank=500):
    """Read-depth ratio: mean depth in CNV region / mean depth in flanking regions.

    Returns NaN when flanking depth is zero (avoids division by zero).
    """
    target_depth = _mean_depth(bam, chrom, start, end)
    left_start = max(0, start - flank)
    left_depth = _mean_depth(bam, chrom, left_start, start)
    right_depth = _mean_depth(bam, chrom, end, end + flank)
    flank_depth = (left_depth + right_depth) / 2.0
    if flank_depth == 0.0:
        return np.nan
    return target_depth / flank_depth


# ── GC content ────────────────────────────────────────────────────────────────

def _gc_content(fasta, chrom, start, end):
    """Fraction of G/C bases in the genomic interval [start, end)."""
    try:
        seq = fasta.fetch(chrom, start, end).upper()
    except (ValueError, KeyError):
        return np.nan
    if not seq:
        return np.nan
    gc = sum(1 for c in seq if c in 'GC')
    return gc / len(seq)


# ── Mappability helpers ───────────────────────────────────────────────────────

def _load_mappability_bed(mappability_file):
    """Load a 4-column mappability BED (chrom, start, end, score)."""
    intervals = defaultdict(list)
    with open(mappability_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start, end, score = int(parts[1]), int(parts[2]), float(parts[3])
            intervals[chrom].append((start, end, score))
    return intervals


def _mean_mappability(chrom, start, end, mappability_intervals):
    """Weighted-mean mappability score over [start, end) from a mappability BED.

    Returns NaN when no intervals overlap or the mappability dict is empty.
    """
    if not mappability_intervals:
        return np.nan
    total_bases = 0
    weighted_sum = 0.0
    for iv_start, iv_end, score in mappability_intervals.get(chrom, []):
        ov_start = max(iv_start, start)
        ov_end = min(iv_end, end)
        if ov_end > ov_start:
            bases = ov_end - ov_start
            total_bases += bases
            weighted_sum += score * bases
    if total_bases == 0:
        return np.nan
    return weighted_sum / total_bases


# ── Main extraction function ──────────────────────────────────────────────────

def extract_normalized_features(
    merged_vcf,
    tool_vcfs,
    indelible_counts,
    merger_mode='survivor',
    sample_id=None,
    bed_file=None,
    bam_file=None,
    reference_fasta=None,
    mappability_file=None,
    rd_flank=500,
):
    """Extract a feature matrix from a merged SV VCF for ML model training/scoring.

    Parameters
    ----------
    merged_vcf : str
        Path to SURVIVOR- or Truvari-merged VCF.
    tool_vcfs : dict[str, str]
        Ordered mapping of caller name -> path to *normalised* per-caller VCF
        (produced by ``normalise_cnv_caller_quality_scores.py``).
        In SURVIVOR mode the dict order must match the SUPP_VEC bit positions.
    indelible_counts : pd.DataFrame
        INDELIBLE small-variant count table with columns:
        Start, Total_SR, Entropy, MAPQ_Avg, Dual_Split.
    merger_mode : {'survivor', 'truvari'}
        How the merged VCF was produced.
    sample_id : str, optional
        Sample identifier (not currently used for feature values but reserved
        for downstream labelling).
    bed_file : str, optional
        BED file of capture target regions used to count probes per CNV.
    bam_file : str, optional
        BAM or CRAM alignment file used to compute read-depth ratio.
    reference_fasta : str, optional
        Indexed FASTA reference used for GC-content calculation and (when
        bam_file is a CRAM) as the CRAM reference.
    mappability_file : str, optional
        4-column BED file (chrom, start, end, score) with mappability scores.
    rd_flank : int
        Flanking bases on each side used for RD-ratio calculation (default 500).

    Returns
    -------
    pd.DataFrame
        One row per merged SV record; columns described in the module docstring.
    """
    # ── load optional annotation sources ─────────────────────────────────
    bed_intervals = _load_bed(bed_file) if bed_file else {}
    mappability_intervals = (
        _load_mappability_bed(mappability_file) if mappability_file else {}
    )

    bam = None
    if bam_file:
        if bam_file.endswith('.cram'):
            bam = pysam.AlignmentFile(
                bam_file, 'rc', reference_filename=reference_fasta
            )
        else:
            bam = pysam.AlignmentFile(bam_file, 'rb')

    fasta = pysam.FastaFile(reference_fasta) if reference_fasta else None

    # ── open VCFs ─────────────────────────────────────────────────────────
    vcf_in = pysam.VariantFile(merged_vcf)
    tools = {k: pysam.VariantFile(v) for k, v in tool_vcfs.items()}
    caller_order = list(tool_vcfs.keys())

    all_records = []

    for record in vcf_in:
        v_data = {
            'chrom':         record.chrom,
            'start':         record.pos,
            'end':           record.stop,
            'chrom_encoded': _encode_chrom(record.chrom),
            'cnv_size':      record.stop - record.pos,
            'size_label':    _cnv_size_label(record.stop - record.pos),
            'cnv_type':      1 if 'DUP' in str(record.info.get('SVTYPE', '')) else 0,
        }

        # ── GC content ────────────────────────────────────────────────────
        v_data['gc_content'] = (
            _gc_content(fasta, record.chrom, record.pos, record.stop)
            if fasta is not None
            else np.nan
        )

        # ── probe count ───────────────────────────────────────────────────
        v_data['n_probes'] = (
            _count_probes(record.chrom, record.pos, record.stop, bed_intervals)
            if bed_intervals
            else np.nan
        )

        # ── mappability ───────────────────────────────────────────────────
        v_data['mappability'] = _mean_mappability(
            record.chrom, record.pos, record.stop, mappability_intervals
        )

        # ── RD ratio ──────────────────────────────────────────────────────
        v_data['rd_ratio'] = (
            _rd_ratio(bam, record.chrom, record.pos, record.stop, rd_flank)
            if bam is not None
            else np.nan
        )

        # ── per-caller quality features ───────────────────────────────────
        if merger_mode == 'survivor':
            supp_vec = str(record.info.get('SUPP_VEC', '0' * len(caller_order)))
            v_data['concordance'] = sum(int(x) for x in supp_vec)

            for i, bit in enumerate(supp_vec):
                caller_name = caller_order[i] if i < len(caller_order) else f'tool_{i}'
                v_data[f'is_{caller_name}'] = int(bit)

                if bit == '1' and caller_name in tools:
                    matches = list(
                        tools[caller_name].fetch(
                            record.chrom, record.pos - 10, record.pos + 10
                        )
                    )
                    if matches:
                        orig = matches[0]
                        # QUAL_norm: normalised QUAL written by
                        # normalise_cnv_caller_quality_scores.py
                        v_data[f'qual_norm_{caller_name}'] = (
                            orig.qual if orig.qual is not None else np.nan
                        )

                        # Caller-specific secondary INFO metrics
                        if caller_name == 'xhmm' and 'RD' in orig.info:
                            v_data['xhmm_rd'] = orig.info['RD']
                        if caller_name == 'cnvkit':
                            if 'weight' in orig.info:
                                v_data['cnvkit_weight'] = orig.info['weight']
                            if 'log2' in orig.info:
                                v_data['cnvkit_log2'] = orig.info['log2']
                        if caller_name == 'dragen':
                            if 'SM' in orig.info:
                                v_data['dragen_sm'] = orig.info['SM']
                            if 'SD' in orig.info:
                                v_data['dragen_sd'] = orig.info['SD']
                    else:
                        v_data[f'qual_norm_{caller_name}'] = np.nan
                else:
                    v_data[f'qual_norm_{caller_name}'] = np.nan

        elif merger_mode == 'truvari':
            var_id = str(record.id) if record.id else ''
            for caller_name in caller_order:
                is_supp = 1 if caller_name.upper() in var_id.upper() else 0
                v_data[f'is_{caller_name}'] = is_supp

                if is_supp and caller_name in tools:
                    matches = list(
                        tools[caller_name].fetch(
                            record.chrom, record.pos - 10, record.pos + 10
                        )
                    )
                    if matches:
                        orig = matches[0]
                        v_data[f'qual_norm_{caller_name}'] = (
                            orig.qual if orig.qual is not None else np.nan
                        )

                        if caller_name == 'xhmm' and 'RD' in orig.info:
                            v_data['xhmm_rd'] = orig.info['RD']
                        if caller_name == 'cnvkit':
                            if 'weight' in orig.info:
                                v_data['cnvkit_weight'] = orig.info['weight']
                            if 'log2' in orig.info:
                                v_data['cnvkit_log2'] = orig.info['log2']
                        if caller_name == 'dragen':
                            if 'SM' in orig.info:
                                v_data['dragen_sm'] = orig.info['SM']
                            if 'SD' in orig.info:
                                v_data['dragen_sd'] = orig.info['SD']
                    else:
                        v_data[f'qual_norm_{caller_name}'] = np.nan
                else:
                    v_data[f'qual_norm_{caller_name}'] = np.nan

        # ── aggregate quality summaries across all supporting callers ─────
        # Inspired by CN-Learn's per-caller binary flags but enriched with
        # QUAL_norm: one summary captures the "best" quality signal; the
        # other captures the average quality among all agreeing callers.
        _qual_norm_vals = [
            v_data[f'qual_norm_{cn}']
            for cn in caller_order
            if v_data.get(f'is_{cn}', 0) == 1
            and not np.isnan(v_data.get(f'qual_norm_{cn}', np.nan))
        ]
        v_data['max_qual_norm'] = max(_qual_norm_vals) if _qual_norm_vals else np.nan
        v_data['mean_qual_norm_supported'] = (
            float(np.mean(_qual_norm_vals)) if _qual_norm_vals else np.nan
        )

        # ── INDELIBLE split-read counts ───────────────────────────────────
        indelible_data = indelible_counts[indelible_counts['Start'] == record.pos]
        if not indelible_data.empty:
            v_data['total_sr']   = indelible_data['Total_SR'].values[0]
            v_data['sr_entropy'] = indelible_data['Entropy'].values[0]
            v_data['mapq_avg']   = indelible_data['MAPQ_Avg'].values[0]
            v_data['dual_split'] = indelible_data['Dual_Split'].values[0]
        else:
            v_data['total_sr']   = np.nan
            v_data['sr_entropy'] = np.nan
            v_data['mapq_avg']   = np.nan
            v_data['dual_split'] = np.nan

        all_records.append(v_data)

    # ── cleanup ───────────────────────────────────────────────────────────
    if bam is not None:
        bam.close()
    if fasta is not None:
        fasta.close()

    return pd.DataFrame(all_records)


# Example usage:
# tool_vcfs = {                               # order matches SURVIVOR SUPP_VEC bits
#     'canoes':    'sample_CANOES.normalised.vcf.gz',
#     'clamms':    'sample_CLAMMS.normalised.vcf.gz',
#     'xhmm':      'sample_XHMM.normalised.vcf.gz',
#     'gatk':      'sample_GCNV.normalised.vcf.gz',
#     'cnvkit':    'sample_CNVKIT.normalised.vcf.gz',
#     'dragen':    'sample_DRAGEN.normalised.vcf.gz',
#     'indelible': 'sample_INDELIBLE.normalised.vcf.gz',
# }
# indelible_counts = pd.read_csv('indelible_counts.tsv', sep='\t')
# df = extract_normalized_features(
#     'merged.vcf',
#     tool_vcfs,
#     indelible_counts,
#     merger_mode='survivor',
#     bed_file='capture_targets.bed',
#     bam_file='sample.cram',
#     reference_fasta='GRCh38.fa',
#     mappability_file='mappability_100mer.bed',
# )
