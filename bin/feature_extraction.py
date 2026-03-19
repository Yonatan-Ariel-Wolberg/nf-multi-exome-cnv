"""
Feature extraction for CNV machine-learning models.

Extracts a rich feature matrix from a SURVIVOR- or Truvari-merged SV VCF,
or from a single-caller VCF (``merger_mode='single'``), together with
per-caller normalised VCFs (produced by
``normalise_cnv_caller_quality_scores.py``) and optional genomic annotation
files (capture BED, BAM/CRAM, reference FASTA, mappability BED).

Feature design is informed by:
  - CN-Learn (Pounraja et al. 2019, https://github.com/girirajanlab/CN_Learn),
    adapted to this pipeline's additional callers (XHMM, GATK-gCNV, CNVkit,
    DRAGEN Germline, INDELIBLE) and the QUAL_norm quality-score normalisation.
  - pd3/cnv-paper-2020 (Danecek et al. 2020,
    https://github.com/pd3/cnv-paper-2020), whose random forest used
    concordance, per-caller quality, read-depth fraction, probe counts, GC,
    mappability, and probe-level log2 ratio (L2R) statistics (l2r_mean,
    l2r_dev, l2r_flank_diff) as the most discriminating features. The nflank
    (flanking probe count) and L2R features from that paper are implemented
    here as n_probes_flank, l2r_mean, l2r_dev, and l2r_flank_diff.

Partial-caller / missing-caller behaviour
------------------------------------------
No single caller is required.  Features degrade gracefully:

  - Structural, genomic, and L2R features (chrom, size, GC, mappability,
    n_probes, n_probes_flank, rd_ratio, l2r_*) require only the merged VCF
    plus optional BAM/BED/FASTA/mappability files.  They do not depend on
    which callers were run.

  - Per-caller flags (is_{caller}) and quality scores (qual_norm_{caller}) are
    set to 0 / NaN for any caller absent from ``tool_vcfs``.  Only the callers
    whose normalised VCFs are provided contribute real values.

  - Caller-specific secondary INFO metrics (xhmm_rd, cnvkit_weight, etc.) are
    NaN when that caller is not in ``tool_vcfs``.

  - INDELIBLE split-read features (total_sr, sr_entropy, mapq_avg, dual_split)
    are NaN when ``indelible_counts`` is None / empty; INDELIBLE does not need
    to be run for the other features to be extracted.

  - In SURVIVOR mode the SUPP_VEC bit-to-caller mapping is auto-detected from
    ``##SAMPLE`` header lines; callers not provided in ``tool_vcfs`` still
    receive an ``is_`` column (derived from SUPP_VEC) but their
    ``qual_norm_`` column is NaN.

  - In single mode (``merger_mode='single'``) the input VCF comes from one
    caller only.  Every record is treated as supported (concordance = 1) and
    is_{caller} = 1 for all callers listed in ``tool_vcfs``; qual_norm values
    are looked up from the corresponding normalised VCF.  Use this mode when
    running feature extraction without a prior multi-caller merge step.

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

Concordance / overlap  (CN-Learn: NUM_OVERLAPS; pd3: concordance)
    concordance -- number of callers supporting this event

Per-caller flags  (CN-Learn: caller_list binary indicators; pd3: per-caller cols)
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

Genomic annotations  (CN-Learn: RD_PROP, GC, MAP, NUM_TARGETS;
                      pd3: read-depth-fraction, GC, mappable, number-of-probes)
    n_probes        -- number of capture-target probes overlapping the CNV
                       (NaN when bed_file is absent)
    n_probes_flank  -- number of capture-target probes in the flanking region
                       (outside CNV, within half-CNV-length on each side);
                       mirrors pd3's nflank feature.
                       (NaN when bed_file is absent)
    rd_ratio        -- mean read depth in CNV / mean read depth in flanking
                       regions (NaN when bam_file is absent)
    gc_content      -- GC fraction of the CNV interval
                       (NaN when reference_fasta is absent)
    mappability     -- weighted-mean mappability score over the CNV interval
                       (NaN when mappability_file is absent)

Probe-level log2 ratio statistics  (pd3: l2r_mean, l2r_dev, l2r_flank_diff)
    Computed per-probe from BAM depths normalised by flanking-probe baseline.
    All three are NaN when bed_file or bam_file is absent.
    l2r_mean       -- mean log2(probe_depth / flank_baseline) across CNV probes;
                      negative for DEL, positive for DUP.
    l2r_dev        -- variance of per-probe log2 ratios within the CNV;
                      low for clean CNV signals, high for noisy/false calls.
    l2r_flank_diff -- l2r_mean minus the mean flanking-probe log2 ratio;
                      captures how much the CNV region deviates from its
                      immediate neighbourhood (robust to whole-genome shifts).

Caller-specific secondary metrics  (NaN when not present in the VCF)
    xhmm_rd        -- XHMM RD Z-score (INFO field)
    cnvkit_weight  -- CNVkit bin weight (INFO field)
    cnvkit_log2    -- CNVkit log2 ratio (INFO field)
    dragen_sm      -- DRAGEN sample median (INFO field)
    dragen_sd      -- DRAGEN sample SD (INFO field)

INDELIBLE split-read metrics  (NaN when not matched in indelible_counts)
    total_sr, sr_entropy, mapq_avg, dual_split
"""

import argparse
import re

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


def _count_probes_flank(chrom, start, end, bed_intervals):
    """Count BED probes in the flanking region around [start, end).

    The flank window extends half the CNV length on each side.  Probes that
    overlap the CNV itself are excluded so only the immediately neighbouring
    targets are counted.  This mirrors the ``nflank`` feature from
    pd3/cnv-paper-2020 (Danecek et al. 2020).

    Parameters
    ----------
    chrom : str
    start, end : int
        Half-open CNV coordinates.
    bed_intervals : dict[str, list[tuple[int, int]]]
        Pre-loaded BED intervals as returned by ``_load_bed``.

    Returns
    -------
    int
        Number of BED probes in the flanking region.
    """
    cnv_len = end - start
    flank = max(1, cnv_len // 2)
    left_start = max(0, start - flank)
    right_end = end + flank
    count = 0
    for iv_start, iv_end in bed_intervals.get(chrom, []):
        in_flank = iv_start < right_end and iv_end > left_start
        in_cnv = iv_start < end and iv_end > start
        if in_flank and not in_cnv:
            count += 1
    return count


# ── SURVIVOR caller-order helper ──────────────────────────────────────────────

# Recognises caller names embedded in VCF file paths from SURVIVOR ##SAMPLE lines.
_CALLER_FILENAME_PATTERNS = [
    ('canoes',    re.compile(r'CANOES',    re.IGNORECASE)),
    ('clamms',    re.compile(r'CLAMMS',    re.IGNORECASE)),
    ('xhmm',      re.compile(r'XHMM',      re.IGNORECASE)),
    ('gatk_gcnv', re.compile(r'GCN[V]|GCNV|GATK', re.IGNORECASE)),
    ('cnvkit',    re.compile(r'CNVKit|CNVKIT', re.IGNORECASE)),
    ('dragen',    re.compile(r'DRAGEN',    re.IGNORECASE)),
    ('indelible', re.compile(r'INDELIBLE', re.IGNORECASE)),
]


def _caller_order_from_survivor_header(vcf_header):
    """Infer SUPP_VEC bit-to-caller mapping from SURVIVOR ##SAMPLE header lines.

    SURVIVOR writes one ``##SAMPLE=<ID=N,File=path>`` meta-information line
    per input VCF, where N is the 0-based bit position in SUPP_VEC.  This
    helper parses those lines and matches file paths against known caller
    filename patterns to produce an ordered list of caller names.

    Parameters
    ----------
    vcf_header : pysam.VariantHeader
        Header object from the SURVIVOR-merged VCF.

    Returns
    -------
    list[str]
        Caller names in SUPP_VEC bit order.  Positions whose filename does not
        match any known caller receive the name ``'tool_N'``.  Returns an
        empty list when no ``##SAMPLE`` lines are found (non-SURVIVOR VCF).
    """
    # pysam exposes raw header records via vcf_header.records
    sample_entries = {}  # bit_index -> file_path
    for rec in vcf_header.records:
        if rec.type == 'GENERIC' and str(rec).startswith('##SAMPLE='):
            raw = str(rec).strip()
            id_match   = re.search(r'ID=(\d+)',     raw)
            file_match = re.search(r'File=([^,>]+)', raw)
            if id_match and file_match:
                sample_entries[int(id_match.group(1))] = file_match.group(1)

    if not sample_entries:
        return []

    order = []
    for idx in sorted(sample_entries):
        file_path = sample_entries[idx]
        name = f'tool_{idx}'
        for caller, pattern in _CALLER_FILENAME_PATTERNS:
            if pattern.search(file_path):
                name = caller
                break
        order.append(name)
    return order


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


def _l2r_stats(bam, chrom, start, end, bed_intervals, n_flank=20):
    """Probe-level log2 ratio statistics for a CNV interval.

    Inspired by the l2r_mean / l2r_dev / l2r_flank_diff features from
    pd3/cnv-paper-2020 (Danecek et al. 2020).

    For each capture probe (BED interval) overlapping the CNV, the per-probe
    read depth is normalised against the mean depth of the *n_flank* nearest
    probes on each side (the flanking baseline).  Log2 ratios are then
    computed per probe and summarised as mean, variance, and the difference
    between the CNV mean and the flanking mean (which is ~0 by definition of
    the baseline, but varies when left and right flanks are unequal).

    Parameters
    ----------
    bam : pysam.AlignmentFile
    chrom : str
    start, end : int
        Half-open CNV coordinates.
    bed_intervals : dict[str, list[tuple[int, int]]]
        Pre-loaded capture BED as returned by ``_load_bed``.
    n_flank : int
        Number of probes to use on each side of the CNV as the baseline.

    Returns
    -------
    tuple[float, float, float]
        (l2r_mean, l2r_dev, l2r_flank_diff).  All three are NaN when
        there are no CNV probes or no flanking probes.
    """
    all_probes = sorted(bed_intervals.get(chrom, []))
    cnv_probes = [(s, e) for s, e in all_probes if s < end and e > start]
    left_probes = [(s, e) for s, e in all_probes if e <= start][-n_flank:]
    right_probes = [(s, e) for s, e in all_probes if s >= end][:n_flank]

    if not cnv_probes or not (left_probes or right_probes):
        return np.nan, np.nan, np.nan

    left_depths = np.array([_mean_depth(bam, chrom, s, e) for s, e in left_probes])
    right_depths = np.array([_mean_depth(bam, chrom, s, e) for s, e in right_probes])
    all_flank_depths = np.concatenate([left_depths, right_depths])
    baseline = float(np.mean(all_flank_depths))
    if baseline == 0.0:
        return np.nan, np.nan, np.nan

    cnv_depths = np.array([_mean_depth(bam, chrom, s, e) for s, e in cnv_probes])
    # Replace zero depths with NaN before log to avoid -inf
    with np.errstate(divide='ignore', invalid='ignore'):
        cnv_l2r = np.where(cnv_depths > 0, np.log2(cnv_depths / baseline), np.nan)
    cnv_l2r = cnv_l2r[np.isfinite(cnv_l2r)]
    if cnv_l2r.size == 0:
        return np.nan, np.nan, np.nan

    l2r_mean = float(np.mean(cnv_l2r))
    l2r_dev = float(np.var(cnv_l2r))

    # Flanking mean l2r (relative to same baseline)
    with np.errstate(divide='ignore', invalid='ignore'):
        left_l2r = np.log2(left_depths / baseline) if left_depths.size > 0 else np.array([])
        right_l2r = np.log2(right_depths / baseline) if right_depths.size > 0 else np.array([])
    flank_l2r = np.concatenate([
        left_l2r[np.isfinite(left_l2r)] if left_l2r.size > 0 else np.array([]),
        right_l2r[np.isfinite(right_l2r)] if right_l2r.size > 0 else np.array([]),
    ])
    mean_flank_l2r = float(np.mean(flank_l2r)) if flank_l2r.size > 0 else 0.0
    l2r_flank_diff = l2r_mean - mean_flank_l2r

    return l2r_mean, l2r_dev, l2r_flank_diff


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
    indelible_counts=None,
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
        Mapping of caller name -> path to *normalised* per-caller VCF
        (produced by ``normalise_cnv_caller_quality_scores.py``).
        Only the callers that were actually run need to be included.
        In SURVIVOR mode the caller order for SUPP_VEC bit mapping is
        auto-detected from the ``##SAMPLE`` header lines; any callers not
        provided still receive an ``is_`` column (from SUPP_VEC) but their
        ``qual_norm_`` column will be NaN.
        Pass an empty dict ``{}`` if no per-caller normalised VCFs are
        available (structural and genomic features are still extracted).
    indelible_counts : pd.DataFrame or None, optional
        INDELIBLE small-variant count table with columns:
        Start, Total_SR, Entropy, MAPQ_Avg, Dual_Split.
        Pass ``None`` (default) when INDELIBLE was not run; the four
        split-read feature columns will be NaN for all records.
    merger_mode : {'survivor', 'truvari', 'single'}
        How the merged VCF was produced.  Use ``'single'`` when running
        feature extraction directly on a single-caller VCF without a prior
        multi-caller merge step.
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
    # ── normalise optional inputs ─────────────────────────────────────────
    if indelible_counts is None:
        indelible_counts = pd.DataFrame(
            columns=['Start', 'Total_SR', 'Entropy', 'MAPQ_Avg', 'Dual_Split']
        )

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

    # In SURVIVOR mode, try to derive the authoritative SUPP_VEC bit order
    # from the ##SAMPLE header lines so that callers not in tool_vcfs still
    # receive an is_{caller} column (NaN qual_norm is expected for those).
    if merger_mode == 'survivor':
        header_order = _caller_order_from_survivor_header(vcf_in.header)
        caller_order = header_order if header_order else list(tool_vcfs.keys())
    else:
        # truvari and single modes: caller order is simply the tool_vcfs keys.
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

        # ── flanking probe count (pd3/cnv-paper-2020: nflank) ─────────────
        v_data['n_probes_flank'] = (
            _count_probes_flank(record.chrom, record.pos, record.stop, bed_intervals)
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

        # ── probe-level L2R statistics (pd3/cnv-paper-2020) ───────────────
        if bam is not None and bed_intervals:
            l2r_mean, l2r_dev, l2r_flank_diff = _l2r_stats(
                bam, record.chrom, record.pos, record.stop, bed_intervals
            )
        else:
            l2r_mean, l2r_dev, l2r_flank_diff = np.nan, np.nan, np.nan
        v_data['l2r_mean'] = l2r_mean
        v_data['l2r_dev'] = l2r_dev
        v_data['l2r_flank_diff'] = l2r_flank_diff

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

            # Concordance for Truvari: sum of is_{caller} flags
            v_data['concordance'] = sum(
                v_data.get(f'is_{cn}', 0) for cn in caller_order
            )

        elif merger_mode == 'single':
            # Single-caller mode: the input VCF is from one caller only.
            # Every record is treated as supported by all callers in tool_vcfs,
            # so is_{caller} = 1 and concordance = 1 for each record.
            for caller_name in caller_order:
                v_data[f'is_{caller_name}'] = 1

                if caller_name in tools:
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

            # In single mode all records are supported by the single caller.
            v_data['concordance'] = 1

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


# ── Command-line interface ────────────────────────────────────────────────────
# Allows the script to be called directly from Nextflow or the shell:
#
#   python feature_extraction.py \
#     --merged_vcf  sample_survivor_union.vcf \
#     --output      sample_features.tsv \
#     [--tool_vcfs  canoes=sample_CANOES.norm.vcf.gz,clamms=sample_CLAMMS.norm.vcf.gz] \
#     [--merger_mode survivor] \
#     [--indelible_counts indelible_counts.tsv] \
#     [--bed_file    capture.bed] \
#     [--bam_file    sample.bam] \
#     [--reference_fasta GRCh38.fa] \
#     [--mappability_file mappability.bed] \
#     [--rd_flank 500] \
#     [--sample_id SAMPLE_001]

def _build_cli_parser():
    p = argparse.ArgumentParser(
        description='Extract ML feature matrix from a SURVIVOR/Truvari merged VCF.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--merged_vcf',        required=True,
                   help='SURVIVOR- or Truvari-merged VCF (bgzipped or plain).')
    p.add_argument('--output',            required=True,
                   help='Output TSV path for the feature matrix.')
    p.add_argument('--tool_vcfs',         default='',
                   help='Comma-separated caller=vcf_path pairs, e.g. '
                        'canoes=sample_CANOES.norm.vcf.gz,'
                        'clamms=sample_CLAMMS.norm.vcf.gz')
    p.add_argument('--merger_mode',       default='survivor',
                   choices=['survivor', 'truvari', 'single'],
                   help='VCF merge strategy that produced merged_vcf.  '
                        'Use "single" when running feature extraction on a '
                        'single-caller VCF without a prior merge step.')
    p.add_argument('--indelible_counts',  default=None,
                   help='INDELIBLE count TSV (Start,Total_SR,Entropy,MAPQ_Avg,'
                        'Dual_Split).  Omit if INDELIBLE was not run.')
    p.add_argument('--bed_file',          default=None,
                   help='Capture-target BED for probe counting.')
    p.add_argument('--bam_file',          default=None,
                   help='BAM or CRAM for read-depth ratio and L2R statistics.')
    p.add_argument('--reference_fasta',   default=None,
                   help='Reference FASTA (required for GC content and CRAMs).')
    p.add_argument('--mappability_file',  default=None,
                   help='4-column mappability BED (chrom,start,end,score).')
    p.add_argument('--rd_flank',          default=500, type=int,
                   help='Flanking bp for read-depth ratio calculation.')
    p.add_argument('--sample_id',         default=None,
                   help='Sample identifier (written to output column sample_id).')
    return p


if __name__ == '__main__':
    args = _build_cli_parser().parse_args()

    # Parse tool_vcfs string: "caller1=path1,caller2=path2"
    tool_vcfs = {}
    if args.tool_vcfs:
        for pair in args.tool_vcfs.split(','):
            if '=' in pair:
                name, path = pair.split('=', 1)
                tool_vcfs[name.strip()] = path.strip()

    indelible_counts = None
    if args.indelible_counts:
        indelible_counts = pd.read_csv(args.indelible_counts, sep='\t')

    df = extract_normalized_features(
        merged_vcf=args.merged_vcf,
        tool_vcfs=tool_vcfs,
        indelible_counts=indelible_counts,
        merger_mode=args.merger_mode,
        sample_id=args.sample_id,
        bed_file=args.bed_file,
        bam_file=args.bam_file,
        reference_fasta=args.reference_fasta,
        mappability_file=args.mappability_file,
        rd_flank=args.rd_flank,
    )

    if args.sample_id:
        df.insert(0, 'sample_id', args.sample_id)

    df.to_csv(args.output, sep='\t', index=False)

