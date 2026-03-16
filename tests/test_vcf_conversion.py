#!/usr/bin/env python3
"""
Tests for bin/canoes_csv_to_vcf.py, bin/clamms_bed_to_vcf.py, and
bin/indelible_tsv_to_vcf.py.

Verifies that each script:
  1. Runs without errors on representative input data.
  2. Produces a valid VCF file with the mandatory header lines.
  3. Includes TOOL=<TOOLNAME> in the INFO field of every data record.
  4. Outputs individual per-sample VCFs (one sample column each).
"""

import os
import sys
import pytest

# Allow importing the scripts directly
BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def fai_file(tmp_path):
    """A minimal FAI file covering chromosomes used in tests."""
    fai = tmp_path / "reference.fa.fai"
    fai.write_text(
        "chr1\t248956422\t112\t70\t71\n"
        "chr2\t242193529\t252513167\t70\t71\n"
        "chr3\t198295559\t498166716\t70\t71\n"
    )
    return str(fai)


@pytest.fixture()
def sample_file(tmp_path):
    """A sample-list file with two sample IDs."""
    sf = tmp_path / "samples.txt"
    sf.write_text("SAMPLE1\nSAMPLE2\n")
    return str(sf)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _read_vcf_records(vcf_path):
    """Return (header_lines, data_lines) from a VCF file."""
    header, data = [], []
    with open(vcf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
            else:
                data.append(line)
    return header, data


# ===========================================================================
# CANOES
# ===========================================================================

@pytest.fixture()
def canoes_csv(tmp_path):
    """Minimal CANOES CSV with one call per sample."""
    csv = tmp_path / "canoes.csv"
    csv.write_text(
        "CNV,INTERVAL,KB,CHR,MID_BP,TARGETS,NUM_TARG,SAMPLE,MLCN,Q_SOME\n"
        "DEL,chr1:1000-2000,1.0,1,1500,GENE1:1000-2000,2,SAMPLE1,1,90.5\n"
        "DUP,chr2:5000-8000,3.0,2,6500,GENE2:5000-8000,3,SAMPLE2,3,75.2\n"
    )
    return str(csv)


class TestCanoesCSVtoVCF:
    """Tests for bin/canoes_csv_to_vcf.py."""

    def test_individual_vcfs_created(self, canoes_csv, sample_file, fai_file, tmp_path):
        """Individual per-sample VCFs are written to the output directory."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        assert os.path.isfile(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        assert os.path.isfile(os.path.join(out, "SAMPLE2_CANOES_output.vcf"))

    def test_only_one_sample_per_vcf(self, canoes_csv, sample_file, fai_file, tmp_path):
        """Each VCF file contains exactly one sample column."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        for sample, fname in (
            ("SAMPLE1", "SAMPLE1_CANOES_output.vcf"),
            ("SAMPLE2", "SAMPLE2_CANOES_output.vcf"),
        ):
            header, _ = _read_vcf_records(os.path.join(out, fname))
            col_header = [l for l in header if l.startswith("#CHROM")][0]
            columns = col_header.split("\t")
            # FORMAT is column index 8; sample columns start at index 9
            sample_cols = columns[9:]
            assert sample_cols == [sample], (
                f"{fname}: expected sample column [{sample}], got {sample_cols}"
            )

    def test_tool_field_in_info(self, canoes_csv, sample_file, fai_file, tmp_path):
        """Every data record in individual VCFs contains TOOL=CANOES in INFO."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        for fname in ("SAMPLE1_CANOES_output.vcf", "SAMPLE2_CANOES_output.vcf"):
            _, data = _read_vcf_records(os.path.join(out, fname))
            assert data, f"{fname} contains no data records"
            for record in data:
                info = record.split("\t")[7]
                assert "TOOL=CANOES" in info, (
                    f"TOOL=CANOES missing from INFO in {fname}: {info}"
                )

    def test_vcf_header_contains_tool_meta(self, canoes_csv, sample_file, fai_file, tmp_path):
        """The VCF header declares the TOOL INFO field."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        tool_meta = [l for l in header if "##INFO=<ID=TOOL" in l]
        assert tool_meta, "##INFO=<ID=TOOL,...> not found in VCF header"

    def test_fileformat_line_present(self, canoes_csv, sample_file, fai_file, tmp_path):
        """The VCF output starts with ##fileformat=VCFv4.x."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        assert header[0].startswith("##fileformat=VCFv4"), (
            f"First header line is not ##fileformat: {header[0]}"
        )

    def test_source_line_is_canoes(self, canoes_csv, sample_file, fai_file, tmp_path):
        """The ##source header line identifies CANOES as the caller."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        source_lines = [l for l in header if l.startswith("##source=")]
        assert source_lines, "No ##source= line in VCF header"
        assert "CANOES" in source_lines[0], (
            f"##source does not mention CANOES: {source_lines[0]}"
        )

    def test_contig_lines_present(self, canoes_csv, sample_file, fai_file, tmp_path):
        """Contig lines derived from the FAI file appear in the header."""
        import canoes_csv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_canoes_data(canoes_csv, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CANOES_output.vcf"))
        contig_lines = [l for l in header if l.startswith("##contig=")]
        assert contig_lines, "No ##contig= lines found in VCF header"


# ===========================================================================
# CLAMMS
# ===========================================================================

@pytest.fixture()
def clamms_bed(tmp_path):
    """
    Minimal CLAMMS BED with 18 tab-separated fields per row (one call each
    for SAMPLE1 and SAMPLE2).
    """
    bed = tmp_path / "clamms.bed"
    rows = [
        # CHROM START  END    INTERVAL             SAMPLE   CNV  MLCN  NUM_WIN  Q_SOME  Q_EXACT  QL_EXT  L_EXT  QR_EXT  R_EXT  QL_CONT  L_CONT  QR_CONT  R_CONT
        "chr1\t1000\t2000\tchr1:1000-2000\tSAMPLE1\tDEL\t1\t5\t600\t50\t30\t800\t25\t1200\t15\t1100\t20\t1300",
        "chr2\t5000\t8000\tchr2:5000-8000\tSAMPLE2\tDUP\t3\t8\t450\t40\t20\t4000\t15\t9000\t10\t5200\t12\t7800",
    ]
    bed.write_text("\n".join(rows) + "\n")
    return str(bed)


class TestCLAMMSBEDtoVCF:
    """Tests for bin/clamms_bed_to_vcf.py."""

    def test_individual_vcfs_created(self, clamms_bed, sample_file, fai_file, tmp_path):
        """Individual per-sample VCFs are written to the output directory."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        assert os.path.isfile(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        assert os.path.isfile(os.path.join(out, "SAMPLE2_CLAMMS_output.vcf"))

    def test_only_one_sample_per_vcf(self, clamms_bed, sample_file, fai_file, tmp_path):
        """Each VCF file contains exactly one sample column."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        for sample, fname in (
            ("SAMPLE1", "SAMPLE1_CLAMMS_output.vcf"),
            ("SAMPLE2", "SAMPLE2_CLAMMS_output.vcf"),
        ):
            header, _ = _read_vcf_records(os.path.join(out, fname))
            col_header = [l for l in header if l.startswith("#CHROM")][0]
            columns = col_header.split("\t")
            sample_cols = columns[9:]
            assert sample_cols == [sample], (
                f"{fname}: expected sample column [{sample}], got {sample_cols}"
            )

    def test_tool_field_in_info(self, clamms_bed, sample_file, fai_file, tmp_path):
        """Every data record in individual VCFs contains TOOL=CLAMMS in INFO."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        for fname in ("SAMPLE1_CLAMMS_output.vcf", "SAMPLE2_CLAMMS_output.vcf"):
            _, data = _read_vcf_records(os.path.join(out, fname))
            assert data, f"{fname} contains no data records"
            for record in data:
                info = record.split("\t")[7]
                assert "TOOL=CLAMMS" in info, (
                    f"TOOL=CLAMMS missing from INFO in {fname}: {info}"
                )

    def test_vcf_header_contains_tool_meta(self, clamms_bed, sample_file, fai_file, tmp_path):
        """The VCF header declares the TOOL INFO field."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        tool_meta = [l for l in header if "##INFO=<ID=TOOL" in l]
        assert tool_meta, "##INFO=<ID=TOOL,...> not found in VCF header"

    def test_fileformat_line_present(self, clamms_bed, sample_file, fai_file, tmp_path):
        """The VCF output starts with ##fileformat=VCFv4.x."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        assert header[0].startswith("##fileformat=VCFv4"), (
            f"First header line is not ##fileformat: {header[0]}"
        )

    def test_source_line_is_clamms(self, clamms_bed, sample_file, fai_file, tmp_path):
        """The ##source header line identifies CLAMMS as the caller."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        source_lines = [l for l in header if l.startswith("##source=")]
        assert source_lines, "No ##source= line in VCF header"
        assert "CLAMMS" in source_lines[0], (
            f"##source does not mention CLAMMS: {source_lines[0]}"
        )

    def test_contig_lines_present(self, clamms_bed, sample_file, fai_file, tmp_path):
        """Contig lines derived from the FAI file appear in the header."""
        import clamms_bed_to_vcf as mod

        out = str(tmp_path / "out")
        mod.process_clamms_data(clamms_bed, sample_file, fai_file, out)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_CLAMMS_output.vcf"))
        contig_lines = [l for l in header if l.startswith("##contig=")]
        assert contig_lines, "No ##contig= lines found in VCF header"


# ===========================================================================
# INDELIBLE
# ===========================================================================

@pytest.fixture()
def indelible_tsv(tmp_path):
    """Minimal INDELIBLE annotated TSV with one PASS and one LowQuality call.

    Uses the actual 37-column INDELIBLE output format produced by
    'indelible.py annotate':
      chrom, position, coverage, insertion_context, deletion_context,
      sr_total, sr_total_long, sr_total_short, sr_long_5, sr_short_5,
      sr_long_3, sr_short_3, sr_entropy, context_entropy,
      entropy_upstream, entropy_downstream, sr_sw_similarity,
      avg_avg_sr_qual, avg_mapq, seq_longest, predicted, prob_N, prob_Y,
      ddg2p, hgnc, hgnc_constrained, exonic, transcripts, exon_numbers,
      maf, blast_hit, blast_strand, blast_identity, blast_dist,
      blast_hgnc, blast_hgnc_constrained, blast_ddg2p
    """
    tsv = tmp_path / "indelible.tsv"
    header = (
        "chrom\tposition\tcoverage\tinsertion_context\tdeletion_context\t"
        "sr_total\tsr_total_long\tsr_total_short\tsr_long_5\tsr_short_5\t"
        "sr_long_3\tsr_short_3\tsr_entropy\tcontext_entropy\t"
        "entropy_upstream\tentropy_downstream\tsr_sw_similarity\t"
        "avg_avg_sr_qual\tavg_mapq\tseq_longest\tpredicted\tprob_N\tprob_Y\t"
        "ddg2p\thgnc\thgnc_constrained\texonic\ttranscripts\texon_numbers\t"
        "maf\tblast_hit\tblast_strand\tblast_identity\tblast_dist\t"
        "blast_hgnc\tblast_hgnc_constrained\tblast_ddg2p"
    )
    # Row 1: predicted=Y → PASS; Row 2: predicted=N → LowQuality
    row1 = (
        "chr1\t1000\t50\t0\t0\t15\t10\t5\t6\t3\t4\t2\t"
        "1.5\t1.2\t1.3\t1.1\t1.0\t35.0\t29.0\t"
        "ACGTACGT\tY\t0.05\t0.95\t"
        "GENE1\tGENE1\t\tFalse\t\t\t\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
    )
    row2 = (
        "chr2\t5000\t30\t0\t0\t8\t6\t2\t3\t1\t3\t1\t"
        "1.7\t1.6\t1.7\t1.5\t1.0\t32.8\t29.0\t"
        "TTGC\tN\t0.70\t0.30\t"
        "NA\tNA\t\tFalse\t\t\t\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
    )
    tsv.write_text(header + "\n" + row1 + "\n" + row2 + "\n")
    return str(tsv)


class TestINDELIBLETSVtoVCF:
    """Tests for bin/indelible_tsv_to_vcf.py."""

    def test_vcf_file_created(self, indelible_tsv, fai_file, tmp_path):
        """A per-sample VCF file is written to the output directory."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        assert os.path.isfile(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))

    def test_only_one_sample_per_vcf(self, indelible_tsv, fai_file, tmp_path):
        """Each INDELIBLE VCF file contains exactly one sample column."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        col_header = [l for l in header if l.startswith("#CHROM")][0]
        columns = col_header.split("\t")
        sample_cols = columns[9:]
        assert sample_cols == ["SAMPLE1"], (
            f"Expected sample column [SAMPLE1], got {sample_cols}"
        )

    def test_tool_field_in_info(self, indelible_tsv, fai_file, tmp_path):
        """Every data record contains TOOL=INDELIBLE in INFO."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        assert data, "INDELIBLE VCF contains no data records"
        for record in data:
            info = record.split("\t")[7]
            assert "TOOL=INDELIBLE" in info, (
                f"TOOL=INDELIBLE missing from INFO: {info}"
            )

    def test_vcf_header_contains_tool_meta(self, indelible_tsv, fai_file, tmp_path):
        """The VCF header declares the TOOL INFO field."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        tool_meta = [l for l in header if "##INFO=<ID=TOOL" in l]
        assert tool_meta, "##INFO=<ID=TOOL,...> not found in VCF header"

    def test_fileformat_line_present(self, indelible_tsv, fai_file, tmp_path):
        """The VCF output starts with ##fileformat=VCFv4.x."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        assert header[0].startswith("##fileformat=VCFv4"), (
            f"First header line is not ##fileformat: {header[0]}"
        )

    def test_source_line_is_indelible(self, indelible_tsv, fai_file, tmp_path):
        """The ##source header line identifies INDELIBLE as the caller."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        source_lines = [l for l in header if l.startswith("##source=")]
        assert source_lines, "No ##source= line in VCF header"
        assert "INDELIBLE" in source_lines[0], (
            f"##source does not mention INDELIBLE: {source_lines[0]}"
        )

    def test_contig_lines_present(self, indelible_tsv, fai_file, tmp_path):
        """Contig lines derived from the FAI file appear in the header."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        header, _ = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        contig_lines = [l for l in header if l.startswith("##contig=")]
        assert contig_lines, "No ##contig= lines found in VCF header"

    def test_pass_filter_for_predicted_y(self, indelible_tsv, fai_file, tmp_path):
        """Rows with predicted=Y receive PASS filter status."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        # First row has predicted=Y → chr1:1000
        first = data[0].split("\t")
        assert first[6] == "PASS", f"Expected PASS filter for predicted=Y, got {first[6]}"

    def test_lowquality_filter_for_predicted_n(self, indelible_tsv, fai_file, tmp_path):
        """Rows with predicted=N receive LowQuality filter status."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file)

        _, data = _read_vcf_records(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
        # Second row has predicted=N → chr2:5000
        second = data[1].split("\t")
        assert second[6] == "LowQuality", (
            f"Expected LowQuality filter for predicted=N, got {second[6]}"
        )

    def test_works_without_fai_file(self, indelible_tsv, tmp_path):
        """Script runs successfully even when no FAI file is supplied."""
        import indelible_tsv_to_vcf as mod

        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_tsv, out, "SAMPLE1", fai_file=None)

        assert os.path.isfile(os.path.join(out, "SAMPLE1_INDELIBLE_output.vcf"))
