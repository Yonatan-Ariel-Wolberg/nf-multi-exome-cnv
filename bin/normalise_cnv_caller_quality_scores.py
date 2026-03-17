import pysam

def standardize_cnv_qual(input_vcf, output_vcf, caller):
    """
    Standardizes CNV quality scores across 7 callers for Truvari integration.
    Valid callers: 'CANOES', 'CLAMMS', 'XHMM', 'GATK', 'DRAGEN', 'CNVKIT', 'INDELIBLE'
    """
    vcf = pysam.VariantFile(input_vcf, "r")
    
    # Add new FORMAT headers to preserve original metrics
    vcf.header.formats.add("OQ", "1", "Float", "Original Quality Score")
    vcf.header.formats.add("OAS", "1", "String", "Original Algorithm Score Metric Used")
    
    out = pysam.VariantFile(output_vcf, "w", header=vcf.header)
    universal_baseline = 100.0
    
    for record in vcf:
        sample = record.samples  # Assuming single-sample VCF
        
        # Safely extract existing QUAL
        orig_qual = record.qual
        if orig_qual is not None:
            sample["OQ"] = orig_qual  # Move to OQ field
            record.qual = None  # Clear the QUAL field to prepare for normalized score
        
        qual_norm = 0.0
        metric_used = "UNKNOWN"
        
        try:
            if caller == "CANOES":
                if "Q_SOME" in record.info:
                    q_some = float(record.info["Q_SOME"])
                    qual_norm = min(1000.0, universal_baseline * (q_some / 80.0))
                    metric_used = "Q_SOME"
        
            elif caller == "CLAMMS":
                if "Q_EXACT" in record.info and "Q_SOME" in record.info:
                    q_exact = float(record.info["Q_EXACT"])
                    q_some = float(record.info["Q_SOME"])
                    if q_exact >= 0.0:
                        qual_norm = min(1000.0, universal_baseline * (q_some / 500.0))
                    metric_used = "Q_SOME"

            elif caller == "XHMM":
                if "SQ" in sample and "EQ" in sample and "NDQ" in sample:
                    sq = float(sample["SQ"])
                    eq = float(sample["EQ"])
                    ndq = float(sample["NDQ"])
                    if eq >= 60.0 and ndq >= 60.0:
                        qual_norm = min(1000.0, universal_baseline * (sq / 60.0))
                    metric_used = "SQ"

            elif caller == "GATK":
                if "QS" in sample and "CN" in sample and "NP" in sample:
                    qs = float(sample["QS"])
                    cn = int(sample["CN"])
                    n_int = float(sample["NP"])
                    
                    t_gatk = None
                    if cn == 0:
                        t_gatk = min(1000.0, max(400.0, 10.0 * n_int))
                    elif cn == 1:
                        t_gatk = min(1000.0, max(100.0, 10.0 * n_int))
                    elif cn > 2:
                        t_gatk = min(400.0, max(50.0, 4.0 * n_int))
                    
                    if t_gatk is not None:
                        qual_norm = min(1000.0, universal_baseline * (qs / t_gatk))
                        metric_used = "QS"

            elif caller == "CNVKIT":
                if "CNQ" in sample:
                    cnq = float(sample["CNQ"])
                    qual_norm = min(1000.0, universal_baseline * (cnq / 20.0))
                    metric_used = "CNQ"

            elif caller == "DRAGEN":
                if orig_qual is not None and "PASS" in record.filter.keys():
                    dragen_qual = float(orig_qual)
                    if dragen_qual <= 10.0:
                        qual_norm = dragen_qual * 10.0
                    else:
                        qual_norm = 100.0 + (dragen_qual - 10.0) * (900.0 / 190.0)
                    metric_used = "QUAL"

            elif caller == "INDELIBLE":
                if all(k in record.info for k in ["sr_total", "avg_mapq", "mum_sr", "dad_sr"]):
                    sr_total = float(record.info["sr_total"])
                    avg_mapq = float(record.info["avg_mapq"])
                    mum_sr = float(record.info["mum_sr"])
                    dad_sr = float(record.info["dad_sr"])
                    
                    if sr_total >= 5.0 and avg_mapq >= 20.0 and mum_sr < 2.0 and dad_sr < 2.0:
                        synthetic_score = sr_total * (avg_mapq / 60.0) * 100.0
                        qual_norm = min(1000.0, universal_baseline * (synthetic_score / 16.67))
                        metric_used = "SYNTHETIC_PHRED"

        except (ValueError, TypeError) as e:
            print(f"Error processing record {record.id}: {e}")  # Logging error

        # Set the normalized quality score in the QUAL field
        record.qual = round(qual_norm, 2)
        sample["OAS"] = metric_used  # Store metric used without overwriting the sample object

        out.write(record)
        
    vcf.close()
    out.close()
