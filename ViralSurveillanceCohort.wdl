version 1.0

import "ViralSurveillanceSingleSample.wdl"

# -------------------------------
# STRUCTS
# -------------------------------
struct Sample {
    String sample_id
    File fastq_R1
    File fastq_R2
}
# -------------------------------
# COHORT WORKFLOW (scatter-gather)
# -------------------------------
workflow ViralSurveillanceCohort {
    input {
        Array[Sample] samples
        
        File host_fasta
        File host_fasta_amb
        File host_fasta_ann
        File host_fasta_bwt
        File host_fasta_pac
        File host_fasta_sa
        
        File viral_db
        File viral_db_amb
        File viral_db_ann
        File viral_db_bwt
        File viral_db_pac
        File viral_db_sa
        
        File kraken2_db
        
        Int min_depth = 10
        Float min_var_af = 0.03
    }
    
    scatter (s in samples) {
        call ViralSurveillanceSingleSample.ViralSurveillanceSingleSample {
            input:
            sample_id=s.sample_id,
            fastq_R1=s.fastq_R1,
            fastq_R2=s.fastq_R2,
            host_fasta      = host_fasta,
            host_fasta_amb  = host_fasta_amb,
            host_fasta_ann  = host_fasta_ann,
            host_fasta_bwt  = host_fasta_bwt,
            host_fasta_pac  = host_fasta_pac,
            host_fasta_sa   = host_fasta_sa,
            viral_db=viral_db,
            viral_db_amb=viral_db_amb,
            viral_db_ann=viral_db_ann,
            viral_db_bwt=viral_db_bwt,
            viral_db_pac=viral_db_pac,
            viral_db_sa=viral_db_sa,
            kraken2_db=kraken2_db,
            min_depth=min_depth,
            min_var_af=min_var_af
        }
    }
    
    # Gather files for global MultiQC
    Array[File] all_fastqc_raw     = select_all(flatten(ViralSurveillanceSingleSample.fastqc_raw_reports))
    Array[File] all_fastqc_trimmed = select_all(flatten(ViralSurveillanceSingleSample.fastqc_trimmed_reports))
    Array[File] all_fastp_json     = ViralSurveillanceSingleSample.fastp_json
    Array[File] all_kraken_report  = ViralSurveillanceSingleSample.kraken_report
    Array[File] all_samtools_stats = ViralSurveillanceSingleSample.samtools_stats
    Array[File] all_samtools_cov   = ViralSurveillanceSingleSample.samtools_coverage

    call ViralSurveillanceSingleSample.MultiQC as GlobalMultiQC {
        input:
        qc_inputs = flatten([
        all_fastqc_raw,
        all_fastqc_trimmed,
        all_fastp_json,
        all_kraken_report,
        all_samtools_stats,
        all_samtools_cov
        ]),
        output_prefix = "viral_surveillance_global"
    }
    
    output {
        # Per-sample outputs
        Array[File] per_sample_multiqc  = ViralSurveillanceSingleSample.multiqc_report
        Array[File] kraken_reports      = ViralSurveillanceSingleSample.kraken_report
        Array[File] vcfs                = ViralSurveillanceSingleSample.vcf
        Array[File] consensus_fastas    = ViralSurveillanceSingleSample.consensus_fasta
        #       Array[File] phylogeny_trees     = ViralSurveillanceSingleSample.phylogeny_tree
        
        # Global QC
        File global_multiqc_report      = GlobalMultiQC.report_html
        File global_multiqc_data        = GlobalMultiQC.report_data
    }
}
