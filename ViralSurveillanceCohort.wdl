version 1.0

import "ViralSurveillanceSingleSample.wdl" as SingleSample

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

        File multiqc_config
        
        Int min_depth = 10
        Float min_var_af = 0.03
        Float host_pct_cutoff = 70.0
    }
    
    scatter (s in samples) {
        call SingleSample.ViralSurveillanceSingleSample {
            input:
            sample_id       = s.sample_id,
            fastq_R1        = s.fastq_R1,
            fastq_R2        = s.fastq_R2,
            host_fasta      = host_fasta,
            host_fasta_amb  = host_fasta_amb,
            host_fasta_ann  = host_fasta_ann,
            host_fasta_bwt  = host_fasta_bwt,
            host_fasta_pac  = host_fasta_pac,
            host_fasta_sa   = host_fasta_sa,
            viral_db        = viral_db,
            viral_db_amb    = viral_db_amb,
            viral_db_ann    = viral_db_ann,
            viral_db_bwt    = viral_db_bwt,
            viral_db_pac    = viral_db_pac,
            viral_db_sa     = viral_db_sa,
            kraken2_db      = kraken2_db,
            min_depth       = min_depth,
            min_var_af      = min_var_af,
            host_pct_cutoff = host_pct_cutoff
        }
    }
    
    # Gather files for global MultiQC
    Array[File] all_fastqc_raw     = select_all(flatten(ViralSurveillanceSingleSample.fastqc_raw_reports))
    Array[File] all_fastqc_trimmed = select_all(flatten(ViralSurveillanceSingleSample.fastqc_trimmed_reports))
    Array[File] all_fastp_json     = ViralSurveillanceSingleSample.fastp_json
    Array[File] all_kraken_report  = ViralSurveillanceSingleSample.kraken_report
    Array[File] all_samtools_stats = ViralSurveillanceSingleSample.samtools_stats
    Array[File] all_samtools_cov   = ViralSurveillanceSingleSample.samtools_coverage
    Array[File] all_bcftools_stats = ViralSurveillanceSingleSample.bcftools_stats
    Array[File] host_contamination_tsvs = ViralSurveillanceSingleSample.host_contamination

    # call ViralSurveillanceSingleSample.MultiQC as GlobalMultiQC {
    #     input:
    #     qc_inputs = flatten([
    #     all_fastqc_raw,
    #     all_fastqc_trimmed,
    #     all_fastp_json,
    #     all_kraken_report,
    #     all_samtools_stats,
    #     all_samtools_cov
    #     ]),
    #     output_prefix = "viral_surveillance_global"
    # }

    call GlobalMultiQC {
        input:
        fastqc_reports         = all_fastqc_raw,
        fastqc_trimmed_reports = all_fastqc_trimmed,
        fastp_reports          = all_fastp_json,
        samtools_stats         = all_samtools_stats,
        samtools_coverage      = all_samtools_cov,
        kraken2_reports        = all_kraken_report,
        host_contamination_tsvs = host_contamination_tsvs,
        bcftools_stats         = all_bcftools_stats,
        multiqc_config         = multiqc_config
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
        #File global_multiqc_data        = GlobalMultiQC.report_data
    }
}

task GlobalMultiQC {
    input {
        Array[File] fastqc_reports
        Array[File] fastqc_trimmed_reports
        Array[File] fastp_reports
        Array[File] samtools_stats
        Array[File] samtools_coverage
        Array[File] kraken2_reports
        Array[File] host_contamination_tsvs
        Array[File] bcftools_stats
        File multiqc_config
    }

    command <<<
        set -euxo pipefail

        mkdir -p multiqc_input/{fastqc,fastp,fastqc_trimmed,samtools_stats,samtools_coverage,kraken2,host_contamination,bcftools_stats}

        # FastQC
        for f in ~{sep=' ' fastqc_reports}; do
            ln -s "$f" multiqc_input/fastqc/
        done

        # Fastp
        for f in ~{sep=' ' fastp_reports}; do
            ln -s "$f" multiqc_input/fastp/
        done

        # FastQC trimmed
        for f in ~{sep=' ' fastqc_trimmed_reports}; do
            ln -s "$f" multiqc_input/fastqc_trimmed/
        done
        # Samtools stats
        for f in ~{sep=' ' samtools_stats}; do
            ln -s "$f" multiqc_input/samtools_stats/
        done

        # Samtools coverage
        for f in ~{sep=' ' samtools_coverage}; do
            ln -s "$f" multiqc_input/samtools_coverage/
        done

        # Kraken2
        for f in ~{sep=' ' kraken2_reports}; do
            ln -s "$f" multiqc_input/kraken2/
        done

        # Host contamination
        for f in ~{sep=' ' host_contamination_tsvs}; do
            ln -s "$f" multiqc_input/host_contamination/
        done

        # vcf stats
        for f in ~{sep=' ' bcftools_stats}; do
            ln -s "$f" multiqc_input/bcftools_stats/
        done

        mkdir -p multiqc_report
        multiqc \
            --config ~{multiqc_config} \
            --outdir multiqc_report \
            multiqc_input
    >>>

    output {
        File report_html = "multiqc_report/multiqc_report.html"
        #File report_data = "multiqc_report/multiqc_data"
    }

    runtime {
        docker: "multiqc/multiqc:v1.33"
        cpu: 2
        memory: "8G"
    }
}
