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
        File human_reference
        File viral_refs_fasta
        File viral_refs_index
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
                human_reference=human_reference,
                viral_refs_fasta=viral_refs_fasta,
                viral_refs_index=viral_refs_index,
                kraken2_db=kraken2_db,
                min_depth=min_depth,
                min_var_af=min_var_af
        }
    }

    # Gather files for global MultiQC
    Array[File] all_fastqc_raw     = flatten(ViralSurveillanceSingleSample.fastqc_raw_reports)
    Array[File] all_fastqc_trimmed = flatten(ViralSurveillanceSingleSample.fastqc_trimmed_reports)
    Array[File] all_fastp_html     = ViralSurveillanceSingleSample.fastp_html
    Array[File] all_kraken_summary = ViralSurveillanceSingleSample.kraken_summary_tsv
    Array[File] all_samtools_stats = ViralSurveillanceSingleSample.samtools_stats
    Array[File] all_samtools_cov   = ViralSurveillanceSingleSample.samtools_coverage

    call ViralSurveillanceSingleSample.MultiQC as GlobalMultiQC {
        input:
            qc_inputs = select_all(flatten([
                all_fastqc_raw,
                all_fastqc_trimmed,
                all_fastp_html,
                all_kraken_summary,
                all_samtools_stats,
                all_samtools_cov
            ])),
            output_prefix = "viral_surveillance_global"
    }

    output {
        # Per-sample outputs
        Array[File] per_sample_multiqc  = ViralSurveillanceSingleSample.multiqc_report
        Array[File] kraken_reports      = ViralSurveillanceSingleSample.kraken_report
        Array[File] variant_vcfs        = ViralSurveillanceSingleSample.variants_vcf
        Array[File] consensus_fastas    = ViralSurveillanceSingleSample.consensus_fasta
        Array[File] phylogeny_trees     = ViralSurveillanceSingleSample.phylogeny_tree

        # Global QC
        File global_multiqc_report      = GlobalMultiQC.report_html
        File global_multiqc_data        = GlobalMultiQC.report_data
    }
}
