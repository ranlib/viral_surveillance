version 1.0

import "AlignVirusRefs.wdl"
import "HostDepletionBwa.wdl"

# -------------------------------
# SINGLE SAMPLE WORKFLOW
# -------------------------------
workflow ViralSurveillanceSingleSample {
    input {
        String sample_id

        File fastq_R1
        File fastq_R2
        
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
        
        Int min_depth
        Float min_var_af

        Float host_pct_cutoff = 70.0
    }
    
    call FastQC as FastQC_Raw {
        input: fastqs = [fastq_R1, fastq_R2]
    }
    
    call FastpTrim {
        input: fastq_R1=fastq_R1, fastq_R2=fastq_R2, sample_id=sample_id
    }
    
    call FastQC as FastQC_Trimmed {
        input: fastqs=[FastpTrim.trimmed_R1, FastpTrim.trimmed_R2]
    }
    
    call HostDepletionBwa.HostDepletionWorkflow {
        input:
        trimmed_R1      = FastpTrim.trimmed_R1,
        trimmed_R2      = FastpTrim.trimmed_R2,
        host_fasta      = host_fasta,
        host_fasta_amb  = host_fasta_amb,
        host_fasta_ann  = host_fasta_ann,
        host_fasta_bwt  = host_fasta_bwt,
        host_fasta_pac  = host_fasta_pac,
        host_fasta_sa   = host_fasta_sa,
        host_pct_cutoff = host_pct_cutoff,
        sample_id       = sample_id
    }
    
    call Kraken2Detect {
      input:
        nohost_R1=HostDepletionWorkflow.nohost_R1,
        nohost_R2=HostDepletionWorkflow.nohost_R2,
        kraken2_db=kraken2_db,
        sample_id=sample_id
    }
    
    call AlignVirusRefs.AlignVirusWorkflow {
        input:
        nohost_R1=HostDepletionWorkflow.nohost_R1,
        nohost_R2=HostDepletionWorkflow.nohost_R2,
        viral_db=viral_db,
        viral_db_amb=viral_db_amb,
        viral_db_ann=viral_db_ann,
        viral_db_bwt=viral_db_bwt,
        viral_db_pac=viral_db_pac,
        viral_db_sa=viral_db_sa,
        sample_id=sample_id
    }
    
    call VariantCalling {
        input:
        virus_bam=AlignVirusWorkflow.virus_bam,
        viral_refs_fasta=viral_db,
        min_depth=min_depth,
        min_var_af=min_var_af,
        sample_id=sample_id
    }
    
    call BgzipAndIndexVcf {
        input:
        vcf = VariantCalling.vcf
    }

    call BcftoolsStats {
        input:
        vcf = VariantCalling.vcf,
        sample_id = sample_id
    }

    call ConsensusGenome {
        input:
        vcf=BgzipAndIndexVcf.vcf_gz,
        tbi=BgzipAndIndexVcf.vcf_tbi,
        reference_fasta=viral_db,
        min_depth=min_depth,
        sample_id=sample_id
    }

    # call Phylogeny {
    #     input:
    #         consensus_fasta=ConsensusGenome.consensus_fasta,
    #         sample_id=sample_id
    # }

    call SamtoolsStats {
        input:
        bam=AlignVirusWorkflow.virus_bam,
        bai=AlignVirusWorkflow.virus_bai,
        sample_id=sample_id
    }
    
    call SamtoolsCoverage {
        input:
        bam=AlignVirusWorkflow.virus_bam,
        bai=AlignVirusWorkflow.virus_bai,
        sample_id=sample_id
    }

    Array[File] fastqc_raw_inputs = FastQC_Raw.zip_reports
    Array[File] fastqc_trimmed_inputs = FastQC_Trimmed.zip_reports
    Array[File] qc_inputs = [FastpTrim.json, Kraken2Detect.kraken_report, SamtoolsStats.stats, SamtoolsCoverage.coverage_tsv, BcftoolsStats.stats]
    call MultiQC as SampleMultiQC {
        input:
        qc_inputs = flatten([fastqc_raw_inputs, fastqc_trimmed_inputs, qc_inputs]),
        output_prefix = sample_id + ".multiqc"
    }
    
    output {
        Array[File] fastqc_raw_reports     = FastQC_Raw.zip_reports
        Array[File] fastqc_trimmed_reports = FastQC_Trimmed.zip_reports
        File fastp_json                    = FastpTrim.json
        File host_contamination            = HostDepletionWorkflow.host_contamination
        File kraken_report                 = Kraken2Detect.kraken_report
        File kraken_summary_tsv            = Kraken2Detect.summary_tsv
        File vcf                           = VariantCalling.vcf
        File bcftools_stats                = BcftoolsStats.stats
        File consensus_fasta               = ConsensusGenome.consensus_fasta
        #File phylogeny_tree               = Phylogeny.tree
        File samtools_stats                = SamtoolsStats.stats
        File samtools_coverage             = SamtoolsCoverage.coverage_tsv
        File multiqc_report                = SampleMultiQC.report_html
    }
}

# -------------------------------
# TASKS
# -------------------------------

task FastQC {
    input {
        Array[File] fastqs
        String prefix = "fastqc"
        Int threads = 32
    }
    
    command {
        mkdir -p fastqc_out
        fastqc --threads ~{threads} --outdir fastqc_out ~{sep=" " fastqs}
    }
    
    output {
        Array[File] html_reports = glob("fastqc_out/*_fastqc.html")
        Array[File] zip_reports  = glob("fastqc_out/*_fastqc.zip")
    }
    
    runtime {
        docker: "staphb/fastqc:0.12.1"
        cpu: 4
        memory: "4G"
    }
}

task FastpTrim {
    input {
        File fastq_R1
        File fastq_R2
        String sample_id
    }
    
    command {
        fastp \
        -i ~{fastq_R1} \
        -I ~{fastq_R2} \
        -o ~{sample_id}_R1_trimmed.fastq.gz \
        -O ~{sample_id}_R2_trimmed.fastq.gz \
        --html ~{sample_id}.fastp.html \
        --json ~{sample_id}.fastp.json
    }
    
    output {
        File trimmed_R1 = "~{sample_id}_R1_trimmed.fastq.gz"
        File trimmed_R2 = "~{sample_id}_R2_trimmed.fastq.gz"
        File html       = "~{sample_id}.fastp.html"
        File json       = "~{sample_id}.fastp.json"
    }
    
    runtime {
        docker: "dbest/fastp:v1.0.1"
        cpu: 4
        memory: "4G"
    }
}

task Kraken2Detect {
    input {
        File nohost_R1
        File nohost_R2
        File kraken2_db
        String sample_id
        Int threads = 1
        Int minimum_base_quality = 20
    }
    
    command <<<
        mkdir -p ${PWD}/kraken
        tar -C ${PWD}/kraken -xvf ~{kraken2_db}
        
        kraken2 \
        --db ./kraken \
        --threads ~{threads} \
        --paired ~{nohost_R1} ~{nohost_R2} \
        --unclassified-out ~{sample_id}.unclassified#.fastq \
        --classified-out ~{sample_id}.classified#.fastq \
        --minimum-base-quality ~{minimum_base_quality} \
        --use-names \
        --gzip-compressed \
        --report ~{sample_id}.kraken.report \
        --output ~{sample_id}.kraken.out
        
        # Generate MultiQC-compatible summary
        awk '$4=="S"{printf "%s\t%s\t%s\t%d\t%.5f\n","~{sample_id}", $6,$4,$2,$1/100}' ~{sample_id}.kraken.report > ~{sample_id}.kraken_summary.tsv
    >>>
    
    output {
        File kraken_report = "~{sample_id}.kraken.report"
        File kraken_output = "~{sample_id}.kraken.out"
        File summary_tsv = "~{sample_id}.kraken_summary.tsv"
    }
    
    runtime {
        docker: "staphb/kraken2:latest"
        cpu: 4
        memory: "16G"
    }
}

task VariantCalling {
    input {
        File virus_bam
        File viral_refs_fasta
        Int min_depth
        Float min_var_af
        String sample_id
    }

    command {
        bcftools mpileup -f ~{viral_refs_fasta} ~{virus_bam} | \
            bcftools call -mv -Ov -o ~{sample_id}.variants.vcf
    }

    output {
        File vcf = "~{sample_id}.variants.vcf"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 4
        memory: "8G"
    }
}

task BgzipVcf {
    input {
        File vcf
    }

    command {
        bgzip -c ~{vcf} > ~{vcf}.gz
    }

    output {
        File vcf_gz = "~{vcf}.gz"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 1
        memory: "1G"
    }
  }
  
task BgzipAndIndexVcf {
    input {
        File vcf
    }

    command {
        bgzip -c ~{vcf} > ~{vcf}.gz
        tabix -p vcf ~{vcf}.gz
    }

    output {
        File vcf_gz = "~{vcf}.gz"
        File vcf_tbi = "~{vcf}.gz.tbi"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 1
        memory: "1G"
    }
}

task BcftoolsStats {
    input {
        File vcf
        String sample_id
        Int threads = 1
    }

    command {
        set -euxo pipefail

        bcftools stats \
            --threads ~{threads} \
            ~{vcf} \
            > ~{sample_id}.bcftools.stats.txt
    }

    output {
        File stats = "~{sample_id}.bcftools.stats.txt"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 1
        memory: "1G"
    }
}

task ConsensusGenome {
    input {
        File vcf
        File tbi
        File reference_fasta
        Int min_depth
        String sample_id
    }

    command <<<
        bcftools consensus -f ~{reference_fasta} -o ~{sample_id}.consensus.fasta ~{vcf}
    >>>
    
    output {
        File consensus_fasta = "~{sample_id}.consensus.fasta"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 4
        memory: "8G"
    }
}

task Phylogeny {
    input {
        File consensus_fasta
        String sample_id
    }

    command {
        FastTree -nt ~{consensus_fasta} > ~{sample_id}.tree
    }

    output {
        File tree = "~{sample_id}.tree"
    }

    runtime {
        docker: "staphb/fasttree:2.2.0"
        cpu: 2
        memory: "4G"
    }
}

task SamtoolsStats {
    input {
        File bam
        File bai
        String sample_id
    }

    command {
        samtools stats ~{bam} > ~{sample_id}.samtools.stats
    }

    output {
        File stats = "~{sample_id}.samtools.stats"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 1
        memory: "2G"
    }
}

task SamtoolsCoverage {
    input {
        File bam
        File bai
        String sample_id
    }

    command {
        samtools coverage -o ~{sample_id}.coverage.tsv ~{bam}
    }

    output {
        File coverage_tsv = "~{sample_id}.coverage.tsv"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 1
        memory: "4G"
    }
}

task MultiQC {
    input {
        Array[File] qc_inputs
        String output_prefix = "multiqc"
    }

    command {
        mkdir -p multiqc_out
        multiqc --force --zip-data-dir --outdir multiqc_out --filename ~{output_prefix}.html ~{sep=" " qc_inputs}
    }

    output {
        File report_html = "multiqc_out/~{output_prefix}.html"
        File report_data = "multiqc_out/~{output_prefix}_data.zip"
    }

    runtime {
        docker: "multiqc/multiqc:v1.33"
        cpu: 2
        memory: "4G"
    }
}
