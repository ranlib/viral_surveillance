version 1.0

# -------------------------------
# TASKS
# -------------------------------

# FastQC
task FastQC {
    input {
        Array[File] fastqs
        String prefix = "fastqc"
    }

    command {
        mkdir -p fastqc_out
        fastqc --threads 4 --outdir fastqc_out ~{sep=" " fastqs}
    }

    output {
        Array[File] html_reports = glob("fastqc_out/*_fastqc.html")
        Array[File] zip_reports  = glob("fastqc_out/*_fastqc.zip")
    }

    runtime {
        docker: "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
        cpu: 4
        memory: "4G"
    }
}

# fastp trimming
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
          -o ~{sample_id}.R1.trimmed.fastq.gz \
          -O ~{sample_id}.R2.trimmed.fastq.gz \
          --html ~{sample_id}.fastp.html \
          --json ~{sample_id}.fastp.json
    }

    output {
        File trimmed_R1 = "~{sample_id}.R1.trimmed.fastq.gz"
        File trimmed_R2 = "~{sample_id}.R2.trimmed.fastq.gz"
        File html       = "~{sample_id}.fastp.html"
    }

    runtime {
        docker: "dbest/fastp:v1.0.1"
        cpu: 4
        memory: "4G"
    }
}

# Host depletion
task HostDepletionBwa {
    input {
      File trimmed_R1
      File trimmed_R2
      File human_reference
      String sample_id
      Int threads = 32
    }

    command {
        mkdir -p host_depletion_out
        bwa mem -t ${hreads} ~{human_reference} ~{trimmed_R1} ~{trimmed_R2} | \
            samtools view -b -f 12 -F 256 - > host_depletion_out/~{sample_id}.nohost.bam
        # Extract no-host FASTQs
        samtools fastq host_depletion_out/~{sample_id}.nohost.bam \
            -1 host_depletion_out/~{sample_id}.R1.nohost.fastq.gz \
            -2 host_depletion_out/~{sample_id}.R2.nohost.fastq.gz
    }

    output {
        File nohost_R1 = "host_depletion_out/~{sample_id}.R1.nohost.fastq.gz"
        File nohost_R2 = "host_depletion_out/~{sample_id}.R2.nohost.fastq.gz"
    }

    runtime {
        docker: "biocontainers/bwa:v0.7.17_cv1"
        cpu: 4
        memory: "16G"
    }
}

task HostDepletionBowtie2 {
    input {
      File trimmed_R1
      File trimmed_R2
      File human_reference
      String sample_id
      Int threads = 32
    }

    command {
      bowtie2 --very-sensitive \
      --threads ~{threads} \
      -x ~{human_reference} \
      -1 ~{trimmed_R1} \
      -2 ~{trimmed_R2} \
      --un-conc-gz ~{sample_id}.nohost.fastq.gz \
      -S host.sam
      
      mv ~{sample_id}.nohost.1.fastq.gz ~{sample_id}.nohost_R1.fastq.gz
      mv ~{sample_id}.nohost.2.fastq.gz ~{sample_id}.nohost_R2.fastq.gz
    }

    output {
        File nohost_R1 = "~{sample_id}.nohost_R1.fastq.gz"
        File nohost_R2 = "~{sample_id}.nohost_R2.fastq.gz"
    }

    runtime {
        docker: "quay.io/biocontainers/bowtie2:2.5.1--py39h2e2f0a9_0"
        memory: "16G"
        cpu: 4
    }
}

# Kraken2 detection + summary
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
        docker: "quay.io/biocontainers/kraken2:2.1.2--pl526hcbee51f_0"
        cpu: 4
        memory: "16G"
    }
}

# Viral alignment
task AlignVirusRefs {
    input {
        File nohost_R1
        File nohost_R2
        File viral_refs_index
        String sample_id
    }

    command {
        bwa mem -t 4 ~{viral_refs_index} ~{nohost_R1} ~{nohost_R2} | \
            samtools sort -o ~{sample_id}.virus.bam
        samtools index ~{sample_id}.virus.bam
    }

    output {
        File virus_bam = "~{sample_id}.virus.bam"
        File virus_bai = "~{sample_id}.virus.bam.bai"
    }

    runtime {
        docker: "biocontainers/bwa:v0.7.17_cv1"
        cpu: 4
        memory: "16G"
    }
}

# Variant calling
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
        File variants_vcf = "~{sample_id}.variants.vcf"
    }

    runtime {
        docker: "biocontainers/bcftools:v1.17-1-deb_cv1"
        cpu: 4
        memory: "8G"
    }
}

# Consensus genome
task ConsensusGenome {
    input {
        File virus_bam
        File viral_refs_fasta
        Int min_depth
        String sample_id
    }

    command {
        bcftools mpileup -f ~{viral_refs_fasta} ~{virus_bam} | \
            bcftools call -c | \
            bcftools consensus -f ~{viral_refs_fasta} -o ~{sample_id}.consensus.fasta
    }

    output {
        File consensus_fasta = "~{sample_id}.consensus.fasta"
    }

    runtime {
        docker: "biocontainers/bcftools:v1.17-1-deb_cv1"
        cpu: 4
        memory: "8G"
    }
}

# Phylogeny (basic)
task Phylogeny {
    input {
        File consensus_fasta
        String sample_id
    }

    command {
        # Example: FastTree
        FastTree -nt ~{consensus_fasta} > ~{sample_id}.tree
    }

    output {
        File tree = "~{sample_id}.tree"
    }

    runtime {
        docker: "biocontainers/fasttree:2.1.11--hdfd78af_0"
        cpu: 2
        memory: "4G"
    }
}

# Samtools stats
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
        docker: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        cpu: 1
        memory: "2G"
    }
}

# Samtools coverage
task SamtoolsCoverage {
    input {
        File bam
        File bai
        String sample_id
    }

    command {
        samtools coverage -a ~{bam} > ~{sample_id}.coverage.tsv
    }

    output {
        File coverage_tsv = "~{sample_id}.coverage.tsv"
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        cpu: 1
        memory: "4G"
    }
}

# MultiQC
task MultiQC {
    input {
        Array[File] qc_inputs
        String output_prefix = "multiqc"
    }

    command {
        mkdir -p multiqc_out
        multiqc --force --outdir multiqc_out --filename ~{output_prefix}.html ~{sep=" " qc_inputs}
    }

    output {
        File report_html = "multiqc_out/~{output_prefix}.html"
        File report_data = "multiqc_out/~{output_prefix}_data.zip"
    }

    runtime {
        docker: "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
        cpu: 2
        memory: "4G"
    }
}

# -------------------------------
# SINGLE SAMPLE WORKFLOW
# -------------------------------
workflow ViralSurveillanceSingleSample {
    input {
        String sample_id
        File fastq_R1
        File fastq_R2
        File human_reference
        File viral_refs_fasta
        File viral_refs_index
        File kraken2_db
        Int min_depth
        Float min_var_af
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

    call HostDepletionBowtie2 {
        input:
            trimmed_R1=FastpTrim.trimmed_R1,
            trimmed_R2=FastpTrim.trimmed_R2,
            human_reference=human_reference,
            sample_id=sample_id
    }

    call Kraken2Detect {
        input:
            nohost_R1=HostDepletion.nohost_R1,
            nohost_R2=HostDepletion.nohost_R2,
            kraken2_db=kraken2_db,
            sample_id=sample_id
    }

    call AlignVirusRefs {
        input:
            nohost_R1=HostDepletion.nohost_R1,
            nohost_R2=HostDepletion.nohost_R2,
            viral_refs_index=viral_refs_index,
            sample_id=sample_id
    }

    call VariantCalling {
        input:
            virus_bam=AlignVirusRefs.virus_bam,
            viral_refs_fasta=viral_refs_fasta,
            min_depth=min_depth,
            min_var_af=min_var_af,
            sample_id=sample_id
    }

    call ConsensusGenome {
        input:
            virus_bam=AlignVirusRefs.virus_bam,
            viral_refs_fasta=viral_refs_fasta,
            min_depth=min_depth,
            sample_id=sample_id
    }

    call Phylogeny {
        input:
            consensus_fasta=ConsensusGenome.consensus_fasta,
            sample_id=sample_id
    }

    call SamtoolsStats {
        input:
            bam=AlignVirusRefs.virus_bam,
            bai=AlignVirusRefs.virus_bai,
            sample_id=sample_id
    }

    call SamtoolsCoverage {
        input:
            bam=AlignVirusRefs.virus_bam,
            bai=AlignVirusRefs.virus_bai,
            sample_id=sample_id
    }

    call MultiQC as SampleMultiQC {
        input:
            qc_inputs = [
                FastQC_Raw.html_reports,
                FastQC_Trimmed.html_reports,
                FastpTrim.html,
                Kraken2Detect.summary_tsv,
                SamtoolsStats.stats,
                SamtoolsCoverage.coverage_tsv
            ],
            output_prefix = sample_id + ".multiqc"
    }

    output {
        File multiqc_report      = SampleMultiQC.report_html
        File kraken_report       = Kraken2Detect.kraken_report
        File variants_vcf        = VariantCalling.variants_vcf
        File consensus_fasta     = ConsensusGenome.consensus_fasta
        File phylogeny_tree      = Phylogeny.tree
        Array[File] fastqc_raw_reports     = FastQC_Raw.html_reports
        Array[File] fastqc_trimmed_reports = FastQC_Trimmed.html_reports
        File fastp_html                    = FastpTrim.html
        File kraken_summary_tsv            = Kraken2Detect.summary_tsv
        File samtools_stats                = SamtoolsStats.stats
        File samtools_coverage             = SamtoolsCoverage.coverage_tsv
    }
}
