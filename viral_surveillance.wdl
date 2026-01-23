version 1.0

workflow ViralSurveillancePanel {

    input {
        Array[File] fastq_R1
        Array[File] fastq_R2
        File human_reference
        File viral_refs_fasta
        File viral_refs_index  # bwa/bowtie2 indices for viral refs
        File kraken2_db
        Int min_depth = 10
        Float min_var_af = 0.03
      }
      
    #
    # Raw read QC
    #
    call FastQC as FastQC_Raw {
        input:
            fastqs = [fastq_R1, fastq_R2],
            prefix = "raw"
    }

    #
    # Trimming
    #
    call FastpTrim {
        input:
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2
    }

    #
    # Post-trim QC
    #
    call FastQC as FastQC_Trimmed {
        input:
            fastqs = [FastpTrim.trimmed_R1, FastpTrim.trimmed_R2],
            prefix = "trimmed"
    }

    call HostDepletion {
        input:
            trimmed_R1=FastpTrim.trimmed_R1,
            trimmed_R2=FastpTrim.trimmed_R2,
            human_reference=human_reference
    }

    call Kraken2Detect {
        input:
            nohost_R1=HostDepletion.nohost_R1,
            nohost_R2=HostDepletion.nohost_R2,
            kraken2_db=kraken2_db
    }

    call AlignVirusRefs {
        input:
            nohost_R1=HostDepletion.nohost_R1,
            nohost_R2=HostDepletion.nohost_R2,
            viral_refs_index=viral_refs_index
    }

    call VariantCalling {
        input:
            virus_bam=AlignVirusRefs.virus_bam,
            viral_refs_fasta=viral_refs_fasta,
            min_depth=min_depth,
            min_var_af=min_var_af
    }

    call ConsensusGenome {
        input:
            virus_bam=AlignVirusRefs.virus_bam,
            viral_refs_fasta=viral_refs_fasta,
            min_depth=min_depth
    }

    call Phylogeny {
        input:
            consensus_fasta=ConsensusGenome.consensus_fasta
    }

    #
    # MultiQC aggregation
    #
    call MultiQC {
        input:
            qc_inputs = [FastQC_Raw.html_reports, FastQC_Trimmed.html_reports, [FastpTrim.html]],
            output_prefix = "viral_surveillance_multiqc"
    }

    output {
      File multiqc_report = MultiQC.report_html
      Array[File] fastqc_raw = FastQC_Raw.html_reports
      Array[File] fastqc_trimmed = FastQC_Trimmed.html_reports
      File kraken_report = Kraken2Detect.report
      File variants_vcf = VariantCalling.variants_vcf
      File consensus_fasta = ConsensusGenome.consensus_fasta
      File phylogeny_tree = Phylogeny.tree
      File virus_bam = AlignVirusRefs.virus_bam
    }
}

task FastQC {
    input {
        Array[File] fastqs
        String prefix = "fastqc"
    }

    command {
        mkdir -p fastqc_out

        fastqc \
          --threads 4 \
          --outdir fastqc_out \
          ~{sep=" " fastqs}
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
  
task MultiQC {
    input {
        Array[File] qc_inputs
        String output_prefix = "multiqc"
    }

    command {
        mkdir -p multiqc_out

        multiqc \
          --force \
          --outdir multiqc_out \
          --filename ~{output_prefix}.html \
          ~{sep=" " qc_inputs}
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

task FastpTrim {
    input {
        Array[File] fastq_R1
        Array[File] fastq_R2
    }

    command {
        fastp \
            -i ~{sep=" -I " fastq_R1} \
            -I ~{sep=" -i " fastq_R2} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 50 \
            --json fastp.json \
            --html fastp.html \
            -o trimmed_R1.fastq.gz \
            -O trimmed_R2.fastq.gz
    }

    output {
        File trimmed_R1 = "trimmed_R1.fastq.gz"
        File trimmed_R2 = "trimmed_R2.fastq.gz"
        File html = "fastp.html"
    }

    runtime {
        docker: "quay.io/biocontainers/fastp:0.23.2--h9ee0642_0"
        memory: "8G"
        cpu: 2
    }
}

task HostDepletion {
    input {
        File trimmed_R1
        File trimmed_R2
        File human_reference
    }

    command {
        bowtie2 --very-sensitive \
            -x ~{human_reference} \
            -1 ~{trimmed_R1} \
            -2 ~{trimmed_R2} \
            --un-conc-gz nohost.fastq.gz \
            -S host.sam

        mv nohost.1.fastq.gz nohost_R1.fastq.gz
        mv nohost.2.fastq.gz nohost_R2.fastq.gz
    }

    output {
        File nohost_R1 = "nohost_R1.fastq.gz"
        File nohost_R2 = "nohost_R2.fastq.gz"
    }

    runtime {
        docker: "quay.io/biocontainers/bowtie2:2.5.1--py39h2e2f0a9_0"
        memory: "16G"
        cpu: 4
    }
}

task Kraken2Detect {
    input {
        File nohost_R1
        File nohost_R2
        File kraken2_db
    }

    command {
        kraken2 \
            --db ~{kraken2_db} \
            --paired ~{nohost_R1} ~{nohost_R2} \
            --report report.kraken \
            --output out.kraken
    }

    output {
        File report = "report.kraken"
        File koutput = "out.kraken"
    }

    runtime {
        docker: "quay.io/biocontainers/kraken2:2.1.2--pl526hcbee51f_0"
        memory: "16G"
        cpu: 4
    }
}

task AlignVirusRefs {
    input {
        File nohost_R1
        File nohost_R2
        File viral_refs_index
    }

    command {
        bwa mem ~{viral_refs_index} ~{nohost_R1} ~{nohost_R2} \
          | samtools sort -o virus.bam

        samtools index virus.bam
    }

    output {
        File virus_bam = "virus.bam"
        File virus_bai = "virus.bam.bai"
    }

    runtime {
        docker: "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
        memory: "16G"
        cpu: 4
    }
}

task VariantCalling {
    input {
        File virus_bam
        File viral_refs_fasta
        Int min_depth
        Float min_var_af
    }

    command {
        samtools mpileup -aa -Q 20 -q 20 -f ~{viral_refs_fasta} ~{virus_bam} \
            > virus.mpileup

        lofreq call-parallel --pp-threads 4 \
            -f ~{viral_refs_fasta} \
            -o variants.vcf virus.mpileup

        # Filter variants by allele frequency and depth
        bcftools view -i "DP>=~{min_depth} && AF>=~{min_var_af}" variants.vcf > filtered.vcf
    }

    output {
        File variants_vcf = "filtered.vcf"
    }

    runtime {
        docker: "quay.io/biocontainers/lofreq:2.1.5.0--py39h9f0ad1d_1"
        memory: "8G"
        cpu: 4
    }
}

task ConsensusGenome {
    input {
        File virus_bam
        File viral_refs_fasta
        Int min_depth
    }

    command {
        samtools mpileup -aa -A -d 0 -B -f ~{viral_refs_fasta} ~{virus_bam} \
            | ivar consensus -p consensus -q 20 -m ~{min_depth}
    }

    output {
        File consensus_fasta = "consensus.fa"
    }

    runtime {
        docker: "quay.io/biocontainers/ivar:1.3.1--h5fd41b6_1"
        memory: "8G"
        cpu: 2
    }
}
task Phylogeny {
    input {
        File consensus_fasta
    }

    command {
        # Align consensus with reference set
        mafft --auto ~{consensus_fasta} > aligned.fa
        iqtree2 -s aligned.fa -m GTR+G -B 1000
    }

    output {
        File tree = "aligned.fa.treefile"
    }

    runtime {
        docker: "quay.io/biocontainers/iqtree2:2.2.2--hf8217df_0"
        memory: "16G"
        cpu: 4
    }
}
