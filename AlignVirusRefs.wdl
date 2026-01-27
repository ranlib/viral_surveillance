version 1.0

workflow AlignVirusWorkflow {
  input {
    File nohost_R1
    File nohost_R2
    File viral_db
    File viral_db_amb
    File viral_db_ann
    File viral_db_bwt
    File viral_db_pac
    File viral_db_sa
    String sample_id
    Int threads = 8
  }

  call BwaMemAlign {
    input:
    fastq_R1 = nohost_R1,
    fastq_R2 = nohost_R2,
    ref_fasta=viral_db,
    ref_fasta_amb=viral_db_amb,
    ref_fasta_ann=viral_db_ann,
    ref_fasta_bwt=viral_db_bwt,
    ref_fasta_pac=viral_db_pac,
    ref_fasta_sa=viral_db_sa,
    sample_id = sample_id,
    threads = threads
  }
  
  call SamToSortedBam {
    input:
    sam = BwaMemAlign.sam,
    sample_id = sample_id,
    threads = threads
  }
  
  call IndexBam {
    input:
    bam = SamToSortedBam.bam
  }
  
  output {
    File virus_bam = SamToSortedBam.bam
    File virus_bai = IndexBam.bai
  }
}

task BwaMemAlign {
    input {
      File fastq_R1
      File fastq_R2
      File ref_fasta
      File ref_fasta_amb
      File ref_fasta_ann
      File ref_fasta_bwt
      File ref_fasta_pac
      File ref_fasta_sa
      String sample_id
      Int threads = 4
    }

    command {
        set -euxo
        
        bwa mem -M -t ~{threads} \
            -o ~{sample_id}.aligned.sam \
            ~{ref_fasta} \
            ~{fastq_R1} \
            ~{fastq_R2}
    }

    output {
        File sam = "~{sample_id}.aligned.sam"
    }

    runtime {
        docker: "staphb/bwa:0.7.19"
        cpu: threads
        memory: "16G"
    }
}

task SamToSortedBam {
    input {
        File sam
        String sample_id
        Int threads = 4
    }

    command {
        samtools view -@ ~{threads} -bS ~{sam} | \
        samtools sort -@ ~{threads} -o ~{sample_id}.sorted.bam
    }

    output {
        File bam = "~{sample_id}.sorted.bam"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: threads
        memory: "8G"
    }
}

task IndexBam {
    input {
        File bam
    }

    command {
        samtools index ~{bam}
    }

    output {
        File bai = "~{bam}.bai"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 1
        memory: "2G"
    }
}
