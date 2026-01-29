version 1.0

workflow HostDepletionWorkflow {
    input {
        File trimmed_R1
        File trimmed_R2

        File host_fasta
        File host_fasta_amb
        File host_fasta_ann
        File host_fasta_bwt
        File host_fasta_pac
        File host_fasta_sa

        String sample_id
        Int threads = 4
    }

    call BwaHostAlign {
        input:
            fastq_R1 = trimmed_R1,
            fastq_R2 = trimmed_R2,

            host_fasta      = host_fasta,
            host_fasta_amb  = host_fasta_amb,
            host_fasta_ann  = host_fasta_ann,
            host_fasta_bwt  = host_fasta_bwt,
            host_fasta_pac  = host_fasta_pac,
            host_fasta_sa   = host_fasta_sa,

            sample_id = sample_id,
            threads   = threads
    }

    # filter sam to keep only reads that did not align to host
    call SamToNoHostBam {
        input:
            sam = BwaHostAlign.sam,
            sample_id = sample_id,
            threads = threads
    }

    call NameSortBam {
        input:
            bam = SamToNoHostBam.bam,
            sample_id = sample_id,
            threads = threads
    }

    call BamToFastqNoHost {
        input:
            bam = NameSortBam.namesorted_bam,
            sample_id = sample_id
    }

    # Optional coordinate-sorted BAM for QC
    call SortBam {
        input:
            bam = SamToNoHostBam.bam,
            sample_id = sample_id,
            threads = threads
    }

    call IndexBam {
        input:
            bam = SortBam.sorted_bam
    }

    output {
        File nohost_R1 = BamToFastqNoHost.nohost_R1
        File nohost_R2 = BamToFastqNoHost.nohost_R2

        File host_coord_bam = SortBam.sorted_bam
        File host_coord_bai = IndexBam.bai
    }
}

task BwaHostAlign {
    input {
        File fastq_R1
        File fastq_R2
        File host_fasta
        File host_fasta_amb
        File host_fasta_ann
        File host_fasta_bwt
        File host_fasta_pac
        File host_fasta_sa

        String sample_id
        Int threads = 4
    }

    command {
        bwa mem -t ~{threads} \
            -o ~{sample_id}.host.sam \
            ~{host_fasta} \
            ~{fastq_R1} \
            ~{fastq_R2}
    }

    output {
        File sam = "~{sample_id}.host.sam"
    }

    runtime {
        docker: "staphb/bwa:0.7.19"
        cpu: threads
        memory: "16G"
    }
}

task SamToNoHostBam {
    input {
        File sam
        String sample_id
        Int threads = 4
    }

    command {
        # -f 12 keep only reads pairs where neither read mapped to the host
        # -F 256 exclude secondary alignments 
        samtools view \
        -@ ~{threads} \
        -o ~{sample_id}.unsorted.bam \
        -bS -f 12 -F 256 ~{sam}
    }

    output {
        File bam = "~{sample_id}.unsorted.bam"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: threads
        memory: "4G"
    }
}

task NameSortBam {
    input {
        File bam
        String sample_id
        Int threads = 4
    }

    command {
        samtools sort \
            -@ ~{threads} \
            -n \
            -o ~{sample_id}.namesorted.bam \
            ~{bam}
    }

    output {
        File namesorted_bam = "~{sample_id}.namesorted.bam"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: threads
        memory: "6G"
    }
}

task BamToFastqNoHost {
    input {
        File bam
        String sample_id
    }

    command {
        samtools fastq \
            -1 ~{sample_id}_R1_nohost.fastq.gz \
            -2 ~{sample_id}_R2_nohost.fastq.gz \
            -0 /dev/null \
            -s /dev/null \
            -n \
            ~{bam}
    }

    output {
        File nohost_R1 = "~{sample_id}_R1_nohost.fastq.gz"
        File nohost_R2 = "~{sample_id}_R2_nohost.fastq.gz"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 2
        memory: "4G"
    }
}

task SortBam {
    input {
        File bam
        String sample_id
        Int threads = 4
    }

    command {
        samtools sort \
            -@ ~{threads} \
            -o ~{sample_id}.sorted.bam \
            ~{bam}
    }

    output {
        File sorted_bam = "~{sample_id}.sorted.bam"
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
        memory: "1G"
    }
}
