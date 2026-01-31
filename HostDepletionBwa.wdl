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
        Float host_pct_cutoff = 70.0
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
    call SamToBam {
        input:
            sam = BwaHostAlign.sam,
            sample_id = sample_id,
            threads = threads
    }

    call NameSortBam {
        input:
            bam = SamToBam.bam,
            sample_id = sample_id,
            threads = threads
    }

    call BamToFastqNoHost {
        input:
        bam = NameSortBam.namesorted_bam,
        sample_id = sample_id,
        threads = threads
    }

    # coordinate-sorted BAM for QC
    call SortBam {
        input:
            bam = SamToBam.bam,
            sample_id = sample_id,
            threads = threads
    }

    call IndexBam {
        input:
            bam = SortBam.sorted_bam
    }

    call HostContaminationMetrics {
        input:
        bam = SortBam.sorted_bam,
        sample_id = sample_id,
        host_pct_cutoff = host_pct_cutoff
    }

    output {
        File nohost_R1 = BamToFastqNoHost.nohost_R1
        File nohost_R2 = BamToFastqNoHost.nohost_R2

        File host_coord_bam = SortBam.sorted_bam
        File host_coord_bai = IndexBam.bai
        
        File host_contamination = HostContaminationMetrics.contamination_tsv
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

task SamToBam {
    input {
        File sam
        String sample_id
        Int threads = 4
    }

    command {
        samtools view \
        -@ ~{threads} \
        -o ~{sample_id}.unsorted.bam \
        -bS ~{sam}
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
        File bam # name sorted bam file
        String sample_id
        Int threads = 4
    }

    command {
        # -f 12 keep only reads pairs where neither read mapped to the host
        # -F 256 exclude secondary alignments 
        samtools fastq \
            -@ ~{threads} \
            -f 12 -F 256 \
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
        cpu: threads
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

task HostContaminationMetrics {
    input {
        File bam        # coordinate-sorted BAM from host alignment
        String sample_id
        Float host_pct_cutoff = 70.0
    }

    command <<<
        set -euxo pipefail
        # samtools view -c counts "Alignments"
        # remove secondary and supplementary alignment to count reads
        # each read is counted separatedly, for paired-end divide the number by 2
        # maybe better to use the trimmed and nohost fastq files for this
        
        # samtools flag 2304 is combination of 2 filter out flags:
        # 256 (0x100): Secondary alignments (non-primary).
        # 2048 (0x800): Supplementary alignments (chimeric or split reads).

        # -F 2316: Filters out any read that matches any of these bits.
        # This specifically excludes:
        # 1) The read itself being unmapped (4).
        # 2) The mate read being unmapped (8).
        # 3) Non-primary alignments (256).
        # 4) Chimeric/supplementary alignments (2048)

        # count number of alignments
        total_reads=$(samtools view -c -F 2304 ~{bam})
        host_mapped_reads=$(samtools view -c -F 2316 ~{bam})

        # account for paired-end reads, i.e. divide number of aligments by 2
        total_reads=$(awk -v n=$total_reads 'BEGIN{printf "%.2f", n/2}')
        host_mapped_reads=$(awk -v n=$host_mapped_reads 'BEGIN{printf "%.2f", n/2}')

        # calculate host percentage
        if [ "$total_reads" -eq 0 ]; then
            host_pct="0.00"
        else
            host_pct=$(awk -v h="$host_mapped_reads" -v t="$total_reads" 'BEGIN { printf "%.2f", (h/t)*100 }')
        fi

        # set pass/fail flag
        if awk -v n="$host_pct" -v c="~{host_pct_cutoff}" 'BEGIN {exit !(n >= c)}' ; then
            qc_flag="FAIL"
        else
            qc_flag="PASS"
        fi

        # output to tsv file
        # header
        echo -e "sample_id\ttotal_reads\thost_mapped_reads\thost_pct\tqc_flag" \
            > ~{sample_id}.host_contamination.tsv

        # data
        echo -e "~{sample_id}\t${total_reads}\t${host_mapped_reads}\t${host_pct}\t${qc_flag}" \
            >> ~{sample_id}.host_contamination.tsv
    >>>

    output {
        File contamination_tsv = "~{sample_id}.host_contamination.tsv"
    }

    runtime {
        docker: "dbest/samtools:v1.23"
        cpu: 1
        memory: "1G"
    }
}
