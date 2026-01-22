#!/usr/bin/env bash
set -euo pipefail

FASTA=viral_db/refseq/fasta/viral_refseq.fasta
OUTDIR=viral_db/refseq/bowtie2

mkdir -p ${OUTDIR}
cp ${FASTA} ${OUTDIR}/
cd ${OUTDIR}

echo "Building Bowtie2 index..."
bowtie2-build viral_refseq.fasta viral_refseq

echo "Bowtie2 index complete:"
ls -lh viral_refseq.*
