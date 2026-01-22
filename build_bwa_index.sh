#!/usr/bin/env bash
set -euo pipefail

FASTA=viral_db/refseq/fasta/viral_refseq.fasta
OUTDIR=viral_db/refseq/bwa

mkdir -p ${OUTDIR}
cp ${FASTA} ${OUTDIR}/

cd ${OUTDIR}

echo "Building BWA index..."
bwa index viral_refseq.fasta

echo "BWA index complete:"
ls -lh viral_refseq.fasta.*
