#!/usr/bin/env bash

set -euo pipefail

BASE_DIR=$(pwd)/viral_db/refseq
FASTA_DIR=${BASE_DIR}/fasta

mkdir -p ${FASTA_DIR}
cd ${FASTA_DIR}

echo "Downloading NCBI RefSeq viral genomes..."

wget -r -np -nd -A "*genomic.fna.gz" ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/

echo "Decompressing..."
gunzip -f *.gz

echo "Concatenating FASTA files..."
cat *.genomic.fna > viral_refseq.fasta

echo "Removing intermediate files..."
rm -f *.genomic.fna

echo "Viral RefSeq FASTA ready:"
echo "${FASTA_DIR}/viral_refseq.fasta"

