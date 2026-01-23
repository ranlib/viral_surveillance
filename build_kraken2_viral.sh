#!/usr/bin/env bash
set -euo pipefail

DB_DIR=viral_db/kraken2/db
THREADS=8

mkdir -p ${DB_DIR}

echo "Downloading viral library via FTP (no rsync)..."
kraken2-build \
  --download-library viral \
  --use-ftp \
  --db ${DB_DIR}

echo "Downloading taxonomy via FTP..."
kraken2-build \
  --download-taxonomy \
  --use-ftp \
  --db ${DB_DIR}

echo "Building Kraken2 database..."
kraken2-build \
  --build \
  --threads ${THREADS} \
  --db ${DB_DIR}

echo "Kraken2 viral DB build complete"
du -sh ${DB_DIR}
