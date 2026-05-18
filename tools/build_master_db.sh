#!/usr/bin/env bash
# Build a Sentinel-compatible master SNP database from gnomAD v4.1.
#
# Usage:
#   build_master_db.sh exomes [out_dir]    # MAF >= 5%, biallelic SNPs, ~1.8 GB
#   build_master_db.sh genomes [out_dir]   # MAF >= 1%, biallelic SNPs, ~15 GB
#
# Streams the source VCFs directly from gnomAD's public Google Cloud bucket
# and filters with bcftools so the output stays small even though the source
# files are huge. The genomes build takes a few hours over a fast connection.
#
# Requires: bcftools on PATH, internet access to storage.googleapis.com.
set -uo pipefail

MODE="${1:-exomes}"
OUT_DIR="${2:-./snp_db}"
mkdir -p "$OUT_DIR"

case "$MODE" in
    exomes)
        URL_BASE=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites
        FILTER='INFO/AF < 0.05 | INFO/AF > 0.95'
        OUT_NAME=gnomad_v41_exomes_common
        ;;
    genomes)
        URL_BASE=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites
        FILTER='INFO/AF < 0.01 | INFO/AF > 0.99'
        OUT_NAME=gnomad_v41_genomes_maf01
        ;;
    *)
        echo "usage: $0 (exomes|genomes) [out_dir]" >&2
        exit 2
        ;;
esac

LOG="$OUT_DIR/build.log"
: > "$LOG"

MAX_PARALLEL=6
PIDS=()

filter_chrom() {
    local c="$1"
    echo "[$(date +%H:%M:%S)] start chr${c}" >> "$LOG"
    if bcftools view -e "$FILTER" -m2 -M2 -v snps \
            -O z -o "$OUT_DIR/chr${c}.vcf.gz" \
            "${URL_BASE}.chr${c}.vcf.bgz" 2>/dev/null; then
        echo "[$(date +%H:%M:%S)] done chr${c} $(du -sh $OUT_DIR/chr${c}.vcf.gz | cut -f1)" >> "$LOG"
    else
        echo "[$(date +%H:%M:%S)] FAIL chr${c}" >> "$LOG"
    fi
}

for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]; do sleep 5; done
    filter_chrom "$c" &
done
wait
echo "[$(date +%H:%M:%S)] all chromosomes done" >> "$LOG"

echo "[$(date +%H:%M:%S)] merging" >> "$LOG"
bcftools concat -O z -o "$OUT_DIR/${OUT_NAME}.vcf.gz" "$OUT_DIR"/chr{1..22}.vcf.gz 2>>"$LOG" \
    && bcftools index -t "$OUT_DIR/${OUT_NAME}.vcf.gz" 2>>"$LOG" \
    && echo "[$(date +%H:%M:%S)] merge done $(du -sh $OUT_DIR/${OUT_NAME}.vcf.gz | cut -f1)" >> "$LOG" \
    || echo "[$(date +%H:%M:%S)] merge FAILED" >> "$LOG"

echo "[$(date +%H:%M:%S)] ALL_DONE" >> "$LOG"
echo "Build complete. Master DB at: $OUT_DIR/${OUT_NAME}.vcf.gz"
