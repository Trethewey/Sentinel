#!/usr/bin/env bash
# DNA Nexus entrypoint for Sentinel.
# Stages inputs, runs `sentinel run`, uploads outputs.
set -euxo pipefail

main() {
    # Pull inputs to local staging
    mkdir -p inputs/bams outputs

    dx download "$sample_sheet" -o inputs/samples.csv
    dx download "$catalog" -o inputs/catalog.tsv
    if [ -n "${reference_fasta:-}" ]; then
        dx download "$reference_fasta" -o inputs/ref.fa
        REF_ARG="--ref inputs/ref.fa"
    else
        REF_ARG=""
    fi

    # Download all BAMs and their indexes (DNA Nexus passes arrays as $bams[@])
    for fid in "${bams[@]}"; do dx download "$fid" -o inputs/bams/; done
    for fid in "${bam_indexes[@]}"; do dx download "$fid" -o inputs/bams/; done

    STRIP=""
    if [ "${strip_chr_prefix:-false}" = "true" ]; then STRIP="--strip-chr-prefix"; fi

    # Run via the Docker image baked into this app's resources, or assume sentinel is installed
    sentinel run \
        --sample-sheet inputs/samples.csv \
        --bams-dir inputs/bams \
        --catalog inputs/catalog.tsv \
        ${REF_ARG} \
        --out outputs \
        --threads "${threads:-16}" \
        ${STRIP}

    # Upload outputs
    report_html=$(dx upload outputs/results/report.html --brief)
    per_sample_csv=$(dx upload outputs/results/per_sample_report.csv --brief)
    per_pair_csv=$(dx upload outputs/results/scores/per_pair_report.tsv --brief)

    dx-jobutil-add-output report_html    "$report_html"    --class=file
    dx-jobutil-add-output per_sample_csv "$per_sample_csv" --class=file
    dx-jobutil-add-output per_pair_csv   "$per_pair_csv"   --class=file
}
