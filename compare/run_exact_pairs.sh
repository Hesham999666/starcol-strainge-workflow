#!/usr/bin/env bash

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PAIR_TABLE="${1:-${REPO_DIR}/example_metadata/exact_donor_pairs.example.tsv}"

INPUT_DIR="${REPO_DIR}/data/hdf5"
OUTPUT_ROOT="${REPO_DIR}/results/exact_pairs"

mkdir -p "$INPUT_DIR" "$OUTPUT_ROOT"

while IFS=$'\t' read -r \
    comparison_id \
    recipient_id \
    recipient_timepoint \
    donor_id \
    donor_hdf5 \
    recipient_hdf5
do
    [[ "$comparison_id" == "comparison_id" ]] && continue
    [[ -n "$comparison_id" ]] || continue

    echo
    echo "=================================================="
    echo "Comparison: $comparison_id"
    echo "Recipient: $recipient_id"
    echo "Timepoint: $recipient_timepoint"
    echo "Exact donor: $donor_id"
    echo "=================================================="

    pair_output="${OUTPUT_ROOT}/${comparison_id}"

    rm -rf "$INPUT_DIR"
    mkdir -p "$INPUT_DIR"

    cp "$donor_hdf5" "$INPUT_DIR/"
    cp "$recipient_hdf5" "$INPUT_DIR/"

    rm -rf "${REPO_DIR}/results/ComparisonSummed" \
           "${REPO_DIR}/results/Hdf5Compared" \
           "${REPO_DIR}/results/ComparisonPlots"

    cd "$REPO_DIR"

    nextflow run compare/main_strainGE_compare.nf \
        -profile saga \
        -resume

    mkdir -p "$pair_output"

    cp -a results/ComparisonSummed/. "$pair_output/" 2>/dev/null || true
    cp -a results/Hdf5Compared/. "$pair_output/" 2>/dev/null || true

    printf "comparison_id\trecipient_id\trecipient_timepoint\tdonor_id\n" \
        > "${pair_output}/comparison_metadata.tsv"

    printf "%s\t%s\t%s\t%s\n" \
        "$comparison_id" \
        "$recipient_id" \
        "$recipient_timepoint" \
        "$donor_id" \
        >> "${pair_output}/comparison_metadata.tsv"

done < "$PAIR_TABLE"

echo
echo "Exact-donor comparison workflow finished."
