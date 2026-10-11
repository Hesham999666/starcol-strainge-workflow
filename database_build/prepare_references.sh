#!/usr/bin/env bash

set -euo pipefail

# Reconstructed public reference-preparation workflow.
# The original project-specific fasta_adapter helper was not retained.
#
# This script:
#   1. reads draft genome assemblies in FASTA format
#   2. assigns unique contig headers based on the assembly filename
#   3. writes each sequence on a single line
#   4. compresses the prepared reference as .fa.gz
#
# Input example:
#   assemblies/Bifidobacterium_longum_37.fasta
#
# Output example:
#   prepared/Bifidobacterium_longum_37.fa.gz
#
# Contig headers become:
#   >Bifidobacterium_longum_37_c1
#   >Bifidobacterium_longum_37_c2
#   ...

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <assembly_directory> <output_directory>"
    exit 1
fi

INPUT_DIR="$(realpath "$1")"
OUTPUT_DIR="$2"

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"

mapfile -t FASTAS < <(
    find "$INPUT_DIR" -maxdepth 1 -type f \
    \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) \
    | sort
)

if [[ ${#FASTAS[@]} -eq 0 ]]; then
    echo "ERROR: No FASTA assemblies found in $INPUT_DIR"
    exit 1
fi

echo "Assemblies found: ${#FASTAS[@]}"

for fasta in "${FASTAS[@]}"
do
    filename="$(basename "$fasta")"
    sample="${filename%.fasta}"
    sample="${sample%.fna}"
    sample="${sample%.fa}"

    output="${OUTPUT_DIR}/${sample}.fa.gz"

    echo "Preparing: $sample"

    awk -v prefix="$sample" '
        function flush_sequence() {
            if (seq != "") {
                print seq
                seq = ""
            }
        }

        /^>/ {
            flush_sequence()
            n++
            print ">" prefix "_c" n
            next
        }

        {
            gsub(/[[:space:]]/, "", $0)
            seq = seq $0
        }

        END {
            flush_sequence()
        }
    ' "$fasta" | gzip -c > "$output"
done

echo
echo "Prepared references written to:"
echo "$OUTPUT_DIR"
