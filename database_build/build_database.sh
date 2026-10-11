#!/usr/bin/env bash

set -euo pipefail

# Reconstructed public database-build workflow.
# The original project-specific helper script was not retained.
#
# Input:
#   A directory containing prepared reference assemblies as *.fa.gz
#
# Output:
#   Individual StrainGST HDF5 k-mer sets
#   All-vs-all similarity table
#   Pan-genome StrainGST database containing all supplied references

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <reference_fasta_directory> [output_prefix]"
    exit 1
fi

REF_DIR="$(realpath "$1")"
PREFIX="${2:-DonorB}"

echo "Reference directory: $REF_DIR"
echo "Database prefix:     $PREFIX"

mapfile -t FASTAS < <(find "$REF_DIR" -maxdepth 1 -type f -name '*.fa.gz' | sort)

if [[ ${#FASTAS[@]} -eq 0 ]]; then
    echo "ERROR: No .fa.gz reference files found in $REF_DIR"
    exit 1
fi

echo "Reference FASTA files found: ${#FASTAS[@]}"

# ------------------------------------------------------------
# 1. Convert each reference assembly to a StrainGST k-mer set
# ------------------------------------------------------------
for fasta in "${FASTAS[@]}"
do
    hdf5="${fasta%.fa.gz}.hdf5"

    echo "Kmerizing: $(basename "$fasta")"

    straingst kmerize \
        -k 23 \
        -o "$hdf5" \
        "$fasta"
done

mapfile -t HDF5S < <(find "$REF_DIR" -maxdepth 1 -type f -name '*.hdf5' | sort)

echo "HDF5 k-mer sets found: ${#HDF5S[@]}"

# ------------------------------------------------------------
# 2. Calculate all-vs-all reference similarities
# ------------------------------------------------------------
straingst kmersim \
    --all-vs-all \
    -t 4 \
    -S jaccard \
    -S subset \
    "${HDF5S[@]}" \
    > "${PREFIX}_similarities.tsv"

# ------------------------------------------------------------
# 3. Retain all reference k-mer sets
# ------------------------------------------------------------
printf "%s\n" "${HDF5S[@]}" > "${PREFIX}_references.txt"

# ------------------------------------------------------------
# 4. Create the StrainGST pan-genome database
# ------------------------------------------------------------
straingst createdb \
    -f "${PREFIX}_references.txt" \
    -o "${PREFIX}_pan-genome-db.hdf5"

echo
echo "Database build complete."
echo "Similarity table: ${PREFIX}_similarities.tsv"
echo "Reference list:   ${PREFIX}_references.txt"
echo "Pan-genome DB:    ${PREFIX}_pan-genome-db.hdf5"
