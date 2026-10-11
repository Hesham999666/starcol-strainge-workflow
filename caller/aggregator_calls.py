#!/usr/bin/env python3
"""
Version 2025-12-17a

Purpose:
--------
This script processes microbiome strain-level summary tables by merging:

1. A `<sample>.called.tsv` file containing contig-level read mapping metrics
2. A `<sample>.strains.tsv` file containing strain-level k-mer metadata

! Note that Nextflow will do the concatenation of all files produced by this script.

The script:
- Filters contigs based on minimum length and coverage
- Aggregates metrics for each reference genome
- Computes weighted means
- Computes percentages such as callable%, confirmed%, etc.
- Outputs a clean, final TSV file for downstream analysis
"""

# -------------------------------------------------
# Standard library imports
# -------------------------------------------------
from pathlib import Path      # Convenient and safer filesystem paths
import argparse               # For command-line argument parsing

# -------------------------------------------------
# Third-party imports
# -------------------------------------------------
import pandas as pd           # Main table handling library
import numpy as np            # Used for weighted calculations (optional but common)

# -------------------------------------------------
# Command-line parser
# -------------------------------------------------
parser = argparse.ArgumentParser(
    description="Process strain and called TSV files and generate summarized output."
)

# Required input files
parser.add_argument(
    "-strains", required=True,
    help="Path to the input .strains.tsv file"
)

parser.add_argument(
    "-called", required=True,
    help="Path to the input .called.tsv file (contig-level metrics)"
)

parser.add_argument(
    "-reads", required=True,
    help="Path to the input sample_id_reads_in_fastqs.txt file (total reads in fastq.gz files)"
)

# Output path
parser.add_argument(
    "-o", "--output", required=True,
    help="Where to write the merged output TSV file"
)

args = parser.parse_args()

# Convert into Path objects for safe filesystem handling
strains_path = Path(args.strains)
called_path = Path(args.called)
total_reads_path = Path(args.reads)
output_path = Path(args.output)

# -------------------------------------------------
# 1. Read the called.tsv file
# -------------------------------------------------
df = pd.read_table(called_path, sep="\t")

# -------------------------------------------------
# Read total reads from FASTQ count file
# -------------------------------------------------
with open(total_reads_path) as f:
    total_reads_fastq = float(f.read().strip())

# -------------------------------------------------
# 2. Compute unfiltered total length per strain
# -------------------------------------------------
# Before any filtering, we need the sum of all contig lengths per strain (ref)
unfiltered_length_sum = (
    df.groupby("ref")["length"]
      .sum()
      .rename("length_unfiltered")
)

# -------------------------------------------------
# 3. Apply filtering rules
# -------------------------------------------------
# The user keeps only:
#   - contigs >= 1000 bp
#   - contigs with coverage > 0.5
filter1_df = df[df["length"] >= 1000]
filter2_df = filter1_df[filter1_df["coverage"] > 0.5]

# -------------------------------------------------
# 4. Helper function: pick the most common contig "name"
# -------------------------------------------------
def most_frequent_name(series):
    """
    Returns the most frequently occurring value in a pandas Series.
    .mode() returns a sorted list of the most frequent values.
    Using .iloc[0] ensures we pick one even when ties occur.
    """
    return series.mode().iloc[0]

# -------------------------------------------------
# 5. Group and summarize filtered contigs
# -------------------------------------------------
# For each strain (ref), sum or average appropriate columns.
grouped_df = filter2_df.groupby("ref").agg({
    "name": most_frequent_name,   # Most common assigned name
    "length": "sum",              # Sum of kept contig lengths
    "coverage": "mean",           # TEMP placeholder — corrected later
    "uReads": "sum",
    "abundance": "mean",          # TEMP placeholder — corrected later
    "median": "mean",             # TEMP placeholder — corrected later
    "callable": "sum",
    "confirmed": "sum",
    "snps": "sum",
    "multi": "sum",
    "lowmq": "sum",
    "high": "sum",
    "gapCount": "sum",
    "gapLength": "sum"
})

# Attach unfiltered lengths
grouped_df = grouped_df.merge(
    unfiltered_length_sum,
    on="ref",
    how="left"
)

# -------------------------------------------------
# 6. Compute weighted means for coverage/abundance/median
# -------------------------------------------------
# Weighted by contig length, which is more biologically meaningful.
weightedCoverage = (
    (filter2_df["coverage"] * filter2_df["length"])
    .groupby(filter2_df["ref"]).sum()
)

weightedAbundance = (
    (filter2_df["abundance"] * filter2_df["length"])
    .groupby(filter2_df["ref"]).sum()
)

strain_median = (
    filter2_df
    .groupby("ref")["coverage"]
    .median()
)
grouped_df["medianCoverage"] = strain_median

# Replace placeholder means with weighted values
grouped_df["coverage"] = weightedCoverage / grouped_df["length"]
grouped_df["abundance"] = weightedAbundance / grouped_df["length"]
grouped_df["median"] = strain_median

# -------------------------------------------------
# 6b. Explicit length-weighted genome coverage
# -------------------------------------------------
grouped_df["weightedCoverage"] = grouped_df["coverage"]

# -------------------------------------------------
# 6c. Relative abundance among detected reference strains
# -------------------------------------------------
total_weighted_coverage = grouped_df["weightedCoverage"].sum()

grouped_df["abundance_mapped"] = (
    grouped_df["weightedCoverage"] / total_weighted_coverage
)

# -------------------------------------------------
# 6d. Fraction of metagenome covered by reference strains (read-based)
# -------------------------------------------------
mapped_reads_total = grouped_df["uReads"].sum()

fraction_metagenome_mapped = mapped_reads_total / total_reads_fastq

# -------------------------------------------------
# 6e. Abundance relative to the whole metagenome
# -------------------------------------------------
grouped_df["abundance_metagenome"] = (
    grouped_df["abundance_mapped"] * fraction_metagenome_mapped
)

# -------------------------------------------------
# 7. Add sample name (derived from filename)
# -------------------------------------------------
# A file called "sample001.called.tsv" becomes sample name "sample001"
sample_name = called_path.stem.replace(".called", "")
grouped_df["samplename"] = sample_name

# -------------------------------------------------
# 8. Recalculate percentages
# -------------------------------------------------
# Protect against division by zero using .replace(0, np.nan)

grouped_df["callablePct"] = (
    grouped_df["callable"] / grouped_df["length_unfiltered"]
)

grouped_df["confirmedPct"] = (
    grouped_df["confirmed"] / grouped_df["callable"].replace(0, np.nan)
)

grouped_df["snpPct"] = (
    grouped_df["snps"] / grouped_df["callable"].replace(0, np.nan)
)

grouped_df["multiPct"] = (
    grouped_df["multi"] / grouped_df["callable"].replace(0, np.nan)
)

grouped_df["lowmqPct"] = (
    grouped_df["lowmq"] / grouped_df["length_unfiltered"]
)

grouped_df["highPct"] = (
    grouped_df["high"] / grouped_df["length_unfiltered"]
)

# -------------------------------------------------
# 9. Reorder output columns
# -------------------------------------------------
new_order = [
    "samplename", "name",
    "length_unfiltered", "length",
    "weightedCoverage",
    "uReads",
    "abundance_mapped",
    "abundance_metagenome",
    "median",
    "callable", "callablePct",
    "confirmed", "confirmedPct",
    "snps", "snpPct",
    "multi", "multiPct",
    "lowmq", "lowmqPct",
    "high", "highPct",
    "gapCount", "gapLength"
]

df_new_order = grouped_df[new_order].reset_index()

# Rename column for clear output
df_renamed = df_new_order.rename(columns={"ref": "strain"})

# -------------------------------------------------
# 10. Merge with the strains.tsv file
# -------------------------------------------------
strainsTable = pd.read_csv(strains_path, sep="\t")

df_merged = pd.merge(
    strainsTable,
    df_renamed,
    on="strain",
    how="left"
)

# -------------------------------------------------
# 11. Drop deprecated or legacy columns if present
# -------------------------------------------------
deprecated_cols = ["wscore", "spec", "old_rapct", "i"]

for col in deprecated_cols:
    if col in df_merged.columns:
        df_merged = df_merged.drop(columns=[col])

# -------------------------------------------------
# 12. Write the final output file
# -------------------------------------------------
df_merged.to_csv(output_path, sep="\t", index=False)

print(f"\n Output written to: {output_path}")
print(df_merged)
