#!/usr/bin/env python3
"""
Version 2025-12-02a

Summarize contig-level comparison statistics from an input TSV file.

This script:
1. Reads an input .summary.tsv file (per-contig comparison output)
2. Filters out contigs shorter than 1000 bp
3. Groups remaining rows by the reference genome / contig ("ref")
4. Computes:
    - Sums of all count-like metrics
    - Weighted mean of gapJaccardSim (weighted by contig length)
    - Recomputed percentage columns
5. Returns one summarized row per reference
"""

# -------------------------------------------------
# Imports
# -------------------------------------------------
from pathlib import Path
import pandas as pd
import argparse
import sys


# -------------------------------------------------
# Parse command-line arguments
# -------------------------------------------------
parser = argparse.ArgumentParser(
    description="Process summarized comparison of HDF5 files (TSV) and generate a per-reference summary."
)

# Required arguments:
parser.add_argument(
    "-compared",
    required=True,
    help="Path to the contig-level .summary.tsv file."
)
parser.add_argument(
    "-o", "--output",
    required=True,
    help="Path to the output summarized .tsv file."
)

args = parser.parse_args()


# -------------------------------------------------
# Convert inputs to pathlib.Path objects
# -------------------------------------------------
compared_path = Path(args.compared)
output_path = Path(args.output)

# Check input existence early
if not compared_path.exists():
    sys.exit(f"ERROR: Input file not found: {compared_path}")


# -------------------------------------------------
# Load input data
# -------------------------------------------------
"""
pd.read_table() reads tab-separated files.
If your input contains many rows, Pandas automatically uses efficient parsers.
"""
df = pd.read_table(compared_path, sep="\t")

# Basic safety check
required_columns = {
    'sample1','sample2','ref','scaffold','length',
    'common','single','singleAgree','multi','sharedAlleles',
    'variants','commonVariant','variantExact','AnotB','BnotA',
    'Agaps','Bgaps','gapJaccardSim'
}

missing = required_columns - set(df.columns)
if missing:
    sys.exit(f"ERROR: Missing required columns: {', '.join(missing)}")


# -------------------------------------------------
# Filter out small contigs
# -------------------------------------------------
"""
Why filter? Very small contigs can distort weighted averages and
artificially inflate percentages, especially when allele counts are tiny.
"""
filtered_df = df[df['length'] >= 1000]


# -------------------------------------------------
# Helper functions for selecting the "most frequent" labels
# -------------------------------------------------
"""
These functions return the most common value in a column (mode).
`.mode()` returns a Series; `.iloc[0]` extracts the first (safe for ties).
"""

def most_frequent(x):
    """Return the most common value in a Series."""
    return x.mode().iloc[0]


# -------------------------------------------------
# Group by reference and aggregate statistics
# -------------------------------------------------
"""
We group by 'ref' because each row represents a contig.
After grouping, each row in grouped_df summarises all contigs
belonging to the same reference (strain/assembly).

For categorical columns: take most frequent
For numerical count columns: sum them
For gapJaccardSim: we compute weighted mean manually
"""
grouped_df = filtered_df.groupby('ref').agg({
    'sample1': most_frequent,
    'sample2': most_frequent,
    'ref': most_frequent,
    'scaffold': most_frequent,   # Often identical, but included for completeness

    # Sum count-like columns
    'length': 'sum',
    'common': 'sum',
    'single': 'sum',
    'singleAgree': 'sum',
    'multi': 'sum',
    'sharedAlleles': 'sum',
    'variants': 'sum',
    'commonVariant': 'sum',
    'variantExact': 'sum',
    'AnotB': 'sum',
    'BnotA': 'sum',
    'Agaps': 'sum',
    'Bgaps': 'sum',

    # Temporary mean—overwritten later by weighted mean
    'gapJaccardSim': 'mean',
})


# -------------------------------------------------
# Compute weighted mean of gapJaccardSim
# -------------------------------------------------
"""
Weighted mean formula:
    Σ(similarity_i * length_i) / Σ(length_i)

We compute the weighted sum per reference,
then divide by total length per reference (already in grouped_df)
"""
weighted_sum = (
    filtered_df['gapJaccardSim'] * filtered_df['length']
).groupby(filtered_df['ref']).sum()

grouped_df['gapJaccardSim'] = weighted_sum / grouped_df['length']


# -------------------------------------------------
# Recompute percentage fields
# -------------------------------------------------
"""
Important: After summing counts, percentages must be recomputed.
We follow the logic from the original columns.

Some divisions may risk division-by-zero for strange data,
so we use .where(divisor != 0, 0) to avoid errors.
"""

def safe_divide(numer, denom):
    return (numer / denom).where(denom != 0, 0)

grouped_df['commonPct']        = safe_divide(grouped_df['common'],        grouped_df['length'])
grouped_df['singlePct']        = safe_divide(grouped_df['single'],        grouped_df['common'])
grouped_df['singleAgreePct']   = safe_divide(grouped_df['singleAgree'],   grouped_df['single'])
grouped_df['multiPct']         = safe_divide(grouped_df['multi'],         grouped_df['common'])
grouped_df['sharedAllelesPct'] = safe_divide(grouped_df['sharedAlleles'], grouped_df['common'])
grouped_df['variantPct']       = safe_divide(grouped_df['variants'],      grouped_df['common'])
grouped_df['commonVariantPct'] = safe_divide(grouped_df['commonVariant'], grouped_df['variants'])
grouped_df['variantExactPct']  = safe_divide(grouped_df['variantExact'],  grouped_df['variants'])
grouped_df['AnotBpct']         = safe_divide(grouped_df['AnotB'],         grouped_df['variants'])
grouped_df['BnotApct']         = safe_divide(grouped_df['BnotA'],         grouped_df['variants'])


# -------------------------------------------------
# Reorder columns for final output
# -------------------------------------------------
column_order = [
    'sample1', 'sample2', 'ref', 'scaffold', 'length',
    'common', 'commonPct',
    'single', 'singlePct',
    'singleAgree', 'singleAgreePct',
    'multi', 'multiPct',
    'sharedAlleles', 'sharedAllelesPct',
    'variants', 'variantPct',
    'commonVariant', 'commonVariantPct',
    'variantExact', 'variantExactPct',
    'AnotB', 'AnotBpct',
    'BnotA', 'BnotApct',
    'Agaps', 'Bgaps',
    'gapJaccardSim'
]

df_final = grouped_df[column_order]


# -------------------------------------------------
# Write output TSV
# -------------------------------------------------
df_final.to_csv(output_path, sep="\t", index=False)

print(f"\n Output written to: {output_path}\n")
print(df_final)

