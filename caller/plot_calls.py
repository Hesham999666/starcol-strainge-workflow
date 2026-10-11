#!/usr/bin/env python3
"""
Version 2025-12-17a

Process summarized k-mer comparison data and produce stacked bar plots showing:

1. Relative abundance based on k-mer RAPCT (ax1)
2. Relative abundance based on read abundance (ax2)

This script:
- Reads a TSV file containing microbiome strain comparison data
- Cleans sample names and strain labels
- Aggregates abundance metrics per sample
- Produces stacked bar charts for visualization
- Saves plots as PNG and PDF
"""

# ============================================================
# Import required libraries
# ============================================================
# re               = regular expressions for cleaning names
# pathlib.Path     = modern, OS-safe path handling
# argparse         = command-line argument parsing
# pandas           = table/column processing
# numpy            = numerical operations
# matplotlib       = plotting
# sys              = system-level error handling (optional)
# sns              = Statistical plotting
# ============================================================

import re
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys


# ============================================================
# 1. Parse command-line arguments
# ============================================================
parser = argparse.ArgumentParser(
    description=(
        "Process summarized comparison of hdf5-derived TSV files "
        "and generate stacked bar plots for abundance metrics."
    )
)

# Input TSV file
parser.add_argument(
    "-callsSummed",
    required=True,
    help="Path to the input .summary.tsv file containing merged metrics"
)

# Output prefix (PNG + PDF will be written)
parser.add_argument(
    "-o", "--output",
    required=True,
    help="Output filename prefix (no extension needed)"
)

args = parser.parse_args()

# Convert into Path objects for OS-safe use
callsSummed_path = Path(args.callsSummed)
output_path = Path(args.output)


# ============================================================
# 2. Load the input TSV file
# ============================================================
try:
    df = pd.read_table(callsSummed_path, sep="\t")
except Exception as e:
    sys.exit(f"ERROR: Could not read input file: {e}")

# Convert fraction to percent for readability
df["abundance"] = df["abundance_metagenome"] * 100


# ============================================================
# 3. Drop rows missing essential fields
# ============================================================
# Many microbiome datasets contain incomplete rows.
# samplename and name are required for grouping.
df = df.dropna(subset=['samplename', 'name'])


# ============================================================
# 4. Clean sample names — adjust dataset labels
# ============================================================
def adjust_name(name) -> str:
    """
    Safely clean sample names.
    Why this function exists:
        - Some samples contain donor tags (Donor_A, Donor_B...)
        - Others have extra suffixes such as "_1A" or "_MG123"
        - For non-strings or NaN, we return unchanged to avoid errors

    Parameters:
        name (str): A sample name

    Returns:
        str: A cleaned sample name
    """

    # These keywords mean we should NOT strip anything off
    donor_keywords = ("Donor_A", "Donor_B", "Donor_C")

    # If name is not a string (float NA), return it unchanged
    if not isinstance(name, str):
        return name

    # Keep donor names untouched
    if any(k in name for k in donor_keywords):
        return name

    # Otherwise split on "_" and take only the first part
    # Example: "sample001_run1" → "sample001"
    return name.split("_")[0]


# Convert to string dtype (pandas nullable string type)
df['samplename'] = df['samplename'].astype("string")

# Apply cleaning function
df['samplename'] = df['samplename'].apply(adjust_name)


# ============================================================
# 5. Clean strain names for plotting (italicized scientific names)
# ============================================================
def shorten_and_format_label(name: str) -> str:
    """
    Convert strain names into shortened scientific notation for plots.

    Example:
        "Bacteroides_uniformis_111_c1"
        → "$\\mathit{B. uniformis}$ 111"

    This improves readability on x/y-axis and legend labels.

    Steps:
        1. Remove "_c1", "_c2" style suffixes
        2. Split by "_"
        3. Abbreviate genus (Bacteroides → B.)
        4. Format for matplotlib italic/roman style

    Returns:
        A TeX/LaTeX-style string for matplotlib (math mode).
    """

    # Handle unexpected non-string input
    if not isinstance(name, str):
        return str(name)

    # Remove trailing "_c1" or "_c12"
    name = re.sub(r"_c\d+$", "", name)

    parts = name.split("_")

    # Expect at least genus + species
    if len(parts) < 2:
        return name

    genus = parts[0]
    species = parts[1]

    # Remaining parts often contain strain codes
    rest = " ".join(parts[2:])

    # Abbreviate genus ("Bacteroides" → "B.")
    genus_abbrev = f"{genus[0]}."

    # Italicize genus + species
    italic_part = f"{genus_abbrev} {species}"

    # Final formatted label
    if rest:
        return rf"$\mathit{{{italic_part}}}\ \mathrm{{{rest}}}$"
    else:
        return rf"$\mathit{{{italic_part}}}$"


# Apply strain name cleaning
df["name"] = df["name"].apply(shorten_and_format_label)


# ============================================================
# 6. Get unique sample names and strain names in consistent order
# ============================================================
# We want a fixed sample order left-to-right on the plot
sample_names = sorted(df['samplename'].unique())

# Strain order is preserved as they appear in the data
strain_order = df['name'].unique()

# x-axis positions
x = np.arange(len(sample_names))


# ============================================================
# 7. Prepare strain labels and color palette
# ============================================================
formatted_labels = [shorten_and_format_label(s) for s in strain_order]

# Choose a colormap: plt.cm.<color map>()
# colors = sns.color_palette("husl", len(strain_order))
colors = sns.color_palette("tab20") + sns.color_palette("tab10")

# ============================================================
# 8. Create the figure layout
# ============================================================
fig, (ax1, ax2) = plt.subplots(
    1, 2,
    figsize=(20, 11),
    sharex=True  # Both plots share the same x-axis labels
)


# ============================================================
# 9. Plot 1 → RAPCT (k-mer relative abundance)
# ============================================================
bottom = np.zeros(len(sample_names))  # tracks stacked bar bottoms

for (strain_raw, label, color) in zip(strain_order, formatted_labels, colors):

    # Compute abundance for each sample (sum handles duplicates safely)
    rel_abun = np.array([
        df.loc[
            (df['samplename'] == sample) &
            (df['name'] == strain_raw),
            'rapct'
        ].sum()
        for sample in sample_names
    ])

    # Skip strains not present in any sample
    if np.allclose(rel_abun, 0.0):
        continue

    # Create stacked bar layer
    ax1.bar(x, rel_abun, bottom=bottom, label=label, color=color, width=0.8)

    # Update bottom for next strain
    bottom += rel_abun

ax1.set_ylabel("Relative abundance (in %)")
ax1.set_xlabel("Sample")
ax1.set_title("Abundance - K-mer RAPCT")
ax1.set_xticks(x)
ax1.set_xticklabels(sample_names, rotation=45, ha='right')


# ============================================================
# 10. Plot 2 → ABUNDANCE (read-based relative abundance)
# ============================================================
bottom = np.zeros(len(sample_names))  # reset for second plot

for (strain_raw, label, color) in zip(strain_order, formatted_labels, colors):

    rel_abun = np.array([
        df.loc[
            (df['samplename'] == sample) &
            (df['name'] == strain_raw),
            'abundance'
        ].sum()
        for sample in sample_names
    ])

    if np.allclose(rel_abun, 0.0):
        continue

    ax2.bar(x, rel_abun, bottom=bottom, label=label, color=color, width=0.8)
    bottom += rel_abun

ax2.set_ylabel("Relative abundance (in %)")
ax2.set_xlabel("Sample")
ax2.set_title("Abundance - Reads")
ax2.set_xticks(x)
ax2.set_xticklabels(sample_names, rotation=45, ha='right')


# ============================================================
# 11. Shared legend (placed outside the figure)
# ============================================================
handles, labels = ax2.get_legend_handles_labels()

fig.legend(
    handles[::-1], labels[::-1],   # reverse for consistent ordering
    title="Strains",
    bbox_to_anchor=(1.0001, 0.95),
    loc="upper left"
)


# ============================================================
# 12. Save the figure
# ============================================================
png_path = output_path.with_suffix(".png")
pdf_path = output_path.with_suffix(".pdf")

plt.tight_layout()
plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
