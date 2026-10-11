#!/usr/bin/env python3
"""
Version 2025-12-03a

Process a summarized comparison .summary.tsv file and generate
two scatterplots showing gap similarity vs. ACNI (singleAgreePct).

This script is designed to be beginner-friendly, with extensive comments
explaining each step, why it is done, and how to modify it.
"""

# ===============================================================
# Imports
# ===============================================================
# Built-in modules
import re                     # For regular expressions used in renaming
import argparse               # For handling command-line arguments
from pathlib import Path      # For filesystem-safe path handling
import sys                    # Only used for edge-case exits (not required here)

# Third-party modules
import pandas as pd           # For loading/manipulating tabular data
import numpy as np            # Often useful in numerical code (not heavily used here)
import seaborn as sns         # Statistical plotting
import matplotlib.pyplot as plt
import matplotlib.patches as patches  # For drawing the highlight box

# ===============================================================
# Command-line argument parsing
# ===============================================================
# Example usage:
#     python script.py -compSummed input.summary.tsv -o output_plot
#
parser = argparse.ArgumentParser(
    description=(
        "Process summarized comparison of HDF5 comparison TSV files "
        "and generate scatterplots (PNG + PDF)."
    )
)

parser.add_argument(
    "-compSummed",
    required=True,
    help="Path to the input .summary.tsv file"
)

parser.add_argument(
    "-o", "--output",
    required=True,
    help=(
        "Output prefix for plot files. "
        "The script will generate <output>.png and <output>.pdf"
    )
)

# Parse arguments
args = parser.parse_args()

# Convert to Path objects (safe, recommended)
compared_path = Path(args.compSummed)
output_path = Path(args.output)

# ===============================================================
# Load Input Data
# ===============================================================
# pandas automatically handles headers and tab-separated files using sep="\t"
df = pd.read_table(compared_path, sep="\t")

# Convert fraction to percent for readability
df["singleAgreePct"] = df["singleAgreePct"] * 100


# ===============================================================
# Formatting helper functions
# ===============================================================

def shorten_species(name: str) -> str:
    """
    Convert a scaffold name such as:
         Bacteroides_uniformis_224_c1
    into a short italic scientific name:
         B. uniformis
    Uses LaTeX formatting so matplotlib renders italics properly.
    """
    parts = name.split('_')

    if len(parts) < 2:
        return name  # Fallback, unexpected case

    genus = parts[0]
    species = parts[1]

    # Abbreviate genus (e.g. Bacteroides → B.)
    short_genus = genus[0] + "."

    # Return species name in italic
    return rf"$\mathit{{{short_genus}}}\ \mathit{{{species}}}$"


def full_name(name: str) -> str:
    """
    Expand scaffold name into:
        Genus species XXX
    Mainly kept for reference; not used in the main plots.
    """
    parts = name.split('_')

    if len(parts) < 3:
        return name

    genus = parts[0]
    species = parts[1]
    number = parts[2]

    return f"{genus} {species} {number}"


def shorten_and_format_label(name: str) -> str:
    """
    Format a scaffold name into:
         B. uniformis (strain number non-italic)
    Removes trailing _c###, which often refers to contigs.

    Returns string with LaTeX formatting for:
      - italic genus+species
      - normal (roman) strain numbers
    """
    # Remove trailing "_c<number>" if present
    name = re.sub(r"_c\d+$", "", name)

    parts = name.split("_")
    if len(parts) < 2:
        return name

    genus = parts[0]
    species = parts[1]
    rest = " ".join(parts[2:])  # strain or extra ID components

    genus_abbrev = f"{genus[0]}."   # First letter of genus + dot

    italic_part = f"{genus_abbrev} {species}"

    if rest:
        return rf"$\mathit{{{italic_part}}}\ \mathrm{{{rest}}}$"
    else:
        return rf"$\mathit{{{italic_part}}}$"


# ===============================================================
# Apply transformations
# ===============================================================
df["Species"] = df["scaffold"].apply(shorten_species)
df["Strains"] = df["scaffold"].apply(shorten_and_format_label)

# ===============================================================
# Create two side-by-side plots
# ===============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 13))

# ============================================================
# Prepare color palette
# ============================================================
# Get unique species in your dataset
unique_species = df["Species"].unique()

# Create combined palette
combined_palette = sns.color_palette("tab20") + sns.color_palette("tab10")

# Trim to exactly the number of categories
colors = combined_palette[:len(unique_species)]

# ---------------------------------------------------------------
# FIRST PLOT — Full range
# ---------------------------------------------------------------
sns.scatterplot(
    x="gapJaccardSim",
    y="singleAgreePct",
    size="commonPct",
    sizes=(20, 140),
    legend=False,            # hide legend here
    style="Strains",
    hue="Species",
    alpha=0.8,
    edgecolor="w",
    linewidth=0.5,
    data=df,
    ax=ax1
)
ax1.set_xlim(0.90, 1.002)
ax1.set_ylim(97, 100.1)
ax1.set_xlabel("Gap similarity (Jaccard index)")
ax1.set_ylabel("ACNI (single-agree in %)")
ax1.grid(True)
ax1.set_title("Wide range (ACNI 97–100%)")

# ---------------------------------------------------------------
# Zoom box
# ---------------------------------------------------------------
highlight_box = patches.Rectangle(
    xy=(0.970, 99.9),
    width=1.000 - 0.970,
    height=100 - 99.9,
    linewidth=1.5,
    edgecolor="gray",
    facecolor="gray",
    alpha=0.10
)
ax1.add_patch(highlight_box)

# ---------------------------------------------------------------
# SECOND PLOT — Zoomed region
# ---------------------------------------------------------------
sns.scatterplot(
    x="gapJaccardSim",
    y="singleAgreePct",
    size="commonPct",
    sizes=(20, 140),
    legend="brief",          # let seaborn create legend handles
    style="Strains",
    hue="Species",
    alpha=0.8,
    edgecolor="w",
    linewidth=0.5,
    data=df,
    ax=ax2
)
ax2.set_xlim(0.970, 1.0007)
ax2.set_ylim(99.9, 100.002)
ax2.set_xlabel("Gap similarity (Jaccard index)")
ax2.set_ylabel("ACNI (single-agree in %)")
ax2.grid(True)
ax2.set_title("Zoomed View (ACNI 99.9–100%)")

# ---------------------------------------------------------------
# Create a single legend below both plots
# ---------------------------------------------------------------
handles, labels = ax2.get_legend_handles_labels()

# Move legend outside plot area
ax2.legend(
    title="Legend",
    loc="upper left",
    bbox_to_anchor=(1.01, 1.0),
    frameon=True
)



# ===============================================================
# Save output
# ===============================================================
png_path = output_path.with_suffix(".png")
pdf_path = output_path.with_suffix(".pdf")

plt.tight_layout()
plt.savefig(png_path, dpi=300, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
