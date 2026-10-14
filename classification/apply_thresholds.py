#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


CALLER_REQUIRED = [
    "abundance_metagenome",
    "callablePct",
    "confirmedPct",
    "snpPct",
    "gapCount",
    "gapLength",
]

COMPARE_REQUIRED = [
    "commonPct",
    "singleAgreePct",
]


def require_columns(df, required, label):
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise SystemExit(
            f"ERROR: {label} table is missing required columns: "
            + ", ".join(missing)
        )


def classify_caller(path, outdir):
    df = pd.read_csv(path, sep="\t")
    require_columns(df, CALLER_REQUIRED, "Caller")

    # Raw workflow values are fractions (0-1), not percentage points.
    df["broad_detected"] = (
        df["abundance_metagenome"].fillna(0) > 0.0001
    )

    df["strict_caller"] = (
        (df["abundance_metagenome"].fillna(0) > 0.0001)
        & (df["callablePct"].fillna(0) > 0.10)
        & (df["confirmedPct"].fillna(0) >= 0.95)
        & (df["snpPct"].fillna(float("inf")) < 0.01)
        & (df["gapCount"].fillna(float("inf")) == 0)
        & (df["gapLength"].fillna(float("inf")) == 0)
    )

    outfile = outdir / "caller_classified.tsv"
    df.to_csv(outfile, sep="\t", index=False)

    print(f"Caller rows: {len(df)}")
    print(f"Broadly detected: {df['broad_detected'].sum()}")
    print(f"Strict Caller detection: {df['strict_caller'].sum()}")
    print(f"Written: {outfile}")


def classify_compare(path, outdir):
    df = pd.read_csv(path, sep="\t")
    require_columns(df, COMPARE_REQUIRED, "Compare")

    # Raw workflow values are fractions (0-1), not percentage points.
    df["donor_concordant"] = (
        (df["commonPct"].fillna(0) >= 0.10)
        & (df["singleAgreePct"].fillna(0) >= 0.9995)
    )

    outfile = outdir / "compare_classified.tsv"
    df.to_csv(outfile, sep="\t", index=False)

    print(f"Compare rows: {len(df)}")
    print(f"Donor-concordant rows: {df['donor_concordant'].sum()}")
    print(f"Written: {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Apply the STRACOL/REFIT2 StrainGE Caller and Compare "
            "thresholds to aggregated workflow outputs."
        )
    )

    parser.add_argument(
        "--caller",
        help="Aggregated StrainGE Caller TSV."
    )

    parser.add_argument(
        "--compare",
        help="Aggregated StrainGE Compare TSV."
    )

    parser.add_argument(
        "--outdir",
        default="classified",
        help="Output directory (default: classified)."
    )

    args = parser.parse_args()

    if not args.caller and not args.compare:
        parser.error("Provide --caller, --compare, or both.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.caller:
        classify_caller(args.caller, outdir)

    if args.compare:
        classify_compare(args.compare, outdir)


if __name__ == "__main__":
    main()
