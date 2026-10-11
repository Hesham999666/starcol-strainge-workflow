#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Determine T12 persistence among strains with confirmed "
            "T03 engraftment."
        )
    )

    parser.add_argument(
        "--t03",
        required=True,
        help="T03 output from combine_pair.py."
    )

    parser.add_argument(
        "--t12",
        required=True,
        help="T12 output from combine_pair.py."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV."
    )

    args = parser.parse_args()

    t03 = pd.read_csv(args.t03, sep="\t")
    t12 = pd.read_csv(args.t12, sep="\t")

    required = {"strain", "confirmed_donor_strain"}

    missing_t03 = required - set(t03.columns)
    missing_t12 = required - set(t12.columns)

    if missing_t03:
        raise SystemExit(
            "ERROR: T03 table missing columns: "
            + ", ".join(sorted(missing_t03))
        )

    if missing_t12:
        raise SystemExit(
            "ERROR: T12 table missing columns: "
            + ", ".join(sorted(missing_t12))
        )

    t03_small = t03[
        ["strain", "confirmed_donor_strain"]
    ].rename(
        columns={
            "confirmed_donor_strain":
            "confirmed_T03_engraftment"
        }
    )

    t12_small = t12[
        ["strain", "confirmed_donor_strain"]
    ].rename(
        columns={
            "confirmed_donor_strain":
            "confirmed_T12_donor_strain"
        }
    )

    merged = t03_small.merge(
        t12_small,
        on="strain",
        how="outer"
    )

    merged["confirmed_T03_engraftment"] = (
        merged["confirmed_T03_engraftment"]
        .fillna(False)
        .astype(bool)
    )

    merged["confirmed_T12_donor_strain"] = (
        merged["confirmed_T12_donor_strain"]
        .fillna(False)
        .astype(bool)
    )

    merged["persistent_at_T12"] = (
        merged["confirmed_T03_engraftment"]
        & merged["confirmed_T12_donor_strain"]
    )

    outfile = Path(args.output)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(outfile, sep="\t", index=False)

    eligible = merged["confirmed_T03_engraftment"].sum()
    persistent = merged["persistent_at_T12"].sum()

    print(f"Confirmed T03 strains evaluated: {eligible}")
    print(f"Persistent at T12: {persistent}")
    print(f"Written: {outfile}")


if __name__ == "__main__":
    main()
