#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Combine strict StrainGE Caller detection with exact-donor "
            "StrainGE Compare concordance for one recipient timepoint."
        )
    )

    parser.add_argument(
        "--caller",
        required=True,
        help="Caller table produced by apply_thresholds.py."
    )

    parser.add_argument(
        "--compare",
        required=True,
        help=(
            "Compare table produced by apply_thresholds.py for the "
            "corresponding exact donor-recipient pair."
        )
    )

    parser.add_argument(
        "--recipient-sample",
        required=True,
        help="Recipient sample name as stored in the Caller samplename column."
    )

    parser.add_argument(
        "--timepoint",
        required=True,
        choices=["T00", "T03", "T12"],
        help="Recipient sampling timepoint."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV."
    )

    args = parser.parse_args()

    caller = pd.read_csv(args.caller, sep="\t")
    compare = pd.read_csv(args.compare, sep="\t")

    caller_required = {"samplename", "strain", "strict_caller"}
    compare_required = {"ref", "donor_concordant"}

    missing_caller = caller_required - set(caller.columns)
    missing_compare = compare_required - set(compare.columns)

    if missing_caller:
        raise SystemExit(
            "ERROR: Caller table missing columns: "
            + ", ".join(sorted(missing_caller))
        )

    if missing_compare:
        raise SystemExit(
            "ERROR: Compare table missing columns: "
            + ", ".join(sorted(missing_compare))
        )

    # Select the requested recipient sample.
    caller = caller[
        caller["samplename"].astype(str) == args.recipient_sample
    ].copy()

    if caller.empty:
        raise SystemExit(
            f"ERROR: Recipient sample not found in Caller table: "
            f"{args.recipient_sample}"
        )

    # One Caller row per reference strain and one Compare row per reference.
    merged = caller.merge(
        compare[["ref", "donor_concordant"]],
        left_on="strain",
        right_on="ref",
        how="left"
    )

    # No Compare row means donor concordance was not demonstrated.
    merged["donor_concordant"] = (
        merged["donor_concordant"]
        .fillna(False)
        .astype(bool)
    )

    merged["confirmed_donor_strain"] = (
        merged["strict_caller"].astype(bool)
        & merged["donor_concordant"]
    )

    merged.insert(0, "recipient_sample", args.recipient_sample)
    merged.insert(1, "recipient_timepoint", args.timepoint)

    if args.timepoint == "T03":
        merged["confirmed_T03_engraftment"] = (
            merged["confirmed_donor_strain"]
        )

    outfile = Path(args.output)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    merged.to_csv(outfile, sep="\t", index=False)

    print(f"Recipient sample: {args.recipient_sample}")
    print(f"Timepoint: {args.timepoint}")
    print(f"Reference strains evaluated: {len(merged)}")
    print(
        "Confirmed donor-concordant strains: "
        f"{merged['confirmed_donor_strain'].sum()}"
    )
    print(f"Written: {outfile}")


if __name__ == "__main__":
    main()
