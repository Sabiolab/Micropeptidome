#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {"locus", "cds_chr", "cds_starts", "cds_ends", "cds_strand"}


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build a SAF annotation for featureCounts from a locus summary CSV."
    )
    ap.add_argument("--loci_csv", required=True, help="Input CSV with locus positions.")
    ap.add_argument("--out_saf", required=True, help="Output SAF file.")
    args = ap.parse_args()

    df = pd.read_csv(args.loci_csv, usecols=list(REQUIRED_COLUMNS), low_memory=False)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {args.loci_csv}: {sorted(missing)}")

    starts = pd.to_numeric(df["cds_starts"], errors="coerce")
    ends = pd.to_numeric(df["cds_ends"], errors="coerce")
    if starts.isna().any() or ends.isna().any():
        raise SystemExit(f"Non-numeric cds_starts/cds_ends found in {args.loci_csv}")

    saf = pd.DataFrame(
        {
            "GeneID": df["locus"].astype(str),
            "Chr": df["cds_chr"].astype(str),
            "Start": pd.concat([starts, ends], axis=1).min(axis=1).astype(int),
            "End": pd.concat([starts, ends], axis=1).max(axis=1).astype(int),
            "Strand": df["cds_strand"].astype(str),
        }
    )

    saf = saf.drop_duplicates(subset=["GeneID"]).copy()
    saf.to_csv(Path(args.out_saf), sep="\t", index=False)


if __name__ == "__main__":
    main()
