#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


LOCUS_COLS = ["cds_chr", "cds_starts", "cds_ends", "cds_strand"]


def first_nonnull(series: pd.Series):
    s = series.dropna()
    if len(s) == 0:
        return pd.NA
    return s.iloc[0]


def unique_nonnull(series: pd.Series):
    s = [str(x) for x in series.dropna() if str(x) != ""]
    if not s:
        return pd.NA
    return ",".join(sorted(set(s)))


def most_common_nonnull(series: pd.Series):
    s = [str(x) for x in series.dropna() if str(x) != ""]
    if not s:
        return pd.NA
    counts = {}
    for val in s:
        counts[val] = counts.get(val, 0) + 1
    max_count = max(counts.values())
    top = sorted([k for k, v in counts.items() if v == max_count])
    return ",".join(top)


def normalize_locus_df(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in LOCUS_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required locus columns: {missing}")

    df = df.copy()
    df["cds_starts"] = pd.to_numeric(df["cds_starts"], errors="coerce")
    df["cds_ends"] = pd.to_numeric(df["cds_ends"], errors="coerce")
    df["locus"] = (
        df["cds_chr"].astype(str)
        + ":"
        + df["cds_starts"].astype("Int64").astype(str)
        + "-"
        + df["cds_ends"].astype("Int64").astype(str)
        + ":"
        + df["cds_strand"].astype(str)
    )
    return df


def summarize_loci(df: pd.DataFrame, patients: list[str] | None = None) -> pd.DataFrame:
    df = normalize_locus_df(df)

    keep_cols = ["locus"] + LOCUS_COLS
    for extra in [
        "orf_id",
        "sam_probability",
        "classification",
        "aa_seq",
        "length",
        "type",
        "cds_seq",
        "smorf_type",
    ]:
        if extra in df.columns:
            keep_cols.append(extra)
    if "sample" in df.columns:
        keep_cols.append("sample")

    sub = df[keep_cols].copy()

    agg_spec = dict(
        cds_chr=("cds_chr", "first"),
        cds_starts=("cds_starts", "first"),
        cds_ends=("cds_ends", "first"),
        cds_strand=("cds_strand", "first"),
        n_rows=("locus", "size"),
    )

    if "sam_probability" in sub.columns:
        agg_spec["max_prob"] = ("sam_probability", "max")

    if "cds_seq" in sub.columns:
        agg_spec["cds_seq"] = ("cds_seq", first_nonnull)

    if "aa_seq" in sub.columns:
        agg_spec["aa_seq"] = ("aa_seq", first_nonnull)

    if "smorf_type" in sub.columns:
        agg_spec["smorf_type"] = ("smorf_type", most_common_nonnull)
        agg_spec["smorf_types"] = ("smorf_type", unique_nonnull)

    if patients is None:
        if "sample" not in sub.columns:
            raise ValueError("Sample-level aggregation requires a 'sample' column.")
        agg_spec["n_patients"] = ("sample", "nunique")
        agg_spec["patients"] = ("sample", lambda x: ",".join(sorted(set(x))))
    else:
        patient_list = [str(x).strip() for x in patients if str(x).strip()]
        patient_list = list(dict.fromkeys(sorted(patient_list)))
        if not patient_list:
            raise ValueError("At least one patient must be provided for merged-condition mode.")

    agg = sub.groupby("locus", as_index=False).agg(**agg_spec)

    if patients is not None:
        patient_str = ",".join(patient_list)
        agg["n_patients"] = len(patient_list)
        agg["patients"] = patient_str

    return agg.sort_values(
        ["cds_chr", "cds_starts", "cds_ends"],
        ascending=[True, True, True],
    )


def aggregate_directory_mode(
    merged_dir: Path,
    out_prefix: str,
    min_patients: int,
    samples: list[str] | None,
) -> None:
    if samples:
        files = [merged_dir / f"{sample}.merged.csv" for sample in samples]
    else:
        files = sorted(merged_dir.glob("*.merged.csv"))
    if not files:
        raise SystemExit(f"No *.merged.csv files found in: {merged_dir}")

    rows = []
    for file_path in files:
        if not file_path.exists():
            raise FileNotFoundError(f"Missing merged ShortStop table: {file_path}")
        sample = file_path.name.replace(".merged.csv", "")
        df = pd.read_csv(file_path)
        df["sample"] = sample
        rows.append(df)

    agg = summarize_loci(pd.concat(rows, ignore_index=True))

    full_out = Path(f"{out_prefix}.all_loci.csv")
    agg.sort_values(
        ["n_patients", "cds_chr", "cds_starts", "cds_ends"],
        ascending=[False, True, True, True],
    ).to_csv(full_out, index=False)

    shared = agg[agg["n_patients"] >= min_patients].copy()
    shared_out = Path(f"{out_prefix}.shared_ge{min_patients}.csv")
    shared.sort_values(
        ["n_patients", "cds_chr", "cds_starts", "cds_ends"],
        ascending=[False, True, True, True],
    ).to_csv(shared_out, index=False)

    print(f"Total loci: {len(agg)}")
    print(f"[OK] Wrote: {full_out}")
    print(f"[OK] Wrote: {shared_out}")


def aggregate_condition_mode(merged_csv: Path, out_csv: Path, patients: list[str]) -> None:
    if not merged_csv.exists():
        raise FileNotFoundError(f"Missing merged ShortStop table: {merged_csv}")

    df = pd.read_csv(merged_csv)
    agg = summarize_loci(df, patients=patients)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(out_csv, index=False)
    print(f"[OK] Wrote: {out_csv}")


def main():
    ap = argparse.ArgumentParser(
        description="Aggregate merged ShortStop tables by genomic locus."
    )
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--merged_dir", help="Directory containing per-sample merged CSVs.")
    source.add_argument("--merged_csv", help="Single merged CSV for one condition.")

    ap.add_argument("--out_prefix", help="Prefix for directory mode outputs.")
    ap.add_argument("--out_csv", help="Output CSV for single-condition mode.")
    ap.add_argument(
        "--min_patients",
        type=int,
        default=2,
        help="Keep loci observed in at least this many patients in directory mode.",
    )
    ap.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional subset of sample names for directory mode.",
    )
    ap.add_argument(
        "--patients",
        nargs="*",
        default=None,
        help="Patient IDs to stamp onto a condition-level summary.",
    )
    args = ap.parse_args()

    if args.merged_dir:
        if not args.out_prefix:
            raise SystemExit("--out_prefix is required with --merged_dir")
        aggregate_directory_mode(
            merged_dir=Path(args.merged_dir),
            out_prefix=args.out_prefix,
            min_patients=args.min_patients,
            samples=args.samples,
        )
        return

    if not args.out_csv:
        raise SystemExit("--out_csv is required with --merged_csv")
    aggregate_condition_mode(
        merged_csv=Path(args.merged_csv),
        out_csv=Path(args.out_csv),
        patients=args.patients or [],
    )


if __name__ == "__main__":
    main()
