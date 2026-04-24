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


def parse_gtf_attrs(attr_str: str) -> dict[str, str]:
    attrs = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, val = part.split(" ", 1)
        attrs[key] = val.strip().strip('"')
    return attrs


def load_exact_locations(smorf_gtf: Path) -> dict[str, str]:
    coords_by_orf = {}
    if not smorf_gtf.exists():
        raise FileNotFoundError(f"Missing smORF GTF: {smorf_gtf}")

    cds_rows = []
    with open(smorf_gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = parse_gtf_attrs(parts[8])
            orf_id = attrs.get("gene_id") or attrs.get("transcript_id")
            if not orf_id:
                continue
            cds_rows.append(
                {
                    "orf_id": orf_id,
                    "chrom": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
            )

    if not cds_rows:
        return coords_by_orf

    cds_df = pd.DataFrame(cds_rows)
    for orf_id, sub in cds_df.groupby("orf_id", sort=False):
        strand = str(sub["strand"].iloc[0])
        chrom = str(sub["chrom"].iloc[0])
        ascending = strand != "-"
        sub = sub.sort_values(["start", "end"], ascending=ascending)
        blocks = ",".join(f"{int(row.start)}-{int(row.end)}" for row in sub.itertuples())
        coords_by_orf[str(orf_id)] = f"{chrom}:{blocks}:{strand}"

    return coords_by_orf


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


def summarize_loci(df: pd.DataFrame, exact_locations: dict[str, str] | None = None) -> pd.DataFrame:
    df = normalize_locus_df(df)
    if exact_locations is not None and "orf_id" in df.columns:
        df = df.copy()
        df["exact_genomic_location"] = df["orf_id"].astype(str).map(exact_locations)

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
        "exact_genomic_location",
    ]:
        if extra in df.columns:
            keep_cols.append(extra)
    sub = df[keep_cols].copy()

    agg_spec = dict(
        cds_chr=("cds_chr", "first"),
        cds_starts=("cds_starts", "first"),
        cds_ends=("cds_ends", "first"),
        cds_strand=("cds_strand", "first"),
        n_orfs=("locus", "size"),
    )

    if "sam_probability" in sub.columns:
        agg_spec["max_prob"] = ("sam_probability", "max")

    if "orf_id" in sub.columns:
        agg_spec["orf_ids"] = ("orf_id", unique_nonnull)

    if "cds_seq" in sub.columns:
        agg_spec["cds_seq"] = ("cds_seq", first_nonnull)

    if "aa_seq" in sub.columns:
        agg_spec["aa_seq"] = ("aa_seq", first_nonnull)

    if "smorf_type" in sub.columns:
        agg_spec["smorf_type"] = ("smorf_type", most_common_nonnull)
        agg_spec["smorf_types"] = ("smorf_type", unique_nonnull)

    if "exact_genomic_location" in sub.columns:
        agg_spec["exact_genomic_location"] = ("exact_genomic_location", unique_nonnull)

    agg = sub.groupby("locus", as_index=False).agg(**agg_spec)

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
        ["n_orfs", "cds_chr", "cds_starts", "cds_ends"],
        ascending=[False, True, True, True],
    ).to_csv(full_out, index=False)

    shared = agg[agg["n_orfs"] >= min_patients].copy()
    shared_out = Path(f"{out_prefix}.shared_ge{min_patients}.csv")
    shared.sort_values(
        ["n_orfs", "cds_chr", "cds_starts", "cds_ends"],
        ascending=[False, True, True, True],
    ).to_csv(shared_out, index=False)

    print(f"Total loci: {len(agg)}")
    print(f"[OK] Wrote: {full_out}")
    print(f"[OK] Wrote: {shared_out}")


def aggregate_condition_mode(
    merged_csv: Path, out_csv: Path, smorf_gtf: Path | None = None
) -> None:
    if not merged_csv.exists():
        raise FileNotFoundError(f"Missing merged ShortStop table: {merged_csv}")

    df = pd.read_csv(merged_csv)
    exact_locations = load_exact_locations(smorf_gtf) if smorf_gtf else None
    agg = summarize_loci(df, exact_locations=exact_locations)
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
        "--smorf_gtf",
        help="Optional condition-level smORF GTF used to add exact genomic CDS block coordinates.",
    )
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
        smorf_gtf=Path(args.smorf_gtf) if args.smorf_gtf else None,
    )


if __name__ == "__main__":
    main()
