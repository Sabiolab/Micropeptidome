#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def load_isoform_expected_counts(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        raise SystemExit(f"{path} is empty")

    id_col = "transcript_id" if "transcript_id" in df.columns else df.columns[0]

    count_col = None
    for candidate in ("expected_count", "expected_counts"):
        if candidate in df.columns:
            count_col = candidate
            break

    if count_col is None:
        raise SystemExit(
            f"No expected-count column found in {path}. "
            f"Available columns: {', '.join(df.columns)}"
        )

    return df.set_index(id_col)[count_col].astype(float)


def unique_loci_in_order(summary_csv: Path) -> list[str]:
    df = pd.read_csv(summary_csv)
    if "locus" not in df.columns:
        raise SystemExit(f"{summary_csv} must contain a locus column")

    loci = [str(x) for x in df["locus"].tolist() if pd.notna(x)]
    return list(dict.fromkeys(loci))


def export_counts_matrix(summary_csv: Path, rsem_dir: Path, out_tsv: Path) -> None:
    loci = unique_loci_in_order(summary_csv)
    if not loci:
        raise SystemExit(f"No loci found in {summary_csv}")

    counts_by_sample = {}
    for sample_dir in sorted(rsem_dir.glob("*")):
        if not sample_dir.is_dir():
            continue

        sample = sample_dir.name
        isoform_results = sample_dir / f"{sample}.isoforms.results"
        if not isoform_results.exists():
            continue

        counts = load_isoform_expected_counts(isoform_results)
        counts.index = counts.index.astype(str)
        counts_by_sample[sample] = counts.reindex(loci, fill_value=0.0)

    if not counts_by_sample:
        raise SystemExit(f"No *.isoforms.results files found under {rsem_dir}")

    matrix = pd.DataFrame(counts_by_sample, index=loci)
    matrix.index.name = "locus"
    matrix.to_csv(out_tsv, sep="\t", float_format="%.6f")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Export an RSEM expected-count matrix for loci present in a summary CSV."
    )
    ap.add_argument("--loci_csv", required=True, help="CSV containing a locus column.")
    ap.add_argument("--rsem_dir", required=True, help="Directory containing per-sample RSEM outputs.")
    ap.add_argument("--out_tsv", required=True, help="Output TSV count matrix.")
    args = ap.parse_args()

    export_counts_matrix(Path(args.loci_csv), Path(args.rsem_dir), Path(args.out_tsv))


if __name__ == "__main__":
    main()
