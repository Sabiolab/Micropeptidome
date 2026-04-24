#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def load_isoform_tpms(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    id_col = "transcript_id" if "transcript_id" in df.columns else df.columns[0]
    if "TPM" not in df.columns:
        raise SystemExit(f"No TPM column found in {path}")
    return df.set_index(id_col)["TPM"].astype(float)


def add_tpms(summary_csv: Path, rsem_dir: Path, out_csv: Path) -> None:
    summ = pd.read_csv(summary_csv)
    if "locus" not in summ.columns:
        raise SystemExit(f"{summary_csv} must contain a locus column")

    tpm_by_sample = {}
    for sample_dir in sorted(rsem_dir.glob("*")):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        iso = sample_dir / f"{sample}.isoforms.results"
        if iso.exists():
            tpm_by_sample[sample] = load_isoform_tpms(iso)

    sample_list = sorted(tpm_by_sample)

    def row_tpms(row) -> str:
        locus = str(row["locus"])
        values = []
        for sample in sample_list:
            series = tpm_by_sample.get(sample)
            value = 0.0
            if series is not None and locus in series.index:
                value = float(series.loc[locus])
            values.append(f"{value:.6f}")
        return ",".join(values)

    summ["samples"] = ",".join(sample_list)
    summ["tpms"] = summ.apply(row_tpms, axis=1)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    summ.to_csv(out_csv, index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary_csv")
    ap.add_argument("--out_csv")
    ap.add_argument("--all_loci_csv")
    ap.add_argument("--shared_csv")
    ap.add_argument("--out_all_loci_csv")
    ap.add_argument("--out_shared_csv")
    ap.add_argument("--rsem_dir", required=True)
    args = ap.parse_args()

    rsem_dir = Path(args.rsem_dir)

    if args.summary_csv:
        if not args.out_csv:
            raise SystemExit("--out_csv is required with --summary_csv")
        add_tpms(Path(args.summary_csv), rsem_dir, Path(args.out_csv))
        return

    required = [
        args.all_loci_csv,
        args.shared_csv,
        args.out_all_loci_csv,
        args.out_shared_csv,
    ]
    if any(value is None for value in required):
        raise SystemExit(
            "Provide either --summary_csv/--out_csv or the legacy "
            "--all_loci_csv/--shared_csv/--out_all_loci_csv/--out_shared_csv arguments."
        )

    add_tpms(Path(args.all_loci_csv), rsem_dir, Path(args.out_all_loci_csv))
    add_tpms(Path(args.shared_csv), rsem_dir, Path(args.out_shared_csv))


if __name__ == "__main__":
    main()
