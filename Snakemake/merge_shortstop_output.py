#!/usr/bin/env python3
# I run it with this command, but change the directories as needed:
# python merge_shortstop_output.py \
#   --root /storage/scratch01/users/sbarber/Workdir/results_shortstop \
#   --outdir /storage/scratch01/users/sbarber/Workdir/merged_per_sample

import argparse
import os
import re
from pathlib import Path
import pandas as pd


def clean_orf_id(x: str) -> str:
    # Normalize orf_id like: '"""cds.STRG..."""' -> 'cds.STRG...'
    if pd.isna(x):
        return x
    s = str(x).strip()
    # remove surrounding quotes repeatedly
    # e.g. """cds.X""" -> cds.X
    s = re.sub(r'^\"+', '', s)
    s = re.sub(r'\"+$', '', s)
    # also remove single quotes if present
    s = s.strip('"').strip("'")
    return s


def read_table(path: Path) -> pd.DataFrame:
    """
    Read CSV/TSV with delimiter auto-detection. It is ususally comma, but just in case
    """
    return pd.read_csv(path, sep=None, engine="python")

def merge_one_sample(sample_dir: Path, out_dir: Path, min_prob: float | None) -> Path:
    """
    sample_dir should contain:
      sample_dir/shortstop/shortstop_output/predictions/sams.csv
      sample_dir/shortstop/shortstop_output/sequences/unknown_sequences.csv
    """
    pred_path = (
        sample_dir
        / "shortstop"
        / "shortstop_output"
        / "predictions"
        / "sams.csv"
    )
    seq_path = (
        sample_dir
        / "shortstop"
        / "shortstop_output"
        / "sequences"
        / "unknown_sequences.csv"
    )

    if not pred_path.exists():
        raise FileNotFoundError(f"Missing predictions file: {pred_path}")
    if not seq_path.exists():
        raise FileNotFoundError(f"Missing sequences file: {seq_path}")

    pred = read_table(pred_path)
    seq  = read_table(seq_path)

    if "orf_id" not in pred.columns or "orf_id" not in seq.columns:
        raise ValueError(f"orf_id missing in one of the inputs for sample {sample_dir.name}")

    pred["orf_id_clean"] = pred["orf_id"].map(clean_orf_id)
    seq["orf_id_clean"]  = seq["orf_id"].map(clean_orf_id)

    if min_prob is not None and "probability" in pred.columns:
        pred = pred[pred["probability"] >= min_prob].copy()

    merged = pred.merge(
        seq,
        on="orf_id_clean",
        how="inner",
        suffixes=("_pred", "_seq")
    )

    merged.drop( # Remove redundant orf_id columns
        columns=[c for c in ["orf_id_pred", "orf_id_seq"] if c in merged.columns],
        inplace=True
    )

    merged.insert(0, "orf_id", merged["orf_id_clean"])
    merged.drop(columns=["orf_id_clean"], inplace=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{sample_dir.name}.merged.csv"
    merged.to_csv(out_path, index=False)
    return out_path


def main():
    ap = argparse.ArgumentParser(
        description="Merge ShortStop sam_secreted predictions with unknown_sequences by orf_id for each sample."
    )
    ap.add_argument("--root", required=True,
                    help="Root directory containing sample folders (e.g., results_shortstop/)")
    ap.add_argument("--outdir", required=True,
                    help="Where to write per-sample merged CSVs (e.g., merged_per_sample/)")
    ap.add_argument("--min_prob", type=float, default=None,
                    help="Optional: keep only predictions with probability >= min_prob")
    ap.add_argument("--samples", nargs="*", default=None,
                    help="Optional: specific sample folder names to process (default: auto-discover)")
    args = ap.parse_args()

    root = Path(args.root)
    outdir = Path(args.outdir)

    if args.samples:
        samples = [root / s for s in args.samples]
    else:
        # auto-discover: folders directly under root
        samples = [p for p in root.iterdir() if p.is_dir()]

    ok = 0
    for sdir in sorted(samples):
        try:
            out_path = merge_one_sample(sdir, outdir, args.min_prob)
            print(f"[OK] {sdir.name} -> {out_path}")
            ok += 1
        except Exception as e:
            print(f"[FAIL] {sdir.name}: {e}")

    print(f"Done. Successful: {ok}/{len(samples)}")


if __name__ == "__main__":
    main()
