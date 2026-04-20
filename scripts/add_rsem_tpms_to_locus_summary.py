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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all_loci_csv", required=True)
    ap.add_argument("--rsem_dir", required=True)
    ap.add_argument("--tpm_threshold", type=float, default=0.0,
                    help="Min TPM to count a sample as expressing a smORF (default: 0)")
    ap.add_argument("--min_patients", type=int, default=2,
                    help="Min samples expressing a locus to include in shared output (default: 2)")
    ap.add_argument("--out_all_loci_csv", required=True)
    ap.add_argument("--out_shared_csv", required=True)
    ap.add_argument("--out_shared_tpm_csv", required=True)
    args = ap.parse_args()

    rsem_dir = Path(args.rsem_dir)
    summ = pd.read_csv(args.all_loci_csv)
    if "locus" not in summ.columns:
        raise SystemExit(f"{args.all_loci_csv} must contain a 'locus' column")

    # Load all per-sample TPM vectors
    tpm_by_sample = {}
    for sample_dir in sorted(rsem_dir.glob("*")):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        iso = sample_dir / f"{sample}.isoforms.results"
        if iso.exists():
            tpm_by_sample[sample] = load_isoform_tpms(iso)

    if not tpm_by_sample:
        raise SystemExit(f"No isoforms.results files found under {rsem_dir}")

    samples = sorted(tpm_by_sample.keys())

    # Build per-locus TPM matrix
    tpm_records = []
    for locus in summ["locus"]:
        row = {}
        for s in samples:
            v = tpm_by_sample[s].get(locus, 0.0)
            row[s] = float(v)
        tpm_records.append(row)

    tpm_df = pd.DataFrame(tpm_records, index=summ.index)

    # Count samples with TPM > threshold
    expressed = tpm_df > args.tpm_threshold
    summ["n_patients"] = expressed.sum(axis=1)
    summ["patients"] = expressed.apply(
        lambda row: ",".join(s for s in samples if row[s]), axis=1
    )
    summ["tpms"] = tpm_df.apply(
        lambda row: ",".join(f"{row[s]:.6f}" for s in samples), axis=1
    )

    # all_loci with TPMs
    summ.sort_values(
        ["n_patients", "cds_chr", "cds_starts", "cds_ends"],
        ascending=[False, True, True, True]
    ).to_csv(args.out_all_loci_csv, index=False)

    # shared: loci expressed in >= min_patients samples
    shared = summ[summ["n_patients"] >= args.min_patients].copy()
    shared.sort_values(
        ["n_patients", "cds_chr", "cds_starts", "cds_ends"],
        ascending=[False, True, True, True]
    )

    # shared without TPM columns
    shared.drop(columns=["tpms"]).to_csv(args.out_shared_csv, index=False)

    # shared with TPM columns
    shared.to_csv(args.out_shared_tpm_csv, index=False)

    print(f"Total loci: {len(summ)}")
    print(f"Shared loci (>= {args.min_patients} patients, TPM > {args.tpm_threshold}): {len(shared)}")
    for k in [2, 5, 10, 15]:
        n_k = (summ["n_patients"] >= k).sum()
        print(f"  >= {k} patients: {n_k}")
    print(f"[OK] Wrote: {args.out_all_loci_csv}")
    print(f"[OK] Wrote: {args.out_shared_csv}")
    print(f"[OK] Wrote: {args.out_shared_tpm_csv}")


if __name__ == "__main__":
    main()
