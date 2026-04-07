#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path


def run(cmd, *, cwd=None, stdout=None):
    subprocess.run(cmd, cwd=cwd, stdout=stdout, check=True)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Run the StringTie SuperReads helper and produce a sorted BAM."
    )
    ap.add_argument("--read1", required=True)
    ap.add_argument("--read2", required=True)
    ap.add_argument("--install-dir", required=True, help="SuperReads install prefix.")
    ap.add_argument("--hisat-index", required=True, help="HISAT2 genome index prefix.")
    ap.add_argument("--gmap-dir", required=True, help="GMAP index directory.")
    ap.add_argument("--gmap-index", required=True, help="GMAP index name.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bam-out", required=True)
    ap.add_argument("--threads", required=True, type=int)
    ap.add_argument("--frag-len", required=True, type=float)
    ap.add_argument("--frag-std", required=True, type=float)
    ap.add_argument(
        "--hisat-rna-strandness",
        default="",
        help="HISAT2 paired-end strandness (FR/RF) or blank for unstranded.",
    )
    args = ap.parse_args()

    install_dir = Path(args.install_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    bam_out = Path(args.bam_out).resolve()

    create_superreads = install_dir / "bin" / "createSuperReads_RNA"
    assign_reads = install_dir / "bin" / "assign_reads"
    if not create_superreads.exists():
        raise SystemExit(f"Missing createSuperReads_RNA binary: {create_superreads}")
    if not assign_reads.exists():
        raise SystemExit(f"Missing assign_reads binary: {assign_reads}")

    conf_path = outdir / "super_reads.conf"
    conf_path.write_text(
        "\n".join(
            [
                "DATA",
                f"PE = pe {int(args.frag_len)} {int(args.frag_std)} {Path(args.read1).resolve()} {Path(args.read2).resolve()}",
                "END",
                "PARAMETERS",
                "STOP_AFTER_SUPERREADS=1",
                f"NUM_THREADS={args.threads}",
                "GRAPH_KMER_SIZE=auto",
                "END",
                "",
            ]
        )
    )

    work1_dir = outdir / "work1"
    superread_fasta = work1_dir / "superReadSequences.fasta"
    if not superread_fasta.exists():
        run([str(create_superreads), conf_path.name], cwd=outdir)
        run(["./assemble.sh"], cwd=outdir)
    if not superread_fasta.exists():
        raise SystemExit(f"Expected SuperReads FASTA not found: {superread_fasta}")

    sr_sam = outdir / "sr.sam"
    with sr_sam.open("w") as handle:
        run(
            [
                "gmap",
                "-f",
                "samse",
                "-t",
                str(args.threads),
                "-D",
                str(Path(args.gmap_dir).resolve()),
                "-d",
                args.gmap_index,
                str(superread_fasta),
            ],
            stdout=handle,
        )

    sr_quant_sam = outdir / "sr_quant.sam"
    run(
        [
            str(assign_reads),
            "-p",
            str(args.threads),
            "-w",
            str(work1_dir),
            "-s",
            str(sr_sam),
            "-o",
            str(sr_quant_sam),
        ]
    )

    short_sam = outdir / "short.sam"
    hisat_cmd = [
        "hisat2",
        "--dta",
        "-p",
        str(args.threads),
        "-x",
        args.hisat_index,
        "-1",
        args.read1,
        "-2",
        args.read2,
        "-S",
        str(short_sam),
    ]
    if args.hisat_rna_strandness:
        hisat_cmd.extend(["--rna-strandness", args.hisat_rna_strandness])
    run(hisat_cmd)

    sr_quant_bam = outdir / "sr_quant.bam"
    sr_quant_sorted_bam = outdir / "sr_quant.sorted.bam"
    short_bam = outdir / "short.bam"
    short_sorted_bam = outdir / "short.sorted.bam"
    merged_bam = outdir / "sr_merge.unsorted.bam"

    run(["samtools", "view", "-@", str(args.threads), "-b", "-o", str(sr_quant_bam), str(sr_quant_sam)])
    run(["samtools", "sort", "-@", str(args.threads), "-o", str(sr_quant_sorted_bam), str(sr_quant_bam)])
    run(["samtools", "view", "-@", str(args.threads), "-b", "-o", str(short_bam), str(short_sam)])
    run(["samtools", "sort", "-@", str(args.threads), "-o", str(short_sorted_bam), str(short_bam)])
    run(
        [
            "samtools",
            "merge",
            "-f",
            "-@",
            str(args.threads),
            str(merged_bam),
            str(sr_quant_sorted_bam),
            str(short_sorted_bam),
        ]
    )
    run(["samtools", "sort", "-@", str(args.threads), "-o", str(bam_out), str(merged_bam)])
    run(["samtools", "index", "-@", str(args.threads), str(bam_out), f"{bam_out}.bai"])

    for path in [
        sr_sam,
        sr_quant_sam,
        short_sam,
        sr_quant_bam,
        sr_quant_sorted_bam,
        short_bam,
        short_sorted_bam,
        merged_bam,
    ]:
        if path.exists():
            path.unlink()


if __name__ == "__main__":
    main()
