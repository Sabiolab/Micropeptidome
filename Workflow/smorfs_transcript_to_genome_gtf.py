#!/usr/bin/env python3

import argparse # to handle command-line arguments
from collections import defaultdict, namedtuple

Exon = namedtuple("Exon", ["chrom", "start", "end", "strand"]) # object to store exon information

def parse_transcript_exons(gtf_path): # parse GTF to get exons foe each transcript
    exons_by_tx = defaultdict(list) # 

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"): # skip comments
                continue
            cols = line.rstrip("\n").split("\t") # split line into columns
            if len(cols) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = cols # unpack columns
            if feature != "exon": # only process exon features
                continue

            transcript_id = None # initialize transcript_id
            for field in attrs.split(";"): # parse attributes to find transcript_id
                field = field.strip()
                if field.startswith("transcript_id"): # i
                    parts = field.split()
                    if len(parts) > 1:
                        transcript_id = parts[1].strip('"')
                    break

            if transcript_id is None: # if no transcript_id found, skip
                continue

            exons_by_tx[transcript_id].append( # add exon to the list for this transcript
                Exon(chrom=chrom, start=int(start), end=int(end), strand=strand)
            )

    tx_to_exons = {} # final mapping of transcript_id to (strand, sorted exons)
    for tx_id, exons in exons_by_tx.items():
        if not exons:
            continue
        strand = exons[0].strand # assume all exons have the same strand
        exons_sorted = sorted(exons, key=lambda e: e.start) 
        tx_to_exons[tx_id] = (strand, exons_sorted) # store strand and sorted exons

    return tx_to_exons


def map_orf_to_genome(tx_exons, strand, orf_start_tx, orf_end_tx): # map ORF coordinates from transcrip to genome
    segments = []

    if strand == "+":
        ordered_exons = tx_exons
    else:
        ordered_exons = list(reversed(tx_exons))

    tx_pos = 1
    for exon in ordered_exons:
        exon_len = exon.end - exon.start + 1
        exon_tx_start = tx_pos
        exon_tx_end = tx_pos + exon_len - 1

        ov_start = max(orf_start_tx, exon_tx_start)
        ov_end = min(orf_end_tx, exon_tx_end)

        if ov_start <= ov_end:
            offset_start = ov_start - exon_tx_start
            offset_end = ov_end - exon_tx_start

            if strand == "+":
                g_start = exon.start + offset_start
                g_end = exon.start + offset_end
            else:
                g_end = exon.end - offset_start
                g_start = exon.end - offset_end
                if g_start > g_end:
                    g_start, g_end = g_end, g_start

            segments.append((exon.chrom, g_start, g_end, exon.strand))

        tx_pos += exon_len
        if tx_pos > orf_end_tx:
            break

    return segments


def build_smorfs_genomic_gtf(merged_gtf, smorfs_gff3, out_gtf):
    tx_to_exons = parse_transcript_exons(merged_gtf)

    # Store ORF → CDS segments first
    orf_segments = defaultdict(list)
    orf_meta = {}  # (chrom, strand, tx_id)

    with open(smorfs_gff3) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            seqid, source, feature, start, end, score, strand, frame, attrs = cols
            if feature != "CDS":
                continue

            tx_id = seqid

            parent = None
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("Parent="):
                    parent = field.split("=", 1)[1]
                    break

            if tx_id not in tx_to_exons and parent:
                base = parent.split(".p")[0]
                if base in tx_to_exons:
                    tx_id = base

            if tx_id not in tx_to_exons:
                continue

            orf_start_tx = int(start)
            orf_end_tx = int(end)
            strand_tx, exons = tx_to_exons[tx_id]

            segments = map_orf_to_genome(exons, strand_tx, orf_start_tx, orf_end_tx)
            if not segments:
                continue

            orf_id = None
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("ID="):
                    orf_id = field.split("=", 1)[1]
                    break
            if orf_id is None:
                orf_id = f"{tx_id}_orf_{orf_start_tx}_{orf_end_tx}"

            for seg in segments:
                orf_segments[orf_id].append(seg)

            # save meta for transcript line
            orf_meta[orf_id] = (strand_tx, tx_id)

    # Write transcript + CDS
    with open(out_gtf, "w") as fout:
        for orf_id, segments in orf_segments.items():
            strand, tx_id = orf_meta[orf_id]

            # transcript span
            tx_chrom = segments[0][0]
            tx_start = min(s[1] for s in segments)
            tx_end = max(s[2] for s in segments)

            # write transcript line
            tx_attrs = f'gene_id "{orf_id}"; transcript_id "{orf_id}";'
            fout.write("\t".join([
                tx_chrom,
                "smORFmapper",
                "transcript",
                str(tx_start),
                str(tx_end),
                ".",
                strand,
                ".",
                tx_attrs
            ]) + "\n")

            # write CDS segments
            for chrom, g_start, g_end, g_strand in segments:
                cds_attrs = f'gene_id "{orf_id}"; transcript_id "{orf_id}";'
                fout.write("\t".join([
                    chrom,
                    "smORFmapper",
                    "CDS",
                    str(g_start),
                    str(g_end),
                    ".",
                    g_strand,
                    "0",
                    cds_attrs
                ]) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Map smORFs from transcript coords (TransDecoder GFF3) to genomic GTF.")
    ap.add_argument("--merged_gtf", required=True)
    ap.add_argument("--smorfs_gff3", required=True)
    ap.add_argument("--out_gtf", required=True)
    args = ap.parse_args()

    build_smorfs_genomic_gtf(args.merged_gtf, args.smorfs_gff3, args.out_gtf)


if __name__ == "__main__":
    main()
