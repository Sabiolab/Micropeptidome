#!/usr/bin/env python3

import argparse
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


def parse_gff3_attrs(attr_text: str) -> dict[str, str]:
    attrs = {}
    for field in attr_text.strip().split(";"):
        field = field.strip()
        if not field or "=" not in field:
            continue
        key, value = field.split("=", 1)
        attrs[key] = value
    return attrs


def load_transcripts(fasta_path: str) -> dict[str, Seq]:
    transcripts = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        transcripts[record.id] = record.seq
    if not transcripts:
        raise ValueError(f"No transcript sequences found in {fasta_path}")
    return transcripts


def load_orfs(gff3_path: str) -> dict[str, dict]:
    orfs: dict[str, dict] = {}
    cds_by_orf: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)

    with open(gff3_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, _, feature, start, end, _, strand, phase, attrs_text = parts
            attrs = parse_gff3_attrs(attrs_text)

            if feature == "mRNA":
                orf_id = attrs.get("ID")
                if orf_id:
                    orfs.setdefault(orf_id, {})
                    orfs[orf_id]["transcript_id"] = seqid
                    orfs[orf_id]["strand"] = strand

            elif feature == "CDS":
                parent = attrs.get("Parent")
                if not parent:
                    continue
                orf_id = parent.split(",", 1)[0]
                cds_by_orf[orf_id].append((int(start), int(end), strand, phase))
                orfs.setdefault(orf_id, {})
                orfs[orf_id].setdefault("transcript_id", seqid)
                orfs[orf_id].setdefault("strand", strand)

    for orf_id, info in orfs.items():
        cds_parts = cds_by_orf.get(orf_id, [])
        if not cds_parts:
            raise ValueError(f"No CDS features found for ORF {orf_id} in {gff3_path}")
        info["cds_parts"] = cds_parts

    return orfs


def ordered_cds_parts(cds_parts: list[tuple[int, int, str, str]], strand: str) -> list[tuple[int, int, str, str]]:
    return sorted(cds_parts, key=lambda x: x[0], reverse=(strand == "-"))


def extract_coding_sequence(transcript_seq: Seq, cds_parts: list[tuple[int, int, str, str]], strand: str) -> Seq:
    chunks: list[str] = []
    for start, end, _, phase in ordered_cds_parts(cds_parts, strand):
        chunk = transcript_seq[start - 1 : end]
        if strand == "-":
            chunk = chunk.reverse_complement()

        if phase in {"0", "1", "2"}:
            chunk = chunk[int(phase) :]

        chunks.append(str(chunk))

    cds_seq = Seq("".join(chunks))
    if not cds_seq:
        raise ValueError("Extracted empty CDS sequence")
    return cds_seq


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract peptide FASTA from a TransDecoder final GFF3 and transcript FASTA."
    )
    parser.add_argument("--gff3", required=True, help="Final TransDecoder GFF3")
    parser.add_argument("--transcripts", required=True, help="Transcript FASTA used by TransDecoder")
    parser.add_argument("--out_pep", required=True, help="Output peptide FASTA")
    args = parser.parse_args()

    transcripts = load_transcripts(args.transcripts)
    orfs = load_orfs(args.gff3)

    with open(args.out_pep, "w") as out_handle:
        for orf_id in sorted(orfs):
            transcript_id = orfs[orf_id]["transcript_id"]
            strand = orfs[orf_id]["strand"]
            cds_parts = orfs[orf_id]["cds_parts"]

            if transcript_id not in transcripts:
                raise ValueError(
                    f"Transcript {transcript_id} referenced by ORF {orf_id} is missing from {args.transcripts}"
                )

            cds_seq = extract_coding_sequence(transcripts[transcript_id], cds_parts, strand)
            pep_seq = cds_seq.translate(to_stop=False)
            pep_text = str(pep_seq)
            if pep_text.endswith("*"):
                pep_text = pep_text[:-1]

            if not pep_text:
                raise ValueError(f"Extracted empty peptide sequence for ORF {orf_id}")

            out_handle.write(f">{orf_id}\n{pep_text}\n")


if __name__ == "__main__":
    main()
