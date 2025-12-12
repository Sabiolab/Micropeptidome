#!/usr/bin/env python3

# check de novo transcripts GTF for novel transcripts without annotation
# this file just prints the counts of annotated vs novel transcripts, to 
# give a broad overview of the data
import re
from collections import defaultdict

file_path = "merged.gtf" # path to the merged .gtf file that contains all the transcripts : merged.gtf

total_tx = 0
annotated_tx = 0
novel_tx = 0

novel_examples = []

with open(file_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = cols
        if feature != "transcript":
            continue
        
        total_tx += 1
        
        # parse attributes
        attr_dict = {}
        for field in attrs.split(";"):
            field = field.strip()
            if not field:
                continue
            # GTF: key "value"
            if " " in field:
                key, val = field.split(" ", 1)
                attr_dict[key] = val.strip().strip('"')
            elif "=" in field:  # just in case of GFF style
                key, val = field.split("=", 1)
                attr_dict[key] = val.strip().strip('"')
        
        tx_id = attr_dict.get("transcript_id", None)
        ref_gene_id = attr_dict.get("ref_gene_id", None)
        
        if ref_gene_id is None:
            novel_tx += 1
            if len(novel_examples) < 20:
                novel_examples.append({
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "transcript_id": tx_id
                })
        else:
            annotated_tx += 1

print("Total transcripts:", total_tx)
print("Annotated transcripts:", annotated_tx)
print("Novel transcripts:", novel_tx)
print("Example novel transcripts:")
for ex in novel_examples[:3]:
    print(ex)
