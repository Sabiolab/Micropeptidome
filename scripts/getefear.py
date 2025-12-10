#!/usr/bin/env python3
"""
GTF creator for Micropeptidome pipeline
Converts microproteins CSV to GTF format for smORF classification
Using SmProt2_mouse_Ribo.txt as reference for genomic coordinates
Logs Protein_IDs not found in SmProt for later processing
"""

import pandas as pd

print(r"""
          _____           _   _     _                      
 _   _   |  __ \         | | (_)   | |                     
| | | |  | |__) |__ _ __ | |_ _  __| | ___  _ __ ___   ___ 
| | | |  |  ___/ _ \ '_ \| __| |/ _` |/ _ \| '_ ` _ \ / _ \
| |_| |  | |  |  __/ |_) | |_| | (_| | (_) | | | | | |  __/
| ___/   |_|   \___| .__/ \__|_|\__,_|\___/|_| |_| |_|\___|
| |                | |                                       
|_|                |_|                  almoraco 10/12/2025
""")

# Leer CSV de microproteínas
df_micro = pd.read_csv("data/microproteinasSEER.csv", sep="\t", encoding="utf-8-sig")
df_micro.columns = df_micro.columns.str.strip()

# Leer SmProt2_mouse_Ribo.txt
df_smprot = pd.read_csv("data/SmProt_v2_0/SmProt2_mouse_Ribo.txt", sep="\t", encoding="utf-8-sig")
df_smprot.columns = df_smprot.columns.str.strip()

# Hacer merge usando Protein_ID = SmProtID
df_merged = df_micro.merge(df_smprot, left_on="Protein_ID", right_on="SmProtID", how="left", indicator=True)

# Proteínas no encontradas en SmProt
df_missing = df_merged[df_merged["_merge"] == "left_only"]
if not df_missing.empty:
    df_missing[["Protein_ID"]].to_csv("data/proteins_not_in_SmProt.csv", index=False)
    print(
    f"{df_missing['Protein_ID'].nunique()} unique Protein_ID(s) not found in SmProt v2.0.\n"
    f"List saved as data/proteins_not_in_SmProt.csv"
)


# Filtrar solo las que sí se encuentran
df_merged_found = df_merged[df_merged["_merge"] == "both"]

# Crear GTF
gtf = pd.DataFrame({
    "seqname": df_merged_found["Chr"],
    "source": "microproteomics",
    "feature": "microprotein_peptide",
    "start": df_merged_found["Start"],
    "end": df_merged_found["Stop"],
    "score": df_merged_found.get("bitscore", "."),
    "strand": df_merged_found["Strand"],
    "frame": ".",
    "attribute": (
        'gene_id "' + df_merged_found["Protein_ID"] + '"; '
        'peptide "' + df_merged_found["microprotein-derived peptide"] + '"; '
        'qseq "' + df_merged_found["qseq"] + '"; '
        'sseq "' + df_merged_found["sseqid"] + '";'
    )
})

gtf.to_csv("data/microproteins.gtf", sep="\t", index=False, header=False)

print("GTF generated at data/microproteins.gtf")
