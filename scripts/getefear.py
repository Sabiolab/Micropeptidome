#!/usr/bin/env python3

import pandas as pd
import sys
import argparse
from pathlib import Path

print(r"""
          _____           _   _     _                      
 _   _   |  __ \         | | (_)   | |                     
| | | |  | |__) |__ _ __ | |_ _  __| | ___  _ __ ___   ___ 
| | | |  |  ___/ _ \ '_ \| __| |/ _` |/ _ \| '_ ` _ \ / _ \
| |_| |  | |  |  __/ |_) | |_| | (_| | (_) | | | | | |  __/
| ___/   |_|   \___| .__/ \__|_|\__,_|\___/|_| |_| |_|\___|
| |                | |                                       
|_|                |_|                  almoraco 11/12/2025
""")

print(
    "GTF creator for Micropeptidome pipeline\n"
    "Using SmProt2_mouse_Ribo.txt as reference for genomic coordinates\n"
    "Compatible with ShortStop pipeline (requires CDS and transcript features)\n"
    "Logs Protein_IDs not found in SmProt for later processing\n"
)

# Parse arguments
parser = argparse.ArgumentParser(description='Generate GTF for ShortStop from microproteins')
parser.add_argument('--genome', type=str, help='Reference genome FASTA file (optional, for chromosome name detection)')
parser.add_argument('--strip-chr', action='store_true', help='Remove "chr" prefix from chromosome names')
parser.add_argument('--add-chr', action='store_true', help='Add "chr" prefix to chromosome names')
args = parser.parse_args()

# Detect chromosome format from genome if provided
chr_prefix = None
if args.genome:
    print(f"🔍 Detecting chromosome format from genome: {args.genome}")
    try:
        with open(args.genome, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('>'):
                chr_name = first_line[1:].split()[0]  # Get first word after >
                if chr_name.startswith('chr'):
                    chr_prefix = True
                    print(f"   ✓ Genome uses 'chr' prefix (e.g., {chr_name})")
                else:
                    chr_prefix = False
                    print(f"   ✓ Genome does NOT use 'chr' prefix (e.g., {chr_name})")
    except Exception as e:
        print(f"   ⚠️  Could not read genome file: {e}")

# Override with manual flags if provided
if args.strip_chr:
    chr_prefix = False
    print("   🔧 Manual override: REMOVING 'chr' prefix")
elif args.add_chr:
    chr_prefix = True
    print("   🔧 Manual override: ADDING 'chr' prefix")

# ================
# Load microproteins
df_micro = pd.read_csv("data/microproteinasSEER.csv", sep="\t", encoding="utf-8-sig")
df_micro.columns = df_micro.columns.str.strip()

print(f"\n📊 Loaded microproteins CSV:")
print(f"   Total rows: {len(df_micro)}")
print(f"   Unique Protein_IDs: {df_micro['Protein_ID'].nunique()}")

# Drop duplicates, keep highest bitscore
n_duplicates = len(df_micro) - df_micro['Protein_ID'].nunique()
if n_duplicates > 0:
    print(f"   ⚠️  {n_duplicates} duplicate Protein_IDs detected")
    print(f"   Keeping entry with highest bitscore per Protein_ID\n")
    
    if 'bitscore' in df_micro.columns:
        df_micro = df_micro.sort_values('bitscore', ascending=False)
    df_micro = df_micro.drop_duplicates(subset=['Protein_ID'], keep='first')
    
    print(f"   After deduplication: {len(df_micro)} unique proteins\n")

# Load SmProt
df_smprot = pd.read_csv(
    "data/SmProt_v2_0/SmProt2_mouse_Ribo.txt",
    sep="\t",
    encoding="utf-8-sig"
)
df_smprot.columns = df_smprot.columns.str.strip()

print(f"📊 Loaded SmProt reference:")
print(f"   Total entries: {len(df_smprot)}")
print(f"   Unique SmProtIDs: {df_smprot['SmProtID'].nunique()}\n")

# Merge
# ================
df_merged = df_micro.merge(
    df_smprot,
    left_on="Protein_ID",
    right_on="SmProtID",
    how="left",
    indicator=True
)

# Missing
df_missing = df_merged[df_merged["_merge"] == "left_only"]
if not df_missing.empty:
    missing_file = "data/microproteins_not_in_SmProt.csv"
    df_missing[["Protein_ID"]].to_csv(missing_file, index=False)
    print(f"⚠️  {len(df_missing)} Protein_ID(s) not found in SmProt v2.0")
    print(f"   Saved to: {missing_file}")
    print(f"   Sample missing IDs: {df_missing['Protein_ID'].head(5).tolist()}\n")

# Found
df_found = df_merged[df_merged["_merge"] == "both"].copy()
if df_found.empty:
    print("\n❌ ERROR: No matching proteins found!")
    sys.exit(1)

print(f"✅ Successfully matched: {len(df_found)} microproteins\n")

# Coordinate cleaning
# ================
print("🔧 Processing genomic coordinates...")
df_found = df_found.dropna(subset=['Start', 'Stop'])
df_found["Start"] = df_found["Start"].astype(int)
df_found["Stop"] = df_found["Stop"].astype(int)
df_found['size'] = (df_found['Stop'] - df_found['Start']).abs()

# Filter out entries with invalid coordinates
df_found = df_found[df_found['Start'] != df_found['Stop']]
df_found = df_found[df_found['size'] > 0]

print(f"   Entries with valid coordinates: {len(df_found)}")
print(f"   Size range: {df_found['size'].min()} - {df_found['size'].max()} bp\n")

# Adjust chromosome names based on detection or flags
def adjust_chr_name(chr_name, add_prefix):
    """Add or remove 'chr' prefix from chromosome name"""
    chr_str = str(chr_name).strip()
    
    if add_prefix is None:
        # No change
        return chr_str
    elif add_prefix:
        # Add 'chr' if not present
        if not chr_str.startswith('chr'):
            return f"chr{chr_str}"
        return chr_str
    else:
        # Remove 'chr' if present
        if chr_str.startswith('chr'):
            return chr_str[3:]
        return chr_str

# Apply chromosome name adjustment
if chr_prefix is not None:
    print(f"🔧 Adjusting chromosome names...")
    original_example = df_found['Chr'].iloc[0] if len(df_found) > 0 else 'N/A'
    df_found['Chr'] = df_found['Chr'].apply(lambda x: adjust_chr_name(x, chr_prefix))
    adjusted_example = df_found['Chr'].iloc[0] if len(df_found) > 0 else 'N/A'
    print(f"   Example: {original_example} → {adjusted_example}\n")

# ================
# Generate GTF
print("📝 Generating ShortStop-compatible GTF file...")

output_gtf = "data/microproteins.gtf"

with open(output_gtf, 'w') as f:
    # NO headers - ShortStop's pandas parser can't handle them
    # Start directly with GTF data
    
    for idx, row in df_found.iterrows():
        # Extract required fields
        seqname = str(row['Chr']).strip()
        source = "SmProt_v2"
        start = int(row['Start'])
        stop = int(row['Stop'])
        strand = str(row['Strand']).strip() if 'Strand' in row else '+'
        
        # Ensure strand is valid
        if strand not in ['+', '-']:
            strand = '+'
        
        # Ensure proper coordinate order (GTF uses 1-based, inclusive coordinates)
        if start > stop:
            start, stop = stop, start
        
        # Create unique IDs
        transcript_id = f"{row['Protein_ID']}"
        gene_id = f"gene_{row['Protein_ID']}"
        
        # GTF requires 9 fields: seqname, source, feature, start, end, score, strand, frame, attributes
        # ShortStop needs both transcript and CDS features
        
        # Build attributes string
        attributes_transcript = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_name "{row["Protein_ID"]}";'
        attributes_cds = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_name "{row["Protein_ID"]}"; protein_id "{row["Protein_ID"]}";'
        
        # Write transcript feature
        f.write(f"{seqname}\t{source}\ttranscript\t{start}\t{stop}\t.\t{strand}\t.\t{attributes_transcript}\n")
        
        # Write CDS feature (ShortStop requires this)
        # Frame is typically 0 for start codons, but can be calculated if needed
        frame = "0"
        f.write(f"{seqname}\t{source}\tCDS\t{start}\t{stop}\t.\t{strand}\t{frame}\t{attributes_cds}\n")

print(f"✅ GTF file created: {output_gtf}")
print(f"   Total entries: {len(df_found)} microproteins")
print(f"   Features per entry: transcript + CDS (2 lines per microprotein)")
print(f"   Total GTF lines: {len(df_found) * 2}\n")

# Summary statistics
print("📈 Summary Statistics:")
print(f"   Total microproteins processed: {len(df_micro)}")
print(f"   Successfully matched to SmProt: {len(df_found)}")
print(f"   Not found in SmProt: {len(df_missing)}")
print(f"   Success rate: {len(df_found)/len(df_micro)*100:.1f}%")

# Chromosome distribution
if 'Chr' in df_found.columns:
    print(f"\n   Chromosome distribution:")
    chr_counts = df_found['Chr'].value_counts().head(10)
    for chr_name, count in chr_counts.items():
        print(f"      {chr_name}: {count}")

print("\n✅ GTF ready for ShortStop pipeline.")
print(f"   Use: shortstop predict --genome <your_genome.fa> --putative_smorfs_gtf {output_gtf}")
