#############################################
# Snakefile: per-sample ShortStop workflow
# Uses: StringTie, TransDecoder, ShortStop, 
# BlastP, Dr. Brendan Miller's smORF annotator.
#############################################

import csv
from pathlib import Path

configfile: "config.yaml"

# Absolute outdir avoids SLURM working-dir surprises
OUTDIR = str(Path(config["outdir"]).resolve())

# STAR: user may point this at a prebuilt STAR genome index directory.
STAR_INDEX_DIR = str(Path(config.get("star_index_dir", f"{OUTDIR}/STAR_index")).resolve())

# STAR output BAM path
def star_bam(sample: str) -> str:
    return f"{OUTDIR}/star/{sample}.aligned.bam"

def bam_path(wc):
    return star_bam(wc.sample)

RESULTS_SHORTSTOP_DIR = f"{OUTDIR}/results_shortstop"
COHORT_PREFIX_RAW = str(config["cohort_prefix"])
COHORT_PREFIX = (
    COHORT_PREFIX_RAW
    if Path(COHORT_PREFIX_RAW).is_absolute()
    else str((Path(OUTDIR) / COHORT_PREFIX_RAW).resolve())
)
MIN_PATIENTS = int(config.get("min_patients", 2))
MERGED_DIR = config.get("merged_dir", f"{OUTDIR}/merged_per_sample")
RSEM_DIR = config.get("rsem_dir", f"{OUTDIR}/results_rsem_smorf")
RSEM_REF_DIR = f"{RSEM_DIR}/reference"
RSEM_REF_PREFIX = f"{RSEM_REF_DIR}/smorfs"

# Build mapping: sample -> (r1, r2)
UNITS = {}
with open(config["units_csv"], newline="") as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        sample = (row.get("sample") or row.get("name") or "").strip()
        r1 = (row.get("r1") or row.get("fastq_r1") or row.get("read1") or "").strip()
        r2 = (row.get("r2") or row.get("fastq_r2") or row.get("read2") or "").strip()
        if not sample:
            continue
        if not r1 or not r2:
            raise ValueError(f"Missing r1/r2 for sample '{sample}' in {config['units_csv']}")
        if sample in UNITS:
            raise ValueError(f"Duplicate sample '{sample}' in {config['units_csv']}")
        UNITS[sample] = (r1, r2)

SAMPLES = sorted(UNITS.keys())
if not SAMPLES:
    raise ValueError(f"No samples found in {config['units_csv']}. Check units.csv.")

def fastq_r1(wc):
    try:
        return UNITS[wc.sample][0]
    except KeyError as e:
        raise ValueError(f"No FASTQ entry for sample '{wc.sample}' in {config['units_csv']}") from e

def fastq_r2(wc):
    try:
        return UNITS[wc.sample][1]
    except KeyError as e:
        raise ValueError(f"No FASTQ entry for sample '{wc.sample}' in {config['units_csv']}") from e

include: "rules/star_align.smk"
include: "rules/stringtie.smk"
include: "rules/transdecoder.smk"
include: "rules/smorfs.smk"
include: "rules/shortstop.smk"
include: "rules/annotator.smk"
include: "rules/merge.smk"
include: "rules/rsem.smk"
include: "rules/blastp.smk"

### --------------------------------------------------oOo------------------------------------------------- ###

rule all:
    input:
        expand(f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/predict.done", sample=SAMPLES),
        expand(f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/{{sample}}.smorfs_shortstop.gtf", sample=SAMPLES),
        expand(f"{MERGED_DIR}/{{sample}}.merged.csv", sample=SAMPLES),
        f"{COHORT_PREFIX}.all_loci.csv",
        f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.csv",
        f"{COHORT_PREFIX}.all_loci.with_tpms.csv",
        f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.csv",
        f"{COHORT_PREFIX}.all_loci.with_tpms.blastp_human.csv",
        f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.blastp_human.csv"
