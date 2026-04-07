#############################################
# Snakefile: per-sample ShortStop workflow
# Uses: Trim Galore, TopHat2 or StringTie-SuperReads, StringTie, TD2, ShortStop,
# BlastP, Dr. Brendan Miller's smORF annotator.
#############################################

import csv
from pathlib import Path

configfile: "config.yaml"

# Absolute outdir avoids SLURM working-dir surprises
OUTDIR = str(Path(config["outdir"]).resolve())

def resolve_pipeline_path(value: str) -> str:
    path = Path(str(value))
    return str(path if path.is_absolute() else (Path(OUTDIR) / path).resolve())

TRIM_DIR = resolve_pipeline_path(config.get("trim_dir", f"{OUTDIR}/trim_galore"))
TOPHAT_DIR = resolve_pipeline_path(config.get("tophat_dir", f"{OUTDIR}/tophat"))
TOPHAT_INDEX_PREFIX = resolve_pipeline_path(
    config.get("tophat_genome_index_prefix", f"{OUTDIR}/tophat_index/genome")
)
TOPHAT_TRANSCRIPTOME_PREFIX = resolve_pipeline_path(
    config.get("tophat_transcriptome_prefix", f"{OUTDIR}/tophat_index/transcriptome/known")
)
SUPERREADS_DIR = resolve_pipeline_path(config.get("superreads_dir", f"{OUTDIR}/superreads"))
SUPERREADS_INSTALL_DIR = resolve_pipeline_path(
    config.get("superreads_install_dir", f"{OUTDIR}/.deps/superreads")
)
SUPERREADS_HISAT2_INDEX_PREFIX = resolve_pipeline_path(
    config.get("superreads_hisat2_index_prefix", f"{OUTDIR}/superreads_index/hisat2/genome")
)
SUPERREADS_GMAP_DIR = resolve_pipeline_path(
    config.get("superreads_gmap_dir", f"{OUTDIR}/superreads_index/gmap")
)
SUPERREADS_GMAP_NAME = str(config.get("superreads_gmap_name", "genome"))
ASSEMBLY_BAM_SOURCE = str(config.get("stringtie_bam_source", "superreads")).lower()

if ASSEMBLY_BAM_SOURCE not in {"superreads", "tophat"}:
    raise ValueError("config['stringtie_bam_source'] must be one of: superreads, tophat")

def trimmed_r1(sample: str) -> str:
    return f"{TRIM_DIR}/{sample}/{sample}_val_1.fq.gz"

def trimmed_r2(sample: str) -> str:
    return f"{TRIM_DIR}/{sample}/{sample}_val_2.fq.gz"

def tophat_bam(sample: str) -> str:
    return f"{TOPHAT_DIR}/{sample}/accepted_hits.bam"

def superreads_bam(sample: str) -> str:
    return f"{SUPERREADS_DIR}/{sample}/sr_merge.bam"

def assembly_bam(sample: str) -> str:
    if ASSEMBLY_BAM_SOURCE == "tophat":
        return tophat_bam(sample)
    return superreads_bam(sample)

def bam_path(wc):
    return assembly_bam(wc.sample)

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

if ASSEMBLY_BAM_SOURCE == "superreads":
    ASSEMBLY_TARGETS = [
        expand(superreads_bam("{sample}"), sample=SAMPLES),
        expand(f"{superreads_bam('{sample}')}.bai", sample=SAMPLES),
    ]
else:
    ASSEMBLY_TARGETS = [
        expand(tophat_bam("{sample}"), sample=SAMPLES),
        expand(f"{tophat_bam('{sample}')}.bai", sample=SAMPLES),
    ]

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

def trimmed_fastq_r1(wc):
    return trimmed_r1(wc.sample)

def trimmed_fastq_r2(wc):
    return trimmed_r2(wc.sample)

include: "rules/trim_galore.smk"
include: "rules/tophat_align.smk"
include: "rules/superreads.smk"
include: "rules/stringtie.smk"
include: "rules/transdecoder.smk"
include: "rules/smorfs.smk"
include: "rules/shortstop.smk"
include: "rules/annotator.smk"
include: "rules/merge.smk"
include: "rules/rsem.smk"
include: "rules/blastp.smk"

### --------------------------------------------------oOo------------------------------------------------- ###

ALL_TARGETS = [
    *expand(trimmed_r1("{sample}"), sample=SAMPLES),
    *expand(trimmed_r2("{sample}"), sample=SAMPLES),
    *ASSEMBLY_TARGETS,
    *expand(f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/predict.done", sample=SAMPLES),
    *expand(f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/{{sample}}.smorfs_shortstop.gtf", sample=SAMPLES),
    *expand(f"{MERGED_DIR}/{{sample}}.merged.csv", sample=SAMPLES),
    f"{COHORT_PREFIX}.all_loci.csv",
    f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.csv",
    f"{COHORT_PREFIX}.all_loci.with_tpms.csv",
    f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.csv",
    f"{COHORT_PREFIX}.all_loci.with_tpms.blastp_human.csv",
    f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.blastp_human.csv",
    f"{COHORT_PREFIX}.all_loci.blastp_human.expected_counts.tsv",
    f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.blastp_human.expected_counts.tsv",
]

rule all:
    input:
        ALL_TARGETS
