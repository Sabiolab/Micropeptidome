#############################################
# Snakefile: cohort-wide smORF workflow
# Uses: Trim Galore, STAR, StringTie, TD2, ShortStop, BLASTP,
# Dr. Brendan Miller's smORF annotator, and cohort-wide RSEM summaries.
#############################################

import csv
import re
from pathlib import Path

configfile: "config.yaml"


# Absolute outdir avoids SLURM working-dir surprises.
OUTDIR = str(Path(config["outdir"]).resolve())


def resolve_pipeline_path(value: str) -> str:
    path = Path(str(value))
    return str(path if path.is_absolute() else (Path(OUTDIR) / path).resolve())


def resolve_repo_path(value: str) -> str:
    path = Path(str(value))
    return str(path if path.is_absolute() else (Path(workflow.basedir) / path).resolve())


def sanitize_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    label = label.strip("._-")
    if not label:
        raise ValueError(f"Label '{value}' becomes empty after sanitization.")
    return label


config["units_csv"] = resolve_repo_path(config["units_csv"])
metadata_config_value = config.get("sample_metadata_csv", config.get("sample_metadata_tsv"))
if metadata_config_value is None:
    raise ValueError("config must define 'sample_metadata_csv'.")
config["sample_metadata_csv"] = resolve_repo_path(metadata_config_value)
config["genome_fa"] = resolve_repo_path(config["genome_fa"])
config["genome_gtf"] = resolve_repo_path(config["genome_gtf"])
config["human_proteome_fa"] = resolve_repo_path(config["human_proteome_fa"])

TRIM_DIR = resolve_pipeline_path(config.get("trim_dir", f"{OUTDIR}/trim_galore"))
STAR_DIR = resolve_pipeline_path(config.get("star_dir", f"{OUTDIR}/star"))
STAR_INDEX_DIR = resolve_pipeline_path(config.get("star_index_dir", f"{OUTDIR}/star_index"))
USE_PREBUILT_STAR_INDEX = bool(config.get("use_prebuilt_star_index", False))
STRINGTIE_MERGE_DIR = resolve_pipeline_path(
    config.get("stringtie_merge_dir", f"{OUTDIR}/stringtie_merge")
)
RESULTS_SHORTSTOP_DIR = resolve_pipeline_path(
    config.get("results_shortstop_dir", f"{OUTDIR}/results_shortstop")
)
COHORT_RESULTS_ROOT = resolve_pipeline_path(
    config.get(
        "cohort_results_dir",
        config.get("condition_results_dir", f"{OUTDIR}/results_cohort"),
    )
)
COHORT_MERGED_DIR = resolve_pipeline_path(
    config.get(
        "cohort_merged_dir",
        config.get("condition_merged_dir", f"{OUTDIR}/merged_cohort"),
    )
)
RSEM_DIR = resolve_pipeline_path(config.get("rsem_dir", f"{OUTDIR}/results_rsem_smorf"))
HUMAN_BLASTDB_PREFIX = resolve_pipeline_path(
    config.get("human_blastdb_prefix", f"{OUTDIR}/blastdb/human_proteome")
)
PATIENT_ID_COLUMN = str(config.get("patient_id_column", "PatientID"))
COHORT_PREFIX_RAW = str(config["cohort_prefix"])
COHORT_PREFIX = (
    COHORT_PREFIX_RAW
    if Path(COHORT_PREFIX_RAW).is_absolute()
    else str((Path(OUTDIR) / COHORT_PREFIX_RAW).resolve())
)
COHORT_LABEL = sanitize_label(Path(COHORT_PREFIX).name)


def trimmed_r1(sample: str) -> str:
    return f"{TRIM_DIR}/{sample}/{sample}_val_1.fq.gz"


def trimmed_r2(sample: str) -> str:
    return f"{TRIM_DIR}/{sample}/{sample}_val_2.fq.gz"


def star_bam(sample: str) -> str:
    return f"{STAR_DIR}/{sample}/{sample}.Aligned.sortedByCoord.out.bam"


def star_bai(sample: str) -> str:
    return f"{star_bam(sample)}.bai"


def stringtie_gtf(sample: str) -> str:
    return f"{RESULTS_SHORTSTOP_DIR}/{sample}/stringtie/{sample}.gtf"


def cohort_results_dir() -> str:
    return f"{COHORT_RESULTS_ROOT}/{COHORT_LABEL}"


def cohort_shortstop_done() -> str:
    return f"{cohort_results_dir()}/shortstop/predict.done"


def cohort_shortstop_gtf() -> str:
    return f"{cohort_results_dir()}/shortstop/{COHORT_LABEL}.smorfs_shortstop.gtf"


def cohort_merged_csv() -> str:
    return f"{COHORT_MERGED_DIR}/{COHORT_LABEL}.merged.csv"


def cohort_stringtie_merge_gtf() -> str:
    return f"{STRINGTIE_MERGE_DIR}/{COHORT_LABEL}/{COHORT_LABEL}.merged.gtf"


def cohort_all_loci() -> str:
    return f"{COHORT_PREFIX}.all_loci.csv"


def cohort_all_loci_tpm() -> str:
    return f"{COHORT_PREFIX}.all_loci.with_tpms.csv"


def cohort_all_loci_blast() -> str:
    return f"{COHORT_PREFIX}.all_loci.with_tpms.blastp_human.csv"


def cohort_tximport_rds() -> str:
    return f"{COHORT_PREFIX}.all_loci.blastp_human.tximport.rds"


def cohort_rsem_dir() -> str:
    return f"{RSEM_DIR}/{COHORT_LABEL}"


def cohort_rsem_ref_dir() -> str:
    return f"{cohort_rsem_dir()}/reference"


def cohort_rsem_ref_prefix() -> str:
    return f"{cohort_rsem_ref_dir()}/smorfs"


def sample_rsem_bam(sample: str) -> str:
    return f"{cohort_rsem_dir()}/{sample}/{sample}.bowtie2.bam"


def sample_rsem_log(sample: str) -> str:
    return f"{cohort_rsem_dir()}/{sample}/{sample}.bowtie2.log"


def sample_rsem_isoforms(sample: str) -> str:
    return f"{cohort_rsem_dir()}/{sample}/{sample}.isoforms.results"


def sample_rsem_genes(sample: str) -> str:
    return f"{cohort_rsem_dir()}/{sample}/{sample}.genes.results"


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

SAMPLES = sorted(UNITS)
if not SAMPLES:
    raise ValueError(f"No samples found in {config['units_csv']}. Check units.csv.")

metadata_by_patient = {}
with open(config["sample_metadata_csv"], newline="") as fh:
    reader = csv.DictReader(fh, delimiter=",")
    if reader.fieldnames is None:
        raise ValueError(f"{config['sample_metadata_csv']} is empty.")
    if PATIENT_ID_COLUMN not in reader.fieldnames:
        raise ValueError(
            f"{config['sample_metadata_csv']} must contain a '{PATIENT_ID_COLUMN}' column."
        )
    for row in reader:
        patient = (row.get(PATIENT_ID_COLUMN) or "").strip()
        if not patient:
            continue
        if patient in metadata_by_patient:
            raise ValueError(
                f"Duplicate patient '{patient}' in {config['sample_metadata_csv']}."
            )
        metadata_by_patient[patient] = row

missing_metadata = [sample for sample in SAMPLES if sample not in metadata_by_patient]
if missing_metadata:
    missing = ", ".join(missing_metadata[:10])
    raise ValueError(
        "Each sample in units.csv must have metadata in SampleMetadata.csv. "
        f"Missing patients: {missing}"
    )


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


def bam_path(wc):
    return star_bam(wc.sample)


def stringtie_gtfs_for_cohort(_wc):
    return [stringtie_gtf(sample) for sample in SAMPLES]


def rsem_isoforms_for_cohort(_wc):
    return [sample_rsem_isoforms(sample) for sample in SAMPLES]


include: "rules/trim_galore.smk"
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

ALL_TARGETS = [
    *expand(trimmed_r1("{sample}"), sample=SAMPLES),
    *expand(trimmed_r2("{sample}"), sample=SAMPLES),
    *expand(star_bam("{sample}"), sample=SAMPLES),
    *expand(star_bai("{sample}"), sample=SAMPLES),
    cohort_stringtie_merge_gtf(),
    cohort_shortstop_done(),
    cohort_shortstop_gtf(),
    cohort_merged_csv(),
    cohort_all_loci_blast(),
    cohort_tximport_rds(),
    *[sample_rsem_isoforms(sample) for sample in SAMPLES],
    *[sample_rsem_genes(sample) for sample in SAMPLES],
]


rule all:
    input:
        ALL_TARGETS
