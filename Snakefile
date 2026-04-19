#############################################
# Snakefile: per-patient ShortStop workflow
# Uses: Trim Galore, STAR, StringTie, TD2, ShortStop, BLASTP,
# Dr. Brendan Miller's smORF annotator, and condition-specific RSEM summaries.
#############################################

import csv
import re
from collections import defaultdict
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


def sanitize_condition_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    label = label.strip("._-")
    if not label:
        raise ValueError(f"Condition label '{value}' becomes empty after sanitization.")
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
STRINGTIE_MERGE_DIR = resolve_pipeline_path(
    config.get("stringtie_merge_dir", f"{OUTDIR}/stringtie_merge")
)
RESULTS_SHORTSTOP_DIR = resolve_pipeline_path(
    config.get("results_shortstop_dir", f"{OUTDIR}/results_shortstop")
)
CONDITION_RESULTS_DIR = resolve_pipeline_path(
    config.get("condition_results_dir", f"{OUTDIR}/results_conditions")
)
CONDITION_MERGED_DIR = resolve_pipeline_path(
    config.get("condition_merged_dir", f"{OUTDIR}/merged_per_condition")
)
RSEM_DIR = resolve_pipeline_path(config.get("rsem_dir", f"{OUTDIR}/results_rsem_smorf"))
HUMAN_BLASTDB_PREFIX = resolve_pipeline_path(
    config.get("human_blastdb_prefix", f"{OUTDIR}/blastdb/human_proteome")
)
PATIENT_ID_COLUMN = str(config.get("patient_id_column", "PatientID"))
CONDITION_COLUMN = str(config.get("condition_column", "Condition"))
COHORT_PREFIX_RAW = str(config["cohort_prefix"])
COHORT_PREFIX = (
    COHORT_PREFIX_RAW
    if Path(COHORT_PREFIX_RAW).is_absolute()
    else str((Path(OUTDIR) / COHORT_PREFIX_RAW).resolve())
)


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


def condition_results_dir(condition: str) -> str:
    return f"{CONDITION_RESULTS_DIR}/{condition}"


def condition_shortstop_done(condition: str) -> str:
    return f"{condition_results_dir(condition)}/shortstop/predict.done"


def condition_shortstop_gtf(condition: str) -> str:
    return f"{condition_results_dir(condition)}/shortstop/{condition}.smorfs_shortstop.gtf"


def condition_merged_csv(condition: str) -> str:
    return f"{CONDITION_MERGED_DIR}/{condition}.merged.csv"


def condition_prefix(condition: str) -> str:
    return f"{COHORT_PREFIX}.{condition}"


def condition_stringtie_merge_gtf(condition: str) -> str:
    return f"{STRINGTIE_MERGE_DIR}/{condition}/{condition}.merged.gtf"


def condition_all_loci(condition: str) -> str:
    return f"{condition_prefix(condition)}.all_loci.csv"


def condition_all_loci_tpm(condition: str) -> str:
    return f"{condition_prefix(condition)}.all_loci.with_tpms.csv"


def condition_all_loci_blast(condition: str) -> str:
    return f"{condition_prefix(condition)}.all_loci.with_tpms.blastp_human.csv"


def condition_all_loci_counts(condition: str) -> str:
    return f"{condition_prefix(condition)}.all_loci.blastp_human.expected_counts.tsv"


def condition_rsem_ref_dir(condition: str) -> str:
    return f"{RSEM_DIR}/{condition}/reference"


def condition_rsem_ref_prefix(condition: str) -> str:
    return f"{condition_rsem_ref_dir(condition)}/smorfs"


def condition_rsem_bam(condition: str, sample: str) -> str:
    return f"{RSEM_DIR}/{condition}/{sample}/{sample}.bowtie2.bam"


def condition_rsem_log(condition: str, sample: str) -> str:
    return f"{RSEM_DIR}/{condition}/{sample}/{sample}.bowtie2.log"


def condition_rsem_isoforms(condition: str, sample: str) -> str:
    return f"{RSEM_DIR}/{condition}/{sample}/{sample}.isoforms.results"


def condition_rsem_genes(condition: str, sample: str) -> str:
    return f"{RSEM_DIR}/{condition}/{sample}/{sample}.genes.results"


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
    if CONDITION_COLUMN not in reader.fieldnames:
        raise ValueError(
            f"{config['sample_metadata_csv']} must contain a '{CONDITION_COLUMN}' column."
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

condition_key_to_raw = {}
SAMPLE_TO_CONDITION = {}
CONDITION_TO_SAMPLES = defaultdict(list)
for sample in SAMPLES:
    row = metadata_by_patient[sample]
    raw_condition = (row.get(CONDITION_COLUMN) or "").strip()
    if not raw_condition:
        raise ValueError(
            f"Patient '{sample}' has an empty '{CONDITION_COLUMN}' value in "
            f"{config['sample_metadata_csv']}."
        )
    condition_key = sanitize_condition_label(raw_condition)
    existing = condition_key_to_raw.get(condition_key)
    if existing is not None and existing != raw_condition:
        raise ValueError(
            "Different raw condition labels collapse to the same sanitized key: "
            f"'{existing}' and '{raw_condition}' -> '{condition_key}'."
        )
    condition_key_to_raw[condition_key] = raw_condition
    SAMPLE_TO_CONDITION[sample] = condition_key
    CONDITION_TO_SAMPLES[condition_key].append(sample)

for condition in CONDITION_TO_SAMPLES:
    CONDITION_TO_SAMPLES[condition] = sorted(CONDITION_TO_SAMPLES[condition])

CONDITIONS = sorted(CONDITION_TO_SAMPLES)
CONDITION_SAMPLE_PAIRS = [
    (condition, sample)
    for condition in CONDITIONS
    for sample in CONDITION_TO_SAMPLES[condition]
]


def samples_for_condition(condition: str) -> list[str]:
    try:
        return CONDITION_TO_SAMPLES[condition]
    except KeyError as e:
        raise ValueError(f"Unknown condition '{condition}'.") from e


def validate_sample_condition(sample: str, condition: str) -> None:
    actual_condition = SAMPLE_TO_CONDITION.get(sample)
    if actual_condition != condition:
        raise ValueError(
            f"Sample '{sample}' belongs to condition '{actual_condition}', not '{condition}'."
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


def trimmed_fastq_r1_for_condition(wc):
    validate_sample_condition(wc.sample, wc.condition)
    return trimmed_r1(wc.sample)


def trimmed_fastq_r2_for_condition(wc):
    validate_sample_condition(wc.sample, wc.condition)
    return trimmed_r2(wc.sample)


def bam_path(wc):
    return star_bam(wc.sample)


def stringtie_gtfs_for_condition(wc):
    return [stringtie_gtf(sample) for sample in samples_for_condition(wc.condition)]


def rsem_isoforms_for_condition(wc):
    return [
        condition_rsem_isoforms(wc.condition, sample)
        for sample in samples_for_condition(wc.condition)
    ]


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
    *[condition_stringtie_merge_gtf(condition) for condition in CONDITIONS],
    *[condition_shortstop_done(condition) for condition in CONDITIONS],
    *[condition_shortstop_gtf(condition) for condition in CONDITIONS],
    *[condition_merged_csv(condition) for condition in CONDITIONS],
    *[condition_all_loci(condition) for condition in CONDITIONS],
    *[condition_all_loci_tpm(condition) for condition in CONDITIONS],
    *[condition_all_loci_blast(condition) for condition in CONDITIONS],
    *[condition_all_loci_counts(condition) for condition in CONDITIONS],
    *[
        condition_rsem_isoforms(condition, sample)
        for condition, sample in CONDITION_SAMPLE_PAIRS
    ],
    *[
        condition_rsem_genes(condition, sample)
        for condition, sample in CONDITION_SAMPLE_PAIRS
    ],
]


rule all:
    input:
        ALL_TARGETS
