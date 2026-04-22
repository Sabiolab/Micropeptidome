STRAND = config.get("strandedness", "none").lower()

if STRAND not in {"none", "forward", "reverse"}:
    raise ValueError(f"Invalid strandedness: {STRAND}")

STRINGTIE_STRAND_FLAG = {
    "none": "",
    "forward": "--fr",
    "reverse": "--rf",
}[STRAND]

STRINGTIE_NASCENT_FLAG = "-N" if config.get("stringtie_rRNA", False) else ""

rule stringtie_assemble:
    input:
        bam=bam_path,
        ref_gtf=config["genome_gtf"]
    output:
        gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/stringtie/{{sample}}.gtf"
    threads: config.get("threads_stringtie", 8)
    resources:
        mem_mb=16000,
        runtime=120
    params:
        strand_flag=STRINGTIE_STRAND_FLAG,
        nascent_flag=STRINGTIE_NASCENT_FLAG,
        min_len_nt=int(config.get("stringtie_min_transcript_nt", 150))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/stringtie"

        stringtie "{input.bam}" \
          -G "{input.ref_gtf}" \
          -o "{output.gtf}" \
          -p {threads} \
          {params.nascent_flag} \
          {params.strand_flag} \
          -m {params.min_len_nt}
        """


rule stringtie_merge_cohort:
    input:
        gtfs=stringtie_gtfs_for_cohort,
        ref_gtf=config["genome_gtf"]
    output:
        gtf=cohort_stringtie_merge_gtf()
    threads: config.get("threads_stringtie", 8)
    resources:
        mem_mb=16000,
        runtime=120
    params:
        merge_dir=f"{STRINGTIE_MERGE_DIR}/{COHORT_LABEL}",
        gtf_list=f"{STRINGTIE_MERGE_DIR}/{COHORT_LABEL}/{COHORT_LABEL}.assemblies.txt"
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.merge_dir}"

        printf '%s\n' {input.gtfs:q} > "{params.gtf_list}"

        stringtie --merge \
          -G "{input.ref_gtf}" \
          -o "{output.gtf}" \
          -p {threads} \
          "{params.gtf_list}"
        """
