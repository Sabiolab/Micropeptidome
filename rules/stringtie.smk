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
        min_len_nt=(int(config["min_aa"]) + 1) * 3
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


rule stringtie_merge_condition:
    input:
        gtfs=stringtie_gtfs_for_condition,
        ref_gtf=config["genome_gtf"]
    output:
        gtf=f"{STRINGTIE_MERGE_DIR}/{{condition}}/{{condition}}.merged.gtf"
    threads: config.get("threads_stringtie", 8)
    resources:
        mem_mb=16000,
        runtime=120
    params:
        gtf_list=lambda wc: f"{STRINGTIE_MERGE_DIR}/{wc.condition}/{wc.condition}.assemblies.txt"
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{STRINGTIE_MERGE_DIR}/{wildcards.condition}"

        printf '%s\n' {input.gtfs:q} > "{params.gtf_list}"

        stringtie --merge \
          -G "{input.ref_gtf}" \
          -o "{output.gtf}" \
          -p {threads} \
          "{params.gtf_list}"
        """
