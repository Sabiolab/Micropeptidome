STRAND = config.get("strandedness", "none").lower()

if STRAND not in {"none", "forward", "reverse"}:
    raise ValueError(f"Invalid strandedness: {STRAND}")

STRINGTIE_STRAND_FLAG = {
    "none": "",
    "forward": "--fr",
    "reverse": "--rf",
}[STRAND]

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
        strand_flag=STRINGTIE_STRAND_FLAG
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
          -m 100 \
          -c 1.5 \
          -f 0.005 \
          {params.strand_flag}
        """

MERGED_GTF = f"{OUTDIR}/stringtie_merged/merged.gtf"

rule stringtie_merge:
    input:
        gtfs=expand(f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/stringtie/{{sample}}.gtf", sample=SAMPLES),
        ref_gtf=config["genome_gtf"]
    output:
        gtf=MERGED_GTF
    threads: config.get("threads_stringtie", 8)
    resources:
        mem_mb=16000,
        runtime=120
    params:
        strand_flag=STRINGTIE_STRAND_FLAG
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/stringtie_merged"

        stringtie --merge \
          -G "{input.ref_gtf}" \
          -o "{output.gtf}" \
          -m 100 \
          -c 1.5 \
          -f 0.005 \
          -i \
          {params.strand_flag} \
          {input.gtfs}

        test -s "{output.gtf}"
        """