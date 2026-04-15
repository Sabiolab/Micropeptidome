TOPHAT_LIBRARY_TYPE = {
    "none": "fr-unstranded",
    "forward": "fr-secondstrand",
    "reverse": "fr-firststrand",
}[config.get("strandedness", "none").lower()]


rule bowtie2_genome_index:
    input:
        genome=config["genome_fa"]
    output:
        idx1=f"{TOPHAT_INDEX_PREFIX}.1.bt2"
    threads: max(1, int(config.get("threads_tophat_index", 4)))
    resources:
        mem_mb=int(config.get("mem_tophat_index_mb", 32000)),
        runtime=int(config.get("runtime_tophat_index_min", 360))
    conda:
        "../envs/tophat.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{TOPHAT_INDEX_PREFIX}")"

        bowtie2-build \
          "{input.genome}" \
          "{TOPHAT_INDEX_PREFIX}"

        test -s "{output.idx1}"
        """


rule tophat_transcriptome_index:
    input:
        genome_index=f"{TOPHAT_INDEX_PREFIX}.1.bt2",
        gtf=config["genome_gtf"]
    output:
        done=f"{TOPHAT_TRANSCRIPTOME_PREFIX}.done"
    threads: max(1, int(config.get("threads_tophat_index", 4)))
    resources:
        mem_mb=int(config.get("mem_tophat_index_mb", 32000)),
        runtime=int(config.get("runtime_tophat_index_min", 360))
    params:
        outdir=lambda wc: f"{Path(TOPHAT_TRANSCRIPTOME_PREFIX).parent}/build"
    conda:
        "../envs/tophat.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{TOPHAT_TRANSCRIPTOME_PREFIX}")"
        mkdir -p "{params.outdir}"

        tophat \
          -p {threads} \
          -G "{input.gtf}" \
          --transcriptome-index "{TOPHAT_TRANSCRIPTOME_PREFIX}" \
          -o "{params.outdir}" \
          "{TOPHAT_INDEX_PREFIX}"

        touch "{output.done}"
        """


rule tophat2_align:
    input:
        r1=trimmed_fastq_r1,
        r2=trimmed_fastq_r2,
        genome_index=f"{TOPHAT_INDEX_PREFIX}.1.bt2",
        transcriptome_done=f"{TOPHAT_TRANSCRIPTOME_PREFIX}.done",
        gtf=config["genome_gtf"]
    output:
        bam=f"{TOPHAT_DIR}/{{sample}}/accepted_hits.bam",
        bai=f"{TOPHAT_DIR}/{{sample}}/accepted_hits.bam.bai"
    threads: max(1, int(config.get("threads_tophat_align", 8)))
    resources:
        mem_mb=int(config.get("mem_tophat_align_mb", 32000)),
        runtime=int(config.get("runtime_tophat_align_min", 720))
    params:
        outdir=lambda wc: f"{TOPHAT_DIR}/{wc.sample}",
        library_type=TOPHAT_LIBRARY_TYPE,
        mate_inner_dist=int(config.get("tophat_mate_inner_dist", 50)),
        mate_std_dev=int(config.get("tophat_mate_std_dev", 20))
    conda:
        "../envs/tophat.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"

        tophat \
          -p {threads} \
          -o "{params.outdir}" \
          -G "{input.gtf}" \
          --transcriptome-index "{TOPHAT_TRANSCRIPTOME_PREFIX}" \
          --library-type "{params.library_type}" \
          -r {params.mate_inner_dist} \
          --mate-std-dev {params.mate_std_dev} \
          "{TOPHAT_INDEX_PREFIX}" \
          "{input.r1}" "{input.r2}"

        samtools index -@ {threads} "{output.bam}" "{output.bai}"
        test -s "{output.bam}"
        test -s "{output.bai}"
        """
