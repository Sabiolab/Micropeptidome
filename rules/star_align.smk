rule star_genome_index:
    input:
        genome=config["genome_fa"],
        gtf=config["genome_gtf"]
    output:
        sa=f"{STAR_INDEX_DIR}/SA"
    threads: max(1, int(config.get("threads_star_index", 8)))
    resources:
        mem_mb=int(config.get("mem_star_index_mb", 64000)),
        runtime=int(config.get("runtime_star_index_min", 360))
    params:
        sjdb_overhang=int(config.get("star_sjdb_overhang", 149))
    conda:
        "../envs/star.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{STAR_INDEX_DIR}"

        STAR \
          --runMode genomeGenerate \
          --runThreadN {threads} \
          --genomeDir "{STAR_INDEX_DIR}" \
          --genomeFastaFiles "{input.genome}" \
          --sjdbGTFfile "{input.gtf}" \
          --sjdbOverhang {params.sjdb_overhang}

        test -s "{output.sa}"
        """


rule star_align:
    input:
        r1=trimmed_fastq_r1,
        r2=trimmed_fastq_r2,
        index_sa=f"{STAR_INDEX_DIR}/SA"
    output:
        bam=f"{STAR_DIR}/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam",
        bai=f"{STAR_DIR}/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai"
    threads: max(1, int(config.get("threads_star_align", 8)))
    resources:
        mem_mb=int(config.get("mem_star_align_mb", 64000)),
        runtime=int(config.get("runtime_star_align_min", 720))
    params:
        out_prefix=lambda wc: f"{STAR_DIR}/{wc.sample}/{wc.sample}."
    conda:
        "../envs/star.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{STAR_DIR}/{wildcards.sample}"

        STAR \
          --runThreadN {threads} \
          --genomeDir "{STAR_INDEX_DIR}" \
          --readFilesIn "{input.r1}" "{input.r2}" \
          --readFilesCommand zcat \
          --outFileNamePrefix "{params.out_prefix}" \
          --outSAMtype BAM SortedByCoordinate

        samtools index -@ {threads} "{output.bam}" "{output.bai}"
        test -s "{output.bam}"
        test -s "{output.bai}"
        """
