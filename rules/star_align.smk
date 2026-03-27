# STAR genome index + alignment
# Requires in main Snakefile:
#   - OUTDIR, STAR_INDEX_DIR
#   - fastq_r1(), fastq_r2()
#   - config["genome_fa"], config["genome_gtf"]

rule star_genome_index:
    input:
        genome=config["genome_fa"],
        gtf=config["genome_gtf"]
    output:
        # Use a STAR-produced file as the sentinel so users can supply a prebuilt index dir.
        genome=f"{STAR_INDEX_DIR}/Genome"
    threads: max(2, int(config.get("threads_star_index", 8)))
    resources:
        mem_mb=int(config.get("mem_star_index_mb", 64000)),
        runtime=int(config.get("runtime_star_index_min", 720))
    conda:
        "../envs/STAR.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{STAR_INDEX_DIR}"

        # Build STAR genome index
        STAR \
          --runMode genomeGenerate \
          --runThreadN {threads} \
          --genomeDir "{STAR_INDEX_DIR}" \
          --genomeFastaFiles "{input.genome}" \
          --sjdbGTFfile "{input.gtf}" \
          --sjdbOverhang 100
        """

rule star_align:
    input:
        r1=fastq_r1,
        r2=fastq_r2,
        index_genome=f"{STAR_INDEX_DIR}/Genome"
    output:
        bam=f"{OUTDIR}/star/{{sample}}.aligned.bam",
        bai=f"{OUTDIR}/star/{{sample}}.aligned.bam.bai",
        done=f"{OUTDIR}/star/{{sample}}.star_align.done"
    threads: max(2, int(config.get("threads_star_align", 12)))
    resources:
        mem_mb=int(config.get("mem_star_align_mb", 64000)),
        runtime=int(config.get("runtime_star_align_min", 360))
    conda:
        "../envs/STAR.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/star"

        tmpdir="{OUTDIR}/star/{wildcards.sample}.star_tmp"
        mkdir -p "$tmpdir"
        sort_mem_mb=$(( ({resources.mem_mb} * 4 / 10) / {threads} ))
        if [ "$sort_mem_mb" -lt 1024 ]; then
          sort_mem_mb=1024
        fi

        # StringTie needs XS tags, but external sorting gives us explicit memory control.
        STAR \
          --runThreadN {threads} \
          --genomeDir "{STAR_INDEX_DIR}" \
          --readFilesIn "{input.r1}" "{input.r2}" \
          --readFilesCommand zcat \
          --outFileNamePrefix "$tmpdir/" \
          --outSAMtype BAM Unsorted \
          --outSAMstrandField intronMotif \
          --outSAMattributes All

        samtools sort \
          -@ {threads} \
          -m "${{sort_mem_mb}}M" \
          -T "$tmpdir/samtools_sort" \
          -o "{output.bam}" \
          "$tmpdir/Aligned.out.bam"

        # Index canonical BAM
        samtools index -@ {threads} "{output.bam}" "{output.bai}"

        # Cleanup
        rm -rf "$tmpdir"

        touch "{output.done}"
        """