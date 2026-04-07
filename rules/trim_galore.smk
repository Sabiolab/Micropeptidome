rule trim_galore_paired:
    input:
        r1=fastq_r1,
        r2=fastq_r2
    output:
        r1=trimmed_fastq_r1,
        r2=trimmed_fastq_r2
    threads: max(1, int(config.get("threads_trim_galore", 4)))
    resources:
        mem_mb=int(config.get("mem_trim_galore_mb", 8000)),
        runtime=int(config.get("runtime_trim_galore_min", 180))
    conda:
        "../envs/trim_galore.yaml"
    shell:
        r"""
        set -euo pipefail
        outdir="{TRIM_DIR}/{wildcards.sample}"
        mkdir -p "$outdir"

        trim_galore \
          --paired \
          --gzip \
          --cores {threads} \
          --basename "{wildcards.sample}" \
          --output_dir "$outdir" \
          "{input.r1}" "{input.r2}"

        test -s "{output.r1}"
        test -s "{output.r2}"
        """
