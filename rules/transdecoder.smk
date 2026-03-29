rule gffread_transcripts:
    input:
        gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/stringtie/{{sample}}.gtf",
        genome=config["genome_fa"]
    output:
        fa=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=int(config.get("runtime_gffread_min", 60))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/transcripts"

        gffread "{input.gtf}" \
          -g "{input.genome}" \
          -w "{output.fa}"
        """

rule transdecoder_longorfs:
    input:
        fa=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa"
    output:
        done=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transdecoder/longorfs.done"
    threads: 1
    params:
        min_aa=config.get("min_aa", "100")
    resources:
        mem_mb=16000,
        runtime=int(config.get("runtime_transdecoder_longorfs_min", 360))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/transdecoder"

        cd "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/transcripts"
        TD2.LongOrfs -t "{wildcards.sample}.transcripts.fa" -m "{params.min_aa}"

        touch "../transdecoder/longorfs.done"
        """

rule transdecoder_predict:
    input:
        fa=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa",
        longorfs_done=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transdecoder/longorfs.done"
    output:
        pep=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa.TD2.pep",
        gff3=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa.TD2.gff3"
    threads: 1
    resources:
        mem_mb=24000,
        runtime=int(config.get("runtime_transdecoder_predict_min", 180))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        cd "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/transcripts"

        TD2.Predict -t "{wildcards.sample}.transcripts.fa"

        test -s "{wildcards.sample}.transcripts.fa.TD2.pep"
        test -s "{wildcards.sample}.transcripts.fa.TD2.gff3"
        """
