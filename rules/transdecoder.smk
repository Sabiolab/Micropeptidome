rule gffread_transcripts:
    input:
        gtf=f"{STRINGTIE_MERGE_DIR}/{{condition}}/{{condition}}.merged.gtf",
        genome=config["genome_fa"]
    output:
        fa=f"{CONDITION_RESULTS_DIR}/{{condition}}/transcripts/{{condition}}.transcripts.fa"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=int(config.get("runtime_gffread_min", 60))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{CONDITION_RESULTS_DIR}/{wildcards.condition}/transcripts"

        gffread "{input.gtf}" \
          -g "{input.genome}" \
          -w "{output.fa}"
        """

rule transdecoder_longorfs:
    input:
        fa=f"{CONDITION_RESULTS_DIR}/{{condition}}/transcripts/{{condition}}.transcripts.fa"
    output:
        done=f"{CONDITION_RESULTS_DIR}/{{condition}}/transdecoder/longorfs.done"
    threads: 1
    params:
        min_aa=config.get("min_aa", "30")
    resources:
        mem_mb=16000,
        runtime=int(config.get("runtime_transdecoder_longorfs_min", 360))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{CONDITION_RESULTS_DIR}/{wildcards.condition}/transdecoder"

        cd "{CONDITION_RESULTS_DIR}/{wildcards.condition}/transcripts"
        TD2.LongOrfs -t "{wildcards.condition}.transcripts.fa" -m "{params.min_aa}"

        touch "../transdecoder/longorfs.done"
        """

rule transdecoder_predict:
    input:
        fa=f"{CONDITION_RESULTS_DIR}/{{condition}}/transcripts/{{condition}}.transcripts.fa",
        longorfs_done=f"{CONDITION_RESULTS_DIR}/{{condition}}/transdecoder/longorfs.done"
    output:
        pep=f"{CONDITION_RESULTS_DIR}/{{condition}}/transcripts/{{condition}}.transcripts.fa.TD2.pep",
        gff3=f"{CONDITION_RESULTS_DIR}/{{condition}}/transcripts/{{condition}}.transcripts.fa.TD2.gff3"
    threads: 1
    resources:
        mem_mb=24000,
        runtime=int(config.get("runtime_transdecoder_predict_min", 180))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        cd "{CONDITION_RESULTS_DIR}/{wildcards.condition}/transcripts"

        TD2.Predict -t "{wildcards.condition}.transcripts.fa"

        test -s "{wildcards.condition}.transcripts.fa.TD2.pep"
        test -s "{wildcards.condition}.transcripts.fa.TD2.gff3"
        """
