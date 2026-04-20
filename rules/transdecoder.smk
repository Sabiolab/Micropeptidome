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
        min_aa=config.get("min_aa", "30"),
        workdir=lambda wc: f"{CONDITION_RESULTS_DIR}/{wc.condition}/transdecoder",
        staged_fa="td2_input.fa"
    resources:
        mem_mb=int(config.get("mem_transdecoder_longorfs_mb", 32000)),
        runtime=int(config.get("runtime_transdecoder_longorfs_min", 360))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.workdir}"

        cp -f "{input.fa}" "{params.workdir}/{params.staged_fa}"

        cd "{params.workdir}"
        TD2.LongOrfs -t "{params.staged_fa}" -m "{params.min_aa}"

        touch "longorfs.done"
        """

rule transdecoder_predict:
    input:
        fa=f"{CONDITION_RESULTS_DIR}/{{condition}}/transcripts/{{condition}}.transcripts.fa",
        longorfs_done=f"{CONDITION_RESULTS_DIR}/{{condition}}/transdecoder/longorfs.done"
    output:
        pep=f"{CONDITION_RESULTS_DIR}/{{condition}}/transcripts/{{condition}}.transcripts.fa.TD2.pep",
        gff3=f"{CONDITION_RESULTS_DIR}/{{condition}}/transcripts/{{condition}}.transcripts.fa.TD2.gff3"
    threads: 1
    params:
        workdir=lambda wc: f"{CONDITION_RESULTS_DIR}/{wc.condition}/transdecoder",
        staged_fa="td2_input.fa"
    resources:
        mem_mb=int(config.get("mem_transdecoder_predict_mb", 64000)),
        runtime=int(config.get("runtime_transdecoder_predict_min", 180))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.workdir}"
        cp -f "{input.fa}" "{params.workdir}/{params.staged_fa}"

        cd "{params.workdir}"
        TD2.Predict -t "{params.staged_fa}"

        cp -f "{params.staged_fa}.TD2.pep" "{output.pep}"
        cp -f "{params.staged_fa}.TD2.gff3" "{output.gff3}"

        test -s "{output.pep}"
        test -s "{output.gff3}"
        """
