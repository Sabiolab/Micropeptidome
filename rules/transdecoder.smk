rule gffread_transcripts:
    input:
        gtf=cohort_stringtie_merge_gtf(),
        genome=config["genome_fa"]
    output:
        fa=f"{cohort_results_dir()}/transcripts/{COHORT_LABEL}.transcripts.fa"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=int(config.get("runtime_gffread_min", 60))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{cohort_results_dir()}/transcripts"

        gffread "{input.gtf}" \
          -g "{input.genome}" \
          -w "{output.fa}"
        """

rule transdecoder_longorfs:
    input:
        fa=f"{cohort_results_dir()}/transcripts/{COHORT_LABEL}.transcripts.fa"
    output:
        done=f"{cohort_results_dir()}/transdecoder/longorfs.done"
    threads: 1
    params:
        min_aa=config.get("min_aa", "30"),
        workdir=f"{cohort_results_dir()}/transdecoder",
        stem="td2_input",
        staged_fa="td2_input.fa",
        td2_output_dir="td2_input"
    resources:
        mem_mb=int(config.get("mem_transdecoder_longorfs_mb", 32000)),
        runtime=int(config.get("runtime_transdecoder_longorfs_min", 360))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.workdir}"
        mkdir -p "{params.workdir}/{params.td2_output_dir}"

        cp -f "{input.fa}" "{params.workdir}/{params.staged_fa}"

        cd "{params.workdir}"
        TD2.LongOrfs \
          -t "{params.staged_fa}" \
          -O "{params.td2_output_dir}" \
          -m "{params.min_aa}"

        touch "longorfs.done"
        """

rule transdecoder_predict:
    input:
        fa=f"{cohort_results_dir()}/transcripts/{COHORT_LABEL}.transcripts.fa",
        longorfs_done=f"{cohort_results_dir()}/transdecoder/longorfs.done"
    output:
        pep=f"{cohort_results_dir()}/transcripts/{COHORT_LABEL}.transcripts.fa.TD2.pep",
        gff3=f"{cohort_results_dir()}/transcripts/{COHORT_LABEL}.transcripts.fa.TD2.gff3"
    threads: 1
    params:
        transcripts_dir=f"{cohort_results_dir()}/transcripts",
        transcript_basename=f"{COHORT_LABEL}.transcripts.fa",
        td2_output_dir="../transdecoder/td2_input"
    resources:
        mem_mb=int(config.get("mem_transdecoder_predict_mb", 64000)),
        runtime=int(config.get("runtime_transdecoder_predict_min", 180))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        cd "{params.transcripts_dir}"
        TD2.Predict \
          -t "{params.transcript_basename}" \
          -O "{params.td2_output_dir}"

        test -s "{output.pep}"
        test -s "{output.gff3}"
        """
