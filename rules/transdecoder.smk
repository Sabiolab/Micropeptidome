MERGED_TX_DIR  = f"{OUTDIR}/stringtie_merged"
MERGED_TX_FA   = f"{MERGED_TX_DIR}/merged.transcripts.fa"
MERGED_TD2_PEP = f"{MERGED_TX_DIR}/merged.transcripts.fa.TD2.pep"
MERGED_TD2_GFF = f"{MERGED_TX_DIR}/merged.transcripts.fa.TD2.gff3"

rule gffread_transcripts:
    input:
        gtf=MERGED_GTF,
        genome=config["genome_fa"]
    output:
        fa=MERGED_TX_FA
    threads: 1
    resources:
        mem_mb=8000,
        runtime=int(config.get("runtime_gffread_min", 60))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        gffread "{input.gtf}" \
          -g "{input.genome}" \
          -w "{output.fa}"

        test -s "{output.fa}"
        """

rule transdecoder_longorfs:
    input:
        fa=MERGED_TX_FA
    output:
        done=f"{MERGED_TX_DIR}/longorfs.done"
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
        cd "{MERGED_TX_DIR}"

        TD2.LongOrfs -t merged.transcripts.fa -m "{params.min_aa}"

        touch longorfs.done
        """

rule transdecoder_predict:
    input:
        fa=MERGED_TX_FA,
        longorfs_done=f"{MERGED_TX_DIR}/longorfs.done"
    output:
        pep=MERGED_TD2_PEP,
        gff3=MERGED_TD2_GFF
    threads: 1
    resources:
        mem_mb=24000,
        runtime=int(config.get("runtime_transdecoder_predict_min", 180))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        cd "{MERGED_TX_DIR}"

        TD2.Predict -t merged.transcripts.fa

        test -s merged.transcripts.fa.TD2.pep
        test -s merged.transcripts.fa.TD2.gff3
        """
