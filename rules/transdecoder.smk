rule gffread_transcripts:
    input:
        gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/stringtie/{{sample}}.gtf",
        genome=config["genome_fa"]
    output:
        fa=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
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
    resources:
        mem_mb=16000,
        runtime=60
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/transdecoder"

        cd "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/transcripts"
        TransDecoder.LongOrfs -t "{wildcards.sample}.transcripts.fa"

        touch "../transdecoder/longorfs.done"
        """

rule transdecoder_predict:
    input:
        fa=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa",
        longorfs_done=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transdecoder/longorfs.done"
    output:
        pep=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa.transdecoder.pep",
        gff3=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa.transdecoder.gff3"
    threads: 1
    resources:
        mem_mb=24000,
        runtime=120
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        tx_dir="$(dirname "{input.fa}")"
        tx_base="$(basename "{input.fa}")"
        td_dir="$tx_dir/${{tx_base}}.transdecoder_dir"

        # Remove stale Predict checkpoint state only
        rm -rf "$td_dir/__checkpoints_TDpredict"

        # Remove stale final outputs
        rm -f "{output.pep}" "{output.gff3}" \
              "$tx_dir/${{tx_base}}.transdecoder.cds" \
              "$tx_dir/${{tx_base}}.transdecoder.bed"

        # Remove stale final-selection intermediates
        rm -f "$td_dir/longest_orfs.cds.best_candidates.gff3" \
              "$td_dir/longest_orfs.cds.best_candidates.gff3.revised_starts.gff3"

        # Start-codon refinement is fragile on some transcript sets and can abort before
        TransDecoder.Predict \
          --no_refine_starts \
          -t "{input.fa}" \
          -O "$tx_dir"

        test -s "{output.pep}"
        test -s "{output.gff3}"
        """
