rule make_smorf_rsem_inputs:
    input:
        loci_csv=f"{COHORT_PREFIX}.all_loci.csv",
        script=lambda wc: config["make_smorf_rsem_ref_script"]
    output:
        fasta=f"{RSEM_REF_DIR}/smorfs.cds.fa",
        tx2gene=f"{RSEM_REF_DIR}/smorfs.tx2gene.tsv"
    resources:
        mem_mb=32000,
        runtime=240
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RSEM_REF_DIR}"

        python "{input.script}" \
          --loci_csv "{input.loci_csv}" \
          --fasta "{output.fasta}" \
          --tx2gene "{output.tx2gene}"

        test -s "{output.fasta}"
        test -s "{output.tx2gene}"
        """

rule rsem_prepare_smorf_reference:
    input:
        fasta=f"{RSEM_REF_DIR}/smorfs.cds.fa",
        tx2gene=f"{RSEM_REF_DIR}/smorfs.tx2gene.tsv"
    output:
        done=f"{RSEM_REF_DIR}/rsem_ref.done"
    resources:
        mem_mb=32000,
        runtime=240
    conda:
        "../envs/RSEM.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RSEM_REF_DIR}"

        rsem-prepare-reference \
          --transcript-to-gene-map "{input.tx2gene}" \
          --bowtie2 \
          "{input.fasta}" "{RSEM_REF_PREFIX}"

        touch "{output.done}"
        """

rule rsem_align_smorf_bowtie2:
    input:
        r1=fastq_r1,
        r2=fastq_r2,
        ref_done=f"{RSEM_REF_DIR}/rsem_ref.done"
    output:
        bam=f"{RSEM_DIR}/{{sample}}/{{sample}}.bowtie2.bam",
        log=f"{RSEM_DIR}/{{sample}}/{{sample}}.bowtie2.log"
    threads: config.get("threads_rsem_align", 8)
    resources:
        mem_mb=32000,
        runtime=600
    params:
        ref=RSEM_REF_PREFIX
    conda:
        "../envs/RSEM.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RSEM_DIR}/{wildcards.sample}"

        bowtie2 \
          --reorder \
          -q --phred33 --sensitive \
          --dpad 0 --gbar 99999999 \
          --mp 1,1 --np 1 \
          --score-min L,0,-0.1 \
          -I 1 -X 1000 \
          --no-mixed --no-discordant \
          -p {threads} -k 200 \
          -x "{params.ref}" \
          -1 "{input.r1}" -2 "{input.r2}" \
          2> "{output.log}" \
        | samtools view -b -o "{output.bam}" -

        test -s "{output.bam}"
        """

rule rsem_quant_smorf:
    input:
        bam=f"{RSEM_DIR}/{{sample}}/{{sample}}.bowtie2.bam",
        ref_done=f"{RSEM_REF_DIR}/rsem_ref.done"
    output:
        isoforms=f"{RSEM_DIR}/{{sample}}/{{sample}}.isoforms.results",
        genes=f"{RSEM_DIR}/{{sample}}/{{sample}}.genes.results"
    threads: config.get("threads_rsem_em", 8)
    resources:
        mem_mb=32000,
        runtime=600
    params:
        ref=RSEM_REF_PREFIX,
        stranded=config.get("strandedness", "none")
    conda:
        "../envs/RSEM.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RSEM_DIR}/{wildcards.sample}"

        rsem-calculate-expression \
          --paired-end \
          --alignments \
          -p {threads} \
          --strandedness "{params.stranded}" \
          "{input.bam}" \
          "{params.ref}" \
          "{RSEM_DIR}/{wildcards.sample}/{wildcards.sample}"

        test -s "{output.isoforms}"
        test -s "{output.genes}"
        """

rule add_rsem_tpms_to_locus_summary:
    input:
        all_loci=f"{COHORT_PREFIX}.all_loci.csv",
        shared=f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.csv",
        rsem_isoforms=expand(f"{RSEM_DIR}/{{sample}}/{{sample}}.isoforms.results", sample=SAMPLES),
        script=lambda wc: config["add_rsem_tpms_script"]
    output:
        all_loci_tpm=f"{COHORT_PREFIX}.all_loci.with_tpms.csv",
        shared_tpm=f"{COHORT_PREFIX}.shared_ge{MIN_PATIENTS}.with_tpms.csv"
    threads: 1
    resources:
        mem_mb=16000
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        python "{input.script}" \
          --all_loci_csv "{input.all_loci}" \
          --shared_csv "{input.shared}" \
          --rsem_dir "{RSEM_DIR}" \
          --out_all_loci_csv "{output.all_loci_tpm}" \
          --out_shared_csv "{output.shared_tpm}"
        """
