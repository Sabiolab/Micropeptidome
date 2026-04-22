rule make_smorf_rsem_inputs:
    input:
        loci_csv=cohort_all_loci(),
        script=config["make_smorf_rsem_ref_script"]
    output:
        fasta=f"{cohort_rsem_ref_dir()}/smorfs.cds.fa",
        tx2gene=f"{cohort_rsem_ref_dir()}/smorfs.tx2gene.tsv"
    resources:
        mem_mb=32000,
        runtime=240
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{cohort_rsem_ref_dir()}"

        python "{input.script}" \
          --loci_csv "{input.loci_csv}" \
          --fasta "{output.fasta}" \
          --tx2gene "{output.tx2gene}"

        test -s "{output.fasta}"
        test -s "{output.tx2gene}"
        """

rule rsem_prepare_smorf_reference:
    input:
        fasta=f"{cohort_rsem_ref_dir()}/smorfs.cds.fa",
        tx2gene=f"{cohort_rsem_ref_dir()}/smorfs.tx2gene.tsv"
    output:
        done=f"{cohort_rsem_ref_dir()}/rsem_ref.done"
    resources:
        mem_mb=32000,
        runtime=240
    conda:
        "../envs/RSEM.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{cohort_rsem_ref_dir()}"

        rsem-prepare-reference \
          --transcript-to-gene-map "{input.tx2gene}" \
          --bowtie2 \
          "{input.fasta}" "{cohort_rsem_ref_prefix()}"

        touch "{output.done}"
        """

rule rsem_align_smorf_bowtie2:
    input:
        r1=trimmed_fastq_r1,
        r2=trimmed_fastq_r2,
        ref_done=f"{cohort_rsem_ref_dir()}/rsem_ref.done"
    output:
        bam=f"{cohort_rsem_dir()}/{{sample}}/{{sample}}.bowtie2.bam",
        log=f"{cohort_rsem_dir()}/{{sample}}/{{sample}}.bowtie2.log"
    threads: config.get("threads_rsem_align", 8)
    resources:
        mem_mb=32000,
        runtime=600
    params:
        ref=cohort_rsem_ref_prefix()
    conda:
        "../envs/RSEM.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{cohort_rsem_dir()}/{wildcards.sample}"

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
        bam=lambda wc: sample_rsem_bam(wc.sample),
        ref_done=f"{cohort_rsem_ref_dir()}/rsem_ref.done"
    output:
        isoforms=f"{cohort_rsem_dir()}/{{sample}}/{{sample}}.isoforms.results",
        genes=f"{cohort_rsem_dir()}/{{sample}}/{{sample}}.genes.results"
    threads: config.get("threads_rsem_em", 8)
    resources:
        mem_mb=32000,
        runtime=600
    params:
        ref=cohort_rsem_ref_prefix(),
        stranded=config.get("strandedness", "none")
    conda:
        "../envs/RSEM.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{cohort_rsem_dir()}/{wildcards.sample}"

        rsem-calculate-expression \
          --paired-end \
          --alignments \
          -p {threads} \
          --strandedness "{params.stranded}" \
          "{input.bam}" \
          "{params.ref}" \
          "{cohort_rsem_dir()}/{wildcards.sample}/{wildcards.sample}"

        test -s "{output.isoforms}"
        test -s "{output.genes}"
        """

rule add_rsem_tpms_to_locus_summary:
    input:
        all_loci=cohort_all_loci(),
        rsem_isoforms=rsem_isoforms_for_cohort,
        script=config["add_rsem_tpms_script"]
    output:
        all_loci_tpm=temp(cohort_all_loci_tpm())
    threads: 1
    resources:
        mem_mb=16000
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        python "{input.script}" \
          --summary_csv "{input.all_loci}" \
          --rsem_dir "{cohort_rsem_dir()}" \
          --out_csv "{output.all_loci_tpm}"
        """

rule export_tximport_all_loci:
    input:
        loci_csv=cohort_all_loci_blast(),
        tx2gene=f"{cohort_rsem_ref_dir()}/smorfs.tx2gene.tsv",
        rsem_isoforms=rsem_isoforms_for_cohort,
        script=config["export_tximport_script"]
    output:
        tximport_rds=cohort_tximport_rds()
    threads: 1
    resources:
        mem_mb=16000
    conda:
        "../envs/tximport.yaml"
    shell:
        r"""
        set -euo pipefail
        Rscript "{input.script}" \
          --loci_csv "{input.loci_csv}" \
          --tx2gene "{input.tx2gene}" \
          --rsem_dir "{cohort_rsem_dir()}" \
          --out_rds "{output.tximport_rds}"
        """
