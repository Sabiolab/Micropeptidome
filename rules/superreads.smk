SUPERREADS_REPO = str(
    config.get("superreads_repo", "https://github.com/gpertea/stringtie.git")
)
SUPERREADS_REF = str(config.get("superreads_ref", "v3.0.3"))

HISAT2_RNA_STRANDNESS = {
    "none": "",
    "forward": "FR",
    "reverse": "RF",
}[config.get("strandedness", "none").lower()]


rule hisat2_genome_index:
    input:
        genome=config["genome_fa"]
    output:
        idx1=f"{SUPERREADS_HISAT2_INDEX_PREFIX}.1.ht2"
    threads: max(1, int(config.get("threads_superreads_index", 4)))
    resources:
        mem_mb=int(config.get("mem_superreads_index_mb", 32000)),
        runtime=int(config.get("runtime_superreads_index_min", 360))
    conda:
        "../envs/superreads.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{SUPERREADS_HISAT2_INDEX_PREFIX}")"

        hisat2-build \
          "{input.genome}" \
          "{SUPERREADS_HISAT2_INDEX_PREFIX}"

        test -s "{output.idx1}"
        """


rule gmap_genome_index:
    input:
        genome=config["genome_fa"]
    output:
        done=f"{SUPERREADS_GMAP_DIR}/{SUPERREADS_GMAP_NAME}.done"
    threads: max(1, int(config.get("threads_superreads_index", 4)))
    resources:
        mem_mb=int(config.get("mem_superreads_index_mb", 32000)),
        runtime=int(config.get("runtime_superreads_index_min", 360))
    conda:
        "../envs/superreads.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{SUPERREADS_GMAP_DIR}"

        gmap_build \
          -D "{SUPERREADS_GMAP_DIR}" \
          -d "{SUPERREADS_GMAP_NAME}" \
          "{input.genome}"

        touch "{output.done}"
        """


rule install_superreads_module:
    output:
        done=f"{SUPERREADS_INSTALL_DIR}/tools/.installed.done",
        create_superreads=f"{SUPERREADS_INSTALL_DIR}/tools/bin/createSuperReads_RNA",
        assign_reads=f"{SUPERREADS_INSTALL_DIR}/tools/bin/assign_reads"
    params:
        repo=SUPERREADS_REPO,
        ref=SUPERREADS_REF,
        src_dir=f"{SUPERREADS_INSTALL_DIR}/src",
        install_dir=f"{SUPERREADS_INSTALL_DIR}/tools"
    threads: 1
    resources:
        mem_mb=int(config.get("mem_superreads_install_mb", 16000)),
        runtime=int(config.get("runtime_superreads_install_min", 240))
    conda:
        "../envs/superreads.yaml"
    shell:
        r"""
        set -euo pipefail
        rm -rf "{params.src_dir}" "{params.install_dir}"
        mkdir -p "{SUPERREADS_INSTALL_DIR}"

        git clone \
          --branch "{params.ref}" \
          --depth 1 \
          "{params.repo}" \
          "{params.src_dir}"

        cd "{params.src_dir}/SuperReads_RNA"
        DEST="{params.install_dir}" ./install.sh

        test -x "{output.create_superreads}"
        test -x "{output.assign_reads}"
        touch "{output.done}"
        """


rule run_superreads:
    input:
        r1=trimmed_fastq_r1,
        r2=trimmed_fastq_r2,
        install_done=rules.install_superreads_module.output.done,
        create_superreads=rules.install_superreads_module.output.create_superreads,
        assign_reads=rules.install_superreads_module.output.assign_reads,
        hisat_index=f"{SUPERREADS_HISAT2_INDEX_PREFIX}.1.ht2",
        gmap_done=f"{SUPERREADS_GMAP_DIR}/{SUPERREADS_GMAP_NAME}.done",
        script="scripts/run_superreads.py"
    output:
        bam=f"{SUPERREADS_DIR}/{{sample}}/sr_merge.bam",
        bai=f"{SUPERREADS_DIR}/{{sample}}/sr_merge.bam.bai"
    threads: max(1, int(config.get("threads_superreads", 8)))
    resources:
        mem_mb=int(config.get("mem_superreads_mb", 32000)),
        runtime=int(config.get("runtime_superreads_min", 1440))
    params:
        install_dir=f"{SUPERREADS_INSTALL_DIR}/tools",
        outdir=lambda wc: f"{SUPERREADS_DIR}/{wc.sample}",
        hisat_index=SUPERREADS_HISAT2_INDEX_PREFIX,
        gmap_dir=SUPERREADS_GMAP_DIR,
        gmap_name=SUPERREADS_GMAP_NAME,
        frag_len=float(config.get("superreads_frag_len", 100)),
        frag_std=float(config.get("superreads_frag_std", 20)),
        hisat_rna_strandness=HISAT2_RNA_STRANDNESS
    conda:
        "../envs/superreads.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"

        python "{input.script}" \
          --read1 "{input.r1}" \
          --read2 "{input.r2}" \
          --install-dir "{params.install_dir}" \
          --hisat-index "{params.hisat_index}" \
          --gmap-dir "{params.gmap_dir}" \
          --gmap-index "{params.gmap_name}" \
          --outdir "{params.outdir}" \
          --bam-out "{output.bam}" \
          --threads {threads} \
          --frag-len {params.frag_len} \
          --frag-std {params.frag_std} \
          --hisat-rna-strandness "{params.hisat_rna_strandness}"

        test -s "{output.bam}"
        test -s "{output.bai}"
        """
