SUPERREADS_REPO = str(
    config.get("superreads_repo", "https://github.com/gpertea/stringtie.git")
)

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

SUPERREADS_REPO = str(
    config.get("superreads_repo", "https://github.com/gpertea/stringtie.git")
)

rule install_superreads_module:
    output:
        done=f"{SUPERREADS_INSTALL_DIR}/.installed.done",
        create_rna_sr=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA/create_rna_sr.py",
        bin_dir_sentinel=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA/bin/createSuperReads_RNA"
    params:
        repo=SUPERREADS_REPO,
        src_dir=SUPERREADS_INSTALL_DIR,
        sr_dir=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA"
    threads: 4
    resources:
        mem_mb=int(config.get("mem_superreads_install_mb", 16000)),
        runtime=int(config.get("runtime_superreads_install_min", 240))
    conda:
        "../envs/superreads.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.src_dir}"

        if [ ! -d "{params.src_dir}/.git" ]; then
            git clone \
              "{params.repo}" \
              "{params.src_dir}"
        fi

        cd "{params.sr_dir}"
        ./install.sh

        test -f "{output.create_rna_sr}"
        test -x "{output.bin_dir_sentinel}"
        touch "{output.done}"
        """

rule run_superreads:
    input:
        r1=trimmed_fastq_r1,
        r2=trimmed_fastq_r2,
        install_done=rules.install_superreads_module.output.done,
        hisat_index=f"{SUPERREADS_HISAT2_INDEX_PREFIX}.1.ht2",
        gmap_done=f"{SUPERREADS_GMAP_DIR}/{SUPERREADS_GMAP_NAME}.done",
    output:
        bam=f"{SUPERREADS_DIR}/{{sample}}/sr_merge.bam",
        bai=f"{SUPERREADS_DIR}/{{sample}}/sr_merge.bam.bai"
    threads: max(1, int(config.get("threads_superreads", 8)))
    resources:
        mem_mb=int(config.get("mem_superreads_mb", 32000)),
        runtime=int(config.get("runtime_superreads_min", 1440))
    params:
        outdir=lambda wc: f"{SUPERREADS_DIR}/{wc.sample}",
        hisat_index=SUPERREADS_HISAT2_INDEX_PREFIX,
        gmap_dir=SUPERREADS_GMAP_DIR,
        gmap_name=SUPERREADS_GMAP_NAME,
        frag_len=float(config.get("superreads_frag_len", 100)),
        frag_std=float(config.get("superreads_frag_std", 20)),
        script=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA/create_rna_sr.py",
        sr_bin=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA/bin"
    conda:
        "../envs/superreads.yaml"
    shell:
        r"""
        set -euo pipefail
        export PATH="{params.sr_bin}:$PATH"
        mkdir -p "{params.outdir}"

        python "{params.script}" \
          -1 "{input.r1}" \
          -2 "{input.r2}" \
          -H "{params.hisat_index}" \
          -G "{params.gmap_name}" \
          -g "{params.gmap_dir}" \
          -o "{params.outdir}" \
          -p {threads} \
          --frag-len {params.frag_len} \
          --frag-std {params.frag_std}

        test -s "{output.bam}"
        test -s "{output.bai}"
        """