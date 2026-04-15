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
        done=f"{SUPERREADS_INSTALL_DIR}/.installed.done",
        create_rna_sr=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA/create_rna_sr.py",
        bin_dir_sentinel=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA/bin/createSuperReads_RNA"
    params:
        repo=SUPERREADS_REPO,
        ref=SUPERREADS_REF,
        src_dir=SUPERREADS_INSTALL_DIR,
        sr_dir=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA",
        install_log=f"{SUPERREADS_INSTALL_DIR}/SuperReads_RNA/install.log",
        patch_script=str(Path(workflow.basedir) / "scripts" / "patch_superreads_source.py")
    threads: 4
    resources:
        mem_mb=int(config.get("mem_superreads_install_mb", 16000)),
        runtime=int(config.get("runtime_superreads_install_min", 240))
    conda:
        "../envs/superreads.yaml"
    shell:
        r"""
        set -euo pipefail

        if [ -d "{params.src_dir}/.git" ]; then
            git -C "{params.src_dir}" fetch --tags --force
            git -C "{params.src_dir}" checkout "{params.ref}"
            git -C "{params.src_dir}" reset --hard "{params.ref}"
            git -C "{params.src_dir}" clean -fdx
        elif [ -f "{params.sr_dir}/install.sh" ]; then
            echo "Using existing SuperReads_RNA source tree at {params.sr_dir}"
        else
            if [ -d "{params.src_dir}" ] && [ -n "$(find "{params.src_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
                partial_dir="{params.src_dir}.partial.$(date +%Y%m%d%H%M%S).$$"
                echo "Moving incomplete SuperReads dependency directory to $partial_dir" >&2
                mv "{params.src_dir}" "$partial_dir"
            fi
            git clone \
              "{params.repo}" \
              "{params.src_dir}"
            git -C "{params.src_dir}" checkout "{params.ref}"
        fi

        python "{params.patch_script}" "{params.sr_dir}"

        if [ -n "${{CONDA_PREFIX:-}}" ]; then
            export PATH="${{CONDA_PREFIX}}/bin:$PATH"
            if [ -x "${{CONDA_PREFIX}}/bin/m4" ]; then
                export M4="${{CONDA_PREFIX}}/bin/m4"
            fi
        fi

        m4_path="$(command -v m4 || true)"
        m4_version="$(m4 --version 2>/dev/null | head -n 1 || true)"
        if ! printf "%s\n" "$m4_version" | grep -qi "GNU M4"; then
            echo "ERROR: SuperReads build requires GNU m4 >= 1.4 from the Snakemake conda env." >&2
            echo "Resolved m4: ${{m4_path:-not found}}" >&2
            echo "m4 version: ${{m4_version:-unavailable}}" >&2
            echo "CONDA_PREFIX: ${{CONDA_PREFIX:-unset}}" >&2
            echo "Remove the stale SuperReads conda env and rerun so Snakemake rebuilds envs/superreads.yaml." >&2
            exit 1
        fi

        cd "{params.sr_dir}"
        if ! ./install.sh > "{params.install_log}" 2>&1; then
            echo "ERROR: SuperReads install.sh failed. Full log: {params.install_log}" >&2
            echo "----- install.log tail -----" >&2
            tail -n 200 "{params.install_log}" >&2 || true
            echo "----------------------------" >&2
            exit 1
        fi

        if [ ! -f "{output.create_rna_sr}" ]; then
            echo "ERROR: install.sh completed but did not create {output.create_rna_sr}" >&2
            echo "Full log: {params.install_log}" >&2
            exit 1
        fi
        if [ ! -x "{output.bin_dir_sentinel}" ]; then
            echo "ERROR: install.sh completed but did not create executable {output.bin_dir_sentinel}" >&2
            echo "Full log: {params.install_log}" >&2
            exit 1
        fi
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
