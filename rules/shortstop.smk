rule install_shortstop:
    input:
        # Force reinstall if the environment spec changes (Snakemake hashes envs, but this
        # marker file would otherwise prevent re-running in the new env).
        env_spec="envs/smORFs.yaml"
    output:
        done=f"{OUTDIR}/.deps/shortstop_installed.done"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/.deps"
        python -m pip install --no-deps git+https://github.com/brendan-miller-salk/ShortStop.git
        touch "{output.done}"
        """

rule shortstop_predict:
    input:
        genome=config["genome_fa"],
        smorfs_gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/{{sample}}.smorfs_shortstop.gtf",
        shortstop_done=f"{OUTDIR}/.deps/shortstop_installed.done"
    output:
        done=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/predict.done"
    threads: config.get("threads_shortstop", 8)
    resources:
        mem_mb=16000,
        runtime=120
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        cd "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/shortstop"
        mkdir -p shortstop_output

        unset PYTHONPATH || true
        unset LD_PRELOAD || true
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"

        shortstop predict \
          --genome "{input.genome}" \
          --putative_smorfs_gtf "{input.smorfs_gtf}" \
          --outdir shortstop_output \
          --threads {threads}

        touch "{output.done}"
        """