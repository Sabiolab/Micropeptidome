SHORTSTOP_OUTDIR = f"{OUTDIR}/stringtie_merged/shortstop/shortstop_output"

rule install_shortstop:
    input:
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
        smorfs_gtf=MERGED_ANNOTATED_GTF,
        install_done=f"{OUTDIR}/.deps/shortstop_installed.done"
    output:
        done=SHORTSTOP_DONE
    threads: config.get("threads_shortstop", 8)
    params:
        min_prob=config.get("min_prob")
    resources:
        mem_mb=16000,
        runtime=int(config.get("runtime_shortstop_predict_min", 120))
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{SHORTSTOP_OUTDIR}"

        unset PYTHONPATH || true
        unset LD_PRELOAD || true
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"
        cache_root="${{SLURM_TMPDIR:-${{TMPDIR:-/tmp}}}}/shortstop_merged"
        mkdir -p "$cache_root/numba" "$cache_root/xdg"
        export NUMBA_CACHE_DIR="$cache_root/numba"
        export XDG_CACHE_HOME="$cache_root/xdg"

        command -v shortstop >/dev/null 2>&1
        python -c "import shortstop"

        min_prob_flag=""
        if [ "{params.min_prob}" != "None" ] && [ -n "{params.min_prob}" ]; then
          min_prob_flag="--min_prob {params.min_prob}"
        fi

        "$CONDA_PREFIX/bin/shortstop" predict \
          --genome "{input.genome}" \
          --putative_smorfs_gtf "{input.smorfs_gtf}" \
          --outdir "{SHORTSTOP_OUTDIR}" \
          --threads {threads} \
          $min_prob_flag

        touch "{output.done}"
        """
