rule install_shortstop:
    input:
        # Force reinstall if the environment spec changes (Snakemake hashes envs, but this
        # marker file would otherwise prevent re-running in the new env).
        env_spec="envs/shortstop.yaml"
    output:
        done=f"{OUTDIR}/.deps/shortstop_installed.done"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30
    conda:
        "../envs/shortstop.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/.deps"
        python -m pip install --no-deps git+https://github.com/brendan-miller-salk/ShortStop.git
        touch "{output.done}"
        """

rule shortstop_predict:
    input:
        shortstop_installed=rules.install_shortstop.output.done,
        genome=config["genome_fa"],
        smorfs_gtf=cohort_shortstop_gtf()
    output:
        done=cohort_shortstop_done()
    threads: config.get("threads_shortstop", 8)
    resources:
        mem_mb=16000,
        runtime=int(config.get("runtime_shortstop_predict_min", 120))
    conda:
        "../envs/shortstop.yaml"
    shell:
        r"""
        set -euo pipefail

        cd "{cohort_results_dir()}/shortstop"
        mkdir -p shortstop_output

        unset PYTHONPATH || true
        unset LD_PRELOAD || true
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"
        cache_root="${{SLURM_TMPDIR:-${{TMPDIR:-/tmp}}}}/shortstop_{COHORT_LABEL}"
        mkdir -p "$cache_root/numba" "$cache_root/xdg"
        export NUMBA_CACHE_DIR="$cache_root/numba"
        export XDG_CACHE_HOME="$cache_root/xdg"

        test -x "$CONDA_PREFIX/bin/shortstop" || {{
          echo "ShortStop CLI not found in $CONDA_PREFIX/bin" >&2
          python -m pip show ShortStop >&2 || true
          exit 1
        }}
        python -c "import shortstop; print(shortstop.__file__)"

        "$CONDA_PREFIX/bin/shortstop" predict \
          --genome "{input.genome}" \
          --putative_smorfs_gtf "{input.smorfs_gtf}" \
          --outdir shortstop_output \
          --threads {threads}

        touch "{output.done}"
        """
