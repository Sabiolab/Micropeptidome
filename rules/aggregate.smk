MERGED_SHORTSTOP_CSV = f"{MERGED_DIR}/stringtie_merged.merged.csv"

rule process_shortstop_output:
    input:
        predict_done=SHORTSTOP_DONE,
        script=lambda wc: config["merge_script"]
    output:
        merged=MERGED_SHORTSTOP_CSV
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    params:
        root=OUTDIR,
        outdir=MERGED_DIR,
        min_prob=config.get("min_prob", None),
        pred_csv=config.get("pred_csv", "sams.csv")
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"

        MINPROB_ARGS=""
        if [ "{params.min_prob}" != "None" ] && [ -n "{params.min_prob}" ]; then
          MINPROB_ARGS="--min_prob {params.min_prob}"
        fi

        python "{input.script}" \
          --root "{params.root}" \
          --outdir "{params.outdir}" \
          --samples "stringtie_merged" \
          --pred_csv "{params.pred_csv}" \
          $MINPROB_ARGS

        test -s "{output.merged}"
        """

rule aggregate_smorfs_by_locus:
    input:
        merged_csv=MERGED_SHORTSTOP_CSV,
        script=lambda wc: config["aggregate_script"]
    output:
        all_loci=f"{COHORT_PREFIX}.all_loci.csv"
    threads: 1
    resources:
        mem_mb=12000,
        runtime=240
    params:
        merged_dir=MERGED_DIR,
        out_prefix=COHORT_PREFIX
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{params.out_prefix}")"

        python "{input.script}" \
          --merged_dir "{params.merged_dir}" \
          --out_prefix "{params.out_prefix}"

        test -s "{output.all_loci}"
        """
