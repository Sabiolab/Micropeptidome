rule merge_shortstop_output:
    input:
        predict_done=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/predict.done",
        script=lambda wc: config["merge_script"]
    output:
        merged=f"{CONDITION_MERGED_DIR}/{{condition}}.merged.csv"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    params:
        root=CONDITION_RESULTS_DIR,
        outdir=CONDITION_MERGED_DIR,
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
          --samples "{wildcards.condition}" \
          --pred_csv "{params.pred_csv}" \
          $MINPROB_ARGS

        test -s "{output.merged}"
        """

rule aggregate_condition_smorfs_by_locus:
    input:
        merged_csv=f"{CONDITION_MERGED_DIR}/{{condition}}.merged.csv",
        smorf_gtf=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/{{condition}}.smorfs_shortstop.gtf",
        script=lambda wc: config["aggregate_script"]
    output:
        all_loci=f"{COHORT_PREFIX}.{{condition}}.all_loci.csv"
    threads: 1
    resources:
        mem_mb=12000,
        runtime=240
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output.all_loci}")"

        python "{input.script}" \
          --merged_csv "{input.merged_csv}" \
          --smorf_gtf "{input.smorf_gtf}" \
          --out_csv "{output.all_loci}"

        test -s "{output.all_loci}"
        """
