rule merge_shortstop_output:
    input:
        predict_done=cohort_shortstop_done(),
        script=config["merge_script"]
    output:
        merged=cohort_merged_csv()
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    params:
        root=COHORT_RESULTS_ROOT,
        outdir=COHORT_MERGED_DIR,
        min_prob=config.get("min_prob", None),
        pred_csv=config.get("pred_csv", "sams.csv"),
        cohort_label=COHORT_LABEL
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
          --samples "{params.cohort_label}" \
          --pred_csv "{params.pred_csv}" \
          $MINPROB_ARGS

        test -s "{output.merged}"
        """

rule aggregate_cohort_smorfs_by_locus:
    input:
        merged_csv=cohort_merged_csv(),
        smorf_gtf=cohort_shortstop_gtf(),
        script=config["aggregate_script"]
    output:
        all_loci=temp(cohort_all_loci())
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
