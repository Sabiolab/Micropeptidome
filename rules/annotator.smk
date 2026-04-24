rule annotator_smorf_types:
    input:
        smorf_gtf=f"{cohort_results_dir()}/shortstop/{COHORT_LABEL}.smorfs_shortstop.raw.gtf",
        genome_gtf=config["genome_gtf"]
    output:
        intersect=f"{cohort_results_dir()}/shortstop/lineintersect.gtf",
        non_intersect=f"{cohort_results_dir()}/shortstop/linenonintersect.gtf",
        annotations=f"{cohort_results_dir()}/shortstop/Annotations.txt"
    threads: config.get("threads_annotator", 1)
    resources:
        mem_mb=16000,
        runtime=120
    params:
        shortstop_dir=f"{cohort_results_dir()}/shortstop"
    conda:
        "../envs/BedTools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.shortstop_dir}"

        python "scripts/Annotator/Annotator.py" smorf_types \
          --smorf_gtf "{input.smorf_gtf}" \
          --ensembl_gtf "{input.genome_gtf}" \
          --outdir "{params.shortstop_dir}" \
          --intersect_output "{output.intersect}" \
          --non_intersect_output "{output.non_intersect}" \
          --output_file "{output.annotations}" \
          --threads {threads}
        """

rule annotate_smorfs_gtf:
    input:
        smorf_gtf=f"{cohort_results_dir()}/shortstop/{COHORT_LABEL}.smorfs_shortstop.raw.gtf",
        annotations=f"{cohort_results_dir()}/shortstop/Annotations.txt"
    output:
        annotated_gtf=cohort_shortstop_gtf()
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        python "scripts/add_smorf_type_to_gtf.py" \
          --gtf "{input.smorf_gtf}" \
          --annotations "{input.annotations}" \
          --out "{output.annotated_gtf}"
        """
