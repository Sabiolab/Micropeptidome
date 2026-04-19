rule annotator_smorf_types:
    input:
        smorf_gtf=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/{{condition}}.smorfs_shortstop.raw.gtf",
        genome_gtf=config["genome_gtf"]
    output:
        intersect=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/lineintersect.gtf",
        non_intersect=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/linenonintersect.gtf",
        annotations=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/Annotations.txt"
    threads: config.get("threads_annotator", 1)
    resources:
        mem_mb=16000,
        runtime=120
    conda:
        "../envs/BedTools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{CONDITION_RESULTS_DIR}/{wildcards.condition}/shortstop"

        python "scripts/Annotator/Annotator.py" smorf_types \
          --smorf_gtf "{input.smorf_gtf}" \
          --ensembl_gtf "{input.genome_gtf}" \
          --outdir "{CONDITION_RESULTS_DIR}/{wildcards.condition}/shortstop" \
          --intersect_output "{output.intersect}" \
          --non_intersect_output "{output.non_intersect}" \
          --output_file "{output.annotations}" \
          --threads {threads}
        """

rule annotate_smorfs_gtf:
    input:
        smorf_gtf=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/{{condition}}.smorfs_shortstop.raw.gtf",
        annotations=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/Annotations.txt"
    output:
        annotated_gtf=f"{CONDITION_RESULTS_DIR}/{{condition}}/shortstop/{{condition}}.smorfs_shortstop.gtf"
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
