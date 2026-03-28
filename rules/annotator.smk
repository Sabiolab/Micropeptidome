rule annotator_smorf_types:
    input:
        smorf_gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/{{sample}}.smorfs_shortstop.raw.gtf",
        genome_gtf=config["genome_gtf"]
    output:
        intersect=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/lineintersect.gtf",
        non_intersect=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/linenonintersect.gtf",
        annotations=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/Annotations.txt"
    threads: config.get("threads_annotator", 1)
    resources:
        mem_mb=16000,
        runtime=120
    conda:
        "../envs/BedTools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/shortstop"

        python "scripts/Annotator/Annotator.py" smorf_types \
          --smorf_gtf "{input.smorf_gtf}" \
          --ensembl_gtf "{input.genome_gtf}" \
          --outdir "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/shortstop" \
          --intersect_output "{output.intersect}" \
          --non_intersect_output "{output.non_intersect}" \
          --output_file "{output.annotations}" \
          --threads {threads}
        """

rule annotate_smorfs_gtf:
    input:
        smorf_gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/{{sample}}.smorfs_shortstop.raw.gtf",
        annotations=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/Annotations.txt"
    output:
        annotated_gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/{{sample}}.smorfs_shortstop.gtf"
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
