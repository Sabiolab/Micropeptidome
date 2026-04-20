MERGED_ANNOTATOR_DIR     = f"{OUTDIR}/stringtie_merged/shortstop"
MERGED_INTERSECT_GTF     = f"{MERGED_ANNOTATOR_DIR}/lineintersect.gtf"
MERGED_NON_INTERSECT_GTF = f"{MERGED_ANNOTATOR_DIR}/linenonintersect.gtf"
MERGED_ANNOTATIONS_TXT   = f"{MERGED_ANNOTATOR_DIR}/Annotations.txt"
MERGED_ANNOTATED_GTF     = f"{MERGED_ANNOTATOR_DIR}/stringtie_merged.smorfs_shortstop.gtf"

rule annotator_smorf_types:
    input:
        smorf_gtf=MERGED_SMORFS_RAW_GTF,
        genome_gtf=config["genome_gtf"]
    output:
        intersect=MERGED_INTERSECT_GTF,
        non_intersect=MERGED_NON_INTERSECT_GTF,
        annotations=MERGED_ANNOTATIONS_TXT
    threads: config.get("threads_annotator", 1)
    resources:
        mem_mb=16000,
        runtime=120
    conda:
        "../envs/BedTools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{MERGED_ANNOTATOR_DIR}"

        python "scripts/Annotator/Annotator.py" smorf_types \
          --smorf_gtf "{input.smorf_gtf}" \
          --ensembl_gtf "{input.genome_gtf}" \
          --outdir "{MERGED_ANNOTATOR_DIR}" \
          --intersect_output "{output.intersect}" \
          --non_intersect_output "{output.non_intersect}" \
          --output_file "{output.annotations}" \
          --threads {threads}

        test -s "{output.annotations}"
        test -s "{output.intersect}"
        """

rule annotate_smorfs_gtf:
    input:
        smorf_gtf=MERGED_SMORFS_RAW_GTF,
        annotations=MERGED_ANNOTATIONS_TXT
    output:
        annotated_gtf=MERGED_ANNOTATED_GTF
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

        test -s "{output.annotated_gtf}"
        """
