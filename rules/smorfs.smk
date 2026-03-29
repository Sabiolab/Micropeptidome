rule filter_smorfs:
    input:
        pep=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa.TD2.pep",
        gff3=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/transcripts/{{sample}}.transcripts.fa.TD2.gff3",
        script=config["filter_smorf_pep_py"]
    output:
        smorfs_fa=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/smorfs/{{sample}}.smorfs.fa",
        ids=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/smorfs/{{sample}}.smorf_ids.txt",
        smorfs_gff3=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/smorfs/{{sample}}.smorfs.gff3"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/smorfs"

        python "{input.script}" \
          "{input.pep}" \
          --min_len {config[min_aa]} \
          --max_len {config[max_aa]} \
          --out_fasta "{output.smorfs_fa}" \
          --out_ids "{output.ids}"

        if [ ! -s "{output.ids}" ]; then
          echo "No smORFs remained for sample {wildcards.sample} after filtering TD2 peptides to the {config[min_aa]}-{config[max_aa]} aa range." >&2
          exit 1
        fi

        # This line can result in smaller than expected smORFs in the end
        grep -F -f "{output.ids}" "{input.gff3}" > "{output.smorfs_gff3}"

        test -s "{output.smorfs_fa}"
        test -s "{output.ids}"
        test -s "{output.smorfs_gff3}"
        """

rule tx_to_genome_gtf:
    input:
        sample_gtf=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/stringtie/{{sample}}.gtf",
        smorfs_gff3=f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/smorfs/{{sample}}.smorfs.gff3",
        script=config["tx_to_genome_py"]
    output:
        gtf=temp(f"{RESULTS_SHORTSTOP_DIR}/{{sample}}/shortstop/{{sample}}.smorfs_shortstop.raw.gtf")
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{RESULTS_SHORTSTOP_DIR}/{wildcards.sample}/shortstop"

        python "{input.script}" \
          --merged_gtf "{input.sample_gtf}" \
          --smorfs_gff3 "{input.smorfs_gff3}" \
          --out_gtf "{output.gtf}"

        test -s "{output.gtf}"
        """
