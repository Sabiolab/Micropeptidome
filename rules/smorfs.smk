MERGED_SMORFS_FA   = f"{MERGED_TX_DIR}/merged.smorfs.fa"
MERGED_SMORFS_IDS  = f"{MERGED_TX_DIR}/merged.smorf_ids.txt"
MERGED_SMORFS_GFF3 = f"{MERGED_TX_DIR}/merged.smorfs.gff3"
MERGED_SMORFS_RAW_GTF = f"{MERGED_TX_DIR}/merged.smorfs_shortstop.raw.gtf"

rule filter_smorfs:
    input:
        pep=MERGED_TD2_PEP,
        gff3=MERGED_TD2_GFF,
        script=config["filter_smorf_pep_py"]
    output:
        smorfs_fa=MERGED_SMORFS_FA,
        ids=MERGED_SMORFS_IDS,
        smorfs_gff3=MERGED_SMORFS_GFF3
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        python "{input.script}" \
          "{input.pep}" \
          --min_len {config[min_aa]} \
          --max_len {config[max_aa]} \
          --out_fasta "{output.smorfs_fa}" \
          --out_ids "{output.ids}"

        if [ ! -s "{output.ids}" ]; then
          echo "No smORFs remained after filtering merged TD2 peptides to the {config[min_aa]}-{config[max_aa]} aa range." >&2
          exit 1
        fi

        # Exact attribute matching to avoid grep substring collisions
        awk 'FNR==NR {{ ids[$0]=1; next }}
             /^#/ {{ print; next }}
             {{
               n=split($9,a,";");
               for(i=1;i<=n;i++) {{
                 m=split(a[i],kv,"=");
                 if(m==2 && (kv[1]=="ID" || kv[1]=="Parent") && kv[2] in ids) {{
                   print; break
                 }}
               }}
             }}' "{output.ids}" "{input.gff3}" > "{output.smorfs_gff3}"

        test -s "{output.smorfs_fa}"
        test -s "{output.ids}"
        test -s "{output.smorfs_gff3}"
        """

rule tx_to_genome_gtf:
    input:
        merged_gtf=MERGED_GTF,
        smorfs_gff3=MERGED_SMORFS_GFF3,
        script=config["tx_to_genome_py"]
    output:
        gtf=MERGED_SMORFS_RAW_GTF
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/smORFs.yaml"
    shell:
        r"""
        set -euo pipefail

        python "{input.script}" \
          --merged_gtf "{input.merged_gtf}" \
          --smorfs_gff3 "{input.smorfs_gff3}" \
          --out_gtf "{output.gtf}"

        test -s "{output.gtf}"
        """
