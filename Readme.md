We need to download a few files that we would be using when running the pipeline. I would recommend downloading the Ensembl GTF annotation and FASTA because, for smORF annotation later, those make a distinction between 5' and 3' UTRs. However, both should work fine.

This pipeline can be applied to both human and mouse data, but the original ShortStop was trained with human data. For more information, please refere to the original ShortStop repo: 

## Genecode

#### Genome.fa FASTA (GRCh38 primary assembly)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

#### Annotation GTF (GRCh38; contigs like "chr1", "chr2"... to match the FASTA above)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
gunzip gencode.v38.primary_assembly.annotation.gtf.gz

## Ensembl

#### Genome FASTA (GRCh38 primary assembly, unmasked)
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#### Annotation GTF (GRCh38; contigs like "1", "2", ... to match the FASTA above)
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz

Then, we should download the proteome FASTA. 

## Proteome

#### Proteome.faa (FASTA)
wget -O human_proteome.faa \
"https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=fasta"

## Generate the Conda environmnets

Create the enviroments the first time:
```{}
 snakemake --use-conda --conda-create-envs-only \
  --conda-frontend conda \
  --conda-prefix /home/sbarber/conda_envs \
  -j 1 --latency-wait 60
```

## Run the pipeline with

Run with:
```{}
snakemake --use-conda --slurm -j 32 \
  --conda-frontend conda \
  --conda-prefix /home/sbarber/conda_envs \
  --rerun-incomplete \
  --latency-wait 60
```

## Further considerations

To use a prebuilt STAR index, set in config.yaml: `star_index_dir: "/path/to/existing/star_index/2.7.10a"` and ensure /path/to/existing/star_index/Genome exists. The default index is set as 2.7.10a, but this can be changed in the `STAR.yaml` to any other version.

To change what ShortStop prediction to use, change `pred_csv: "sams.csv"` to `sams_secreted.csv` or `sams_intracellular.csv`. For more information, check out ShortStop documentation here.

The "Annotator.py" script works better with Ensembl-style GTF annotations, since those make a distinction between `five_prime_utr` and `three_prime_utr`. In Gencode annotations, there is no distinction and both fall back to custom `UTR_ORF` bucket. Regardless of which annotation you want to use, keep it consistent (specially if you are using that annotation for `STAR` alignment)

StringTie takes the STAR-aligned BAM generated from FASTQs and uses it for transcript assembly.

RSEM quant is done on a different reference (the custom smORF transcriptome built by `rsem-prepare-reference --bowtie2`), so the pipeline alignes the FASTQs again with Bowtie2 to that smORF reference and feed the BAM into `rsem-calculate-expression --alignments`.

Those two BAMs are fundamentally different:
STAR BAM: splice-aware alignments to the genome (needed for StringTie).
Bowtie2 BAM: alignments to the smORF transcriptome reference (needed for RSEM quantification on smORFs)


## Troubleshooting

1. Errors building the enviroments

Nuke the enviroment directory
```{}
rm -rf /home/sbarber/conda_envs
```

Clean packages and tarballs and try again
```{}
conda clean --packages --tarballs -y
```

Then, try again!

2. Add_new_error
