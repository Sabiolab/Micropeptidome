# Setup

This pipeline requires Snakemake to run. It is recommended to set up a conda enviroment with snakemake. To install a snakemake in a conda enviroment, use:
```{}
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
snakemake --version
```
To learn more about conda, please visit [Anaconda](https://anaconda.org/channels/anaconda/packages/conda/overview).

To download this pipeline, use:
```{}
git clone https://github.com/SantiBarber/Micropeptidome.git
cd Micropeptidome
```

Then, please download the required proteome/genome references. Please find below a few examples. I would recommend downloading the Ensembl GTF annotation and FASTA since those make a distinction between 5' and 3' UTRs that Genecode annotation does not make. However, both should work fine.

This pipeline can be applied to both human and mouse data. Keep in mind that the original ShortStop was trained with human data. For more information about this, please refere to the original ShortStop repository [here](https://github.com/brendan-miller-salk/ShortStop) and ShortStop paper [here](https://link.springer.com/article/10.1186/s44330-025-00037-4): 

### Genecode

#### Genome.fa FASTA (GRCh38 primary assembly)
```{}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```
#### Annotation GTF (GRCh38; contigs like "chr1", "chr2"... to match the FASTA above)
```{}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
gunzip gencode.v38.primary_assembly.annotation.gtf.gz
```
### Ensembl

#### Genome FASTA (GRCh38 primary assembly, unmasked)
```{}
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
#### Annotation GTF (GRCh38; contigs like "1", "2", ... to match the FASTA above)
```{}
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz
```

### Proteome

#### Proteome.faa (FASTA)
```{}
wget -O human_proteome.faa "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=fasta"
```

# Run μ-Peptidome analysis

## Generate the Conda environmnets

The firts time we set up the pipeline, we need to create the conda enviroments containing the software used by the pipeline. To do this, use the following command below and change the path `/path/to/conda_envs` to the path where you want the enviroments to be created. The command below will create the enviroments for you.
```{}
 snakemake --use-conda --conda-create-envs-only \
  --conda-frontend conda \
  --conda-prefix /path/to/conda_envs \
  -j 1 --latency-wait 60
```

## Run the pipeline

Once the enviroments have been set up, change the settings in the `config.yaml` file and make sure all the variables are set, including the full pathname to the GTF and FASTA annotations mentioned above. Also, remmber to change `/path/to/conda_envs` to the directory were you created the enviroments.

Them, you can run the pipeline with:
```{}
snakemake --use-conda --slurm -j 32 --conda-frontend conda --conda-prefix /path/to/conda_envs \
  --rerun-incomplete \
  --latency-wait 60
```
The last flags are not strictly necessary.

# Further considerations

To use a prebuilt STAR index, set in config.yaml: `star_index_dir: "/path/to/existing/star_index/2.7.10a"` and ensure it exists. The default index is set as 2.7.10a, but this can be changed in the `STAR.yaml` to any other version.

To change what ShortStop prediction to use, change `pred_csv: "sams.csv"` to `sams_secreted.csv` or `sams_intracellular.csv`. For more information, check out ShortStop documentation [here](https://github.com/brendan-miller-salk/ShortStop).

The "Annotator.py" script works better with Ensembl-style GTF annotations since those make a distinction between `five_prime_utr` and `three_prime_utr`. In Gencode annotations, there is no such distinction and both fall back to custom made `UTR_ORF` bucket. Regardless of which annotation you want to use, keep it consistent (specially if you are using that annotation for `STAR` alignment).

StringTie takes the STAR-aligned BAM generated from FASTQs and uses it for transcript assembly using a GTF reference (which can be the same reference mentioned above).

RSEM quant is done on a different reference (the custom smORF transcriptome built by `rsem-prepare-reference --bowtie2`), so the pipeline alignes the FASTQs again with Bowtie2 to that smORF reference and feed the BAM into `rsem-calculate-expression --alignments`. We use bowtie2 because it is lighter for this task, it is built percisely for transcriptome alignment (whereas STAR has a genome-first mentality with splice awarenes that is not necesarily useful here) and STAR multi-mapping can be troublesom for short sequences.

Thus, those two BAMs are fundamentally different:
STAR BAM: splice-aware alignments to the genome (for StringTie).
Bowtie2 BAM: alignments to the smORF transcriptome reference (for RSEM quantification on smORFs).


## Troubleshooting

1. Errors building the enviroments

Nuke the environment directory
```{}
rm -rf /home/sbarber/conda_envs
```

Clean packages and tarballs and try again
```{}
conda clean --packages --tarballs -y
```

Then, try again!

2. SLURM execution may be expressed as either `--slurm` or `--executor slurm` depending on snakemake version. If the first one does not work for you, try the second one.
