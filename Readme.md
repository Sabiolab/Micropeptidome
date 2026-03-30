# Setup

This pipeline requires Snakemake to run. It is recommended to set up a conda enviroment with snakemake. To install a snakemake in a conda enviroment, use:
```
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
snakemake --version
```
To learn more about conda, please visit [Anaconda](https://anaconda.org/channels/anaconda/packages/conda/overview).

To download this pipeline, use:
```
git clone https://github.com/SantiBarber/Micropeptidome.git
```

The work in progress version:
```
git clone --branch Working_TD2 https://github.com/SantiBarber/Micropeptidome.git
```

Please download the required proteome/genome references. I would recommend downloading the Ensembl GTF annotation and FASTA since it makes a distinction between 5' and 3' UTRs that Genecode annotation does not. However, Ensembl uses contigs such as `1, 2...` for chromosomes which messes up ShortStop (`str` vs `int`). I have to corret this in the ShortStop at some point. FOR NOW, use the Genecode annotation.

The original `ShortStop` published by Dr. Brendan Miller was trained with human data. For more information about `ShorStop`, please refere to the original repository [here](https://github.com/brendan-miller-salk/ShortStop) and paper [here](https://link.springer.com/article/10.1186/s44330-025-00037-4).

### Genecode - Recommended for now

#### Genome.fa FASTA (GRCh38 primary assembly)
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
```
#### Annotation GTF (GRCh38; contigs like "chr1", "chr2"... to match the FASTA above)
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
gunzip gencode.v38.primary_assembly.annotation.gtf.gz
```
### Ensembl

#### Genome FASTA (GRCh38 primary assembly, unmasked)
```
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
#### Annotation GTF (GRCh38; contigs like "1", "2", ... to match the FASTA above)
```
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz
```

### Proteome

#### Proteome.faa (FASTA)
```
wget -O human_proteome.faa "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&format=fasta"
```

# Run μ-Peptidome analysis

## Run the pipeline

The firts the pipeline runs, it builds the conda enviroments containing dependencies. Change the path `/path/to/conda_envs` to the path where you want the enviroments to be created. 

You can change the settings in the `config.yaml` file. Importantly, set the full pathname to the GTF and FASTA annotations.

You can create the enviroments and run the pipeline with:
```
snakemake --use-conda --slurm -j 32 --conda-frontend conda --conda-prefix /path/to/conda_envs \
  --rerun-incomplete \
  --latency-wait 60
```
The last flags are not strictly necessary.

## Further considerations

To use a prebuilt STAR index, set in config.yaml: `star_index_dir: "/path/to/existing/star_index/2.7.10a"` and ensure it exists. The default index is set as 2.7.10a, but this can be changed in the `STAR.yaml` to any other version.

To change what ShortStop prediction to use, change `pred_csv: "sams.csv"` to `sams_secreted.csv` or `sams_intracellular.csv`. For more information, check out ShortStop documentation [here](https://github.com/brendan-miller-salk/ShortStop).

The "Annotator.py" script works better with Ensembl-style GTF annotations since those make a distinction between `five_prime_utr` and `three_prime_utr`. In Gencode annotations, there is no such distinction and both fall back to custom made `UTR_ORF` bucket. Regardless of which annotation you want to use, keep it consistent (specially if you are using that annotation for `STAR` alignment).

StringTie takes the STAR-aligned BAM generated from FASTQs and uses it for transcript assembly using the GTF as a reference, but it also reconstructs transcripts that are not present in the reference when there is enough transcriptomic evidence. The pipeline pins StringTie to `3.0.3` and exposes `stringtie_rRNA` in `config.yaml`; set `stringtie_rRNA: true` to pass `-N` for Total RNA / rRNA-depleted libraries, and leave it `false` for polyA-selected libraries.

RSEM quant is done on a different reference (the custom smORF transcriptome built by `rsem-prepare-reference --bowtie2`), so the pipeline alignes the FASTQs again with Bowtie2 to that smORF reference and feed the BAM into `rsem-calculate-expression --alignments`. We use bowtie2 because it is lighter for this task, it is built percisely for transcriptome alignment (whereas STAR has a genome-first mentality with splice awarenes that is not necesarily useful here) and STAR multi-mapping can be troublesom for short sequences.

Thus, those two BAMs are fundamentally different:
STAR BAM: splice-aware alignments to the genome (for StringTie).
Bowtie2 BAM: alignments to the smORF transcriptome reference (for RSEM quantification on smORFs).


## Troubleshooting

1. Errors building the enviroments

Nuke the environment directory
```
rm -rf /path/to/your/conda_envs
```

Clean packages and tarballs and try again
```
conda clean --packages --tarballs -y
```

Then, try again!

If the failure happens while creating `envs/smORFs.yaml` or `envs/shortstop.yaml`, remove the failed environment directory and let Snakemake recreate it. `ShortStop` is installed with `pip --no-deps`, and its Python/ML stack now lives in `envs/shortstop.yaml` so it does not conflict with the newer `StringTie` toolchain in `envs/smORFs.yaml`.

2. SLURM execution may be expressed as either `--slurm` or `--executor slurm` depending on snakemake version. If the first one does not work for you, try the second one.
