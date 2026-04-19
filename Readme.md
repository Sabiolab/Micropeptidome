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

To change what ShortStop prediction to use, change `pred_csv: "sams.csv"` to `sams_secreted.csv` or `sams_intracellular.csv`. For more information, check out ShortStop documentation [here](https://github.com/brendan-miller-salk/ShortStop).

The "Annotator.py" script works better with Ensembl-style GTF annotations since those make a distinction between `five_prime_utr` and `three_prime_utr`. In Gencode annotations, there is no such distinction and both fall back to custom made `UTR_ORF` bucket. Regardless of which annotation you want to use, keep it consistent across STAR index generation, StringTie assembly/merge, and the downstream annotation steps.

The current workflow is:
`FASTQ -> Trim Galore -> STAR BAMs -> StringTie per patient -> StringTie merge per condition -> TD2 -> ShortStop -> condition-level locus summary -> RSEM`.

Trim Galore removes adapter-contaminated and low-quality sequence before alignment. The trimmed paired FASTQs are then aligned with `STAR` to produce per-patient coordinate-sorted genome BAMs. The STAR index is generated inside the workflow and uses `genome_gtf` plus the configurable `star_sjdb_overhang` setting from `config.yaml`.

StringTie first assembles transcripts per patient from the STAR BAMs, then merges those patient GTFs within each metadata-defined condition. TD2, smORF filtering, annotation, and ShortStop all run on the condition-level merged transcriptome, not on per-patient assemblies.

Conditions are defined through `SampleMetadata.tsv`, which must contain a `PatientID` column and a configurable condition column (default `Condition`). The patient IDs in that metadata file must match the sample names in `units.csv`.

RSEM quantification is done per condition on a condition-specific smORF reference built from that condition's `all_loci.csv`. Each patient in the condition is quantified against the same condition reference, and the pipeline emits both TPM-augmented summaries and expected-count matrices per condition.

Thus, the BAMs in this pipeline now have distinct purposes:
STAR BAM: genome alignment input for per-patient StringTie assembly.
Bowtie2/RSEM BAM: alignments to the condition-specific smORF transcriptome for quantification.

The pipeline also exports locus-by-sample matrices from the RSEM `expected_count` field:
`<cohort>.<condition>.all_loci.blastp_human.expected_counts.tsv`

These are count-like RSEM expected counts for the loci present in the final BLAST-annotated summaries, not TPM values and not integer raw read counts.


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
