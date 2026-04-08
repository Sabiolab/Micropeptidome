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

The "Annotator.py" script works better with Ensembl-style GTF annotations since those make a distinction between `five_prime_utr` and `three_prime_utr`. In Gencode annotations, there is no such distinction and both fall back to custom made `UTR_ORF` bucket. Regardless of which annotation you want to use, keep it consistent across TopHat2, SuperReads indexing, and the downstream annotation steps.

The current workflow is:
`FASTQ -> Trim Galore -> TopHat2 BAMs + StringTie SuperReads BAMs -> StringTie -> TD2 -> ShortStop -> locus aggregation -> RSEM`.

Trim Galore removes adapter-contaminated and low-quality sequence before any alignment. The trimmed paired FASTQs are then used in two separate genome-alignment branches:

1. `TopHat2` generates per-sample genome-aligned BAM files (`accepted_hits.bam`).
2. The optional `StringTie` SuperReads helper is installed from the official `gpertea/stringtie` repository, builds super-reads from the same trimmed FASTQs, and produces a merged BAM that can be fed directly into StringTie.

StringTie assembles transcripts from whichever BAM source is selected by `stringtie_bam_source` in `config.yaml` (`superreads` by default, `tophat` as a fallback). The pipeline still pins StringTie to `3.0.3` and exposes `stringtie_rRNA`; set `stringtie_rRNA: true` to pass `-N` for Total RNA / rRNA-depleted libraries, and leave it `false` for polyA-selected libraries.

RSEM quantification is still done on a separate reference: the custom smORF transcriptome built by `rsem-prepare-reference --bowtie2`. The only change is that RSEM now consumes the Trim Galore outputs instead of the raw FASTQs.

Thus, the BAMs in this pipeline now have distinct purposes:
TopHat2 BAM: splice-aware genome alignments, retained as explicit per-sample alignment outputs.
SuperReads BAM: merged short-read + super-read genome alignments for StringTie assembly.
Bowtie2/RSEM BAM: alignments to the custom smORF transcriptome for quantification.

The pipeline also exports locus-by-sample matrices from the RSEM `expected_count` field:
`<cohort>.all_loci.blastp_human.expected_counts.tsv`
`<cohort>.shared_ge<N>.blastp_human.expected_counts.tsv`

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

If the failure happens while creating `envs/superreads.yaml`, note that the SuperReads helper is built from the official StringTie GitHub repository during the workflow. That step requires outbound GitHub access plus a compiler/autotools toolchain inside the environment.

If `install_superreads_module` fails with an autotools error like `autom4te: need GNU m4 1.4 or later: /usr/bin/m4`, remove the failed SuperReads conda environment and rerun that rule so Snakemake rebuilds it with GNU `m4` available inside the env.

2. SLURM execution may be expressed as either `--slurm` or `--executor slurm` depending on snakemake version. If the first one does not work for you, try the second one.
