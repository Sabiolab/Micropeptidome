### SETUP
## Set up with mamba
conda install -n base -c conda-forge mamba -y

## Set up environment
mamba create -n smORFs \
  -c conda-forge -c bioconda \
  python=3.10 \
  fastqc fastp star stringtie gffread samtools transdecoder \
  -y

# for snakemake workflow
mamba install -c conda-forge -c bioconda snakemake -n smORFs -y

## install ShortStop
pip install git+https://github.com/brendan-miller-salk/ShortStop.git

## genome.fa
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

## Annotation gtf
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz
gunzip gencode.v38.primary_assembly.annotation.gtf.gz

## For ShortStop, a positive genome annotation ## DO NOT DOWNLOAD THIS WAY IF USING THE DEMO
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.basic.annotation.gtf.gz
# gunzip gencode.v43.primary_assembly.basic.annotation.gtf.gz

      # Make a copy of the GTF without any lines starting with "#"
      # grep -v '^#' gencode.v43.primary_assembly.basic.annotation.gtf \
      # > gencode.v43.primary_assembly.basic.annotation.noheader.gtf

# ## QC on raw reads (FastQC)
# fastqc sample1_R1.fastq.gz sample1_R2.fastq.gz -o qc_raw/

# ## trimming
# mkdir -p trimmed qc_fastp

# fastp \
#   -i sample1_R1.fastq.gz \
#   -I sample1_R2.fastq.gz \
#   -o trimmed/OB1_R1.trimmed.fq.gz \
#   -O trimmed/OB1_R2.trimmed.fq.gz \
#   --detect_adapter_for_pe \
#   --thread 8 \
#   --html qc_fastp/OB1_fastp.html \
#   --json qc_fastp/OB1_fastp.json

# mv star/sample1_Aligned.sortedByCoord.out.bam star/sample1.bam ## rename

# ## Build STAR genome index rose use a pre-built one
# mkdir -p star_index

# STAR --runThreadN 16 \
#   --runMode genomeGenerate \
#   --genomeDir star_index \
#   --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
#   --sjdbGTFfile gencode.v38.primary_assembly.annotation.gtf \
#   --sjdbOverhang 100

# ## Align trimmed RNA-seq reads (1 per sample)
# mkdir -p star

# STAR --runThreadN 16 \
#   --genomeDir star_index \
#   --readFilesIn trimmed/OB1_R1.trimmed.fq.gz trimmed/OB1_R2.trimmed.fq.gz \
#   --readFilesCommand zcat \
#   --outSAMtype BAM SortedByCoordinate \
#   --outFileNamePrefix star/OB1_
#   ## --outSAMunmapped Within

                                              #### WE CAN START HERE ####

# 1. Use string tie in an already aligned file (.bam) to Assemble transcripts with StringTie (which includes novel transcripts)
# stringtie will produce transript annotations like: chr1 | StringTie | exon ... gene_id | "STRG.1"; transcript_id "STRG.1.1";
# i believe this format is already ShortStop compatible? since it has genomic locations and chromosome location. However, I
# run TransDecoder later to scan for ORFs.

stringtie /storage/scratch01/groups/md/rnaseq_liver_OBECAN/results/mapped/star/OB3/Aligned.sortedByCoord.out.bam \
  -G gencode.v38.primary_assembly.annotation.gtf \
  -o OB3.gtf \
  -p 8

stringtie /storage/scratch01/groups/md/rnaseq_liver_OBECAN/results/mapped/star/OB315/Aligned.sortedByCoord.out.bam \
  -G gencode.v38.primary_assembly.annotation.gtf \
  -o OB315.gtf \
  -p 8

# This .gtf files now contains both annotated and novel transcripts found in each sample.
# I can then run TransDecoder on the files individually or merge them together (the first approach seems to work better)

# If merging files together, use the code below. Else skip.
# This creates a list that can be used to merge all the *.gtf into the same condition 

ls OB*.gtf > mergelist.txt

stringtie --merge \
  -o merged.gtf \
  mergelist.txt

## 2. Check number of novel transcripts vs annotated ones
chmod +x De_novo_transcripts.py
./De_novo_transcripts.py

# -G gencode.v38.primary_assembly.annotation.gtf \ this argument is optional, and can be used to guide merging...
# but for de novo transcripts, its better not to use it.

## 3. Extract transcript sequences into fasta *.fa format: 
# TransDecoder only works in FASTA files, and TransDecoder is what is doing the ORF prediction (not functionally like
# short stop, but the actual ORF prediction). This produces a real nucleotide sequence for each transcript:
# >STRG.1.1
# AGTCAGTCAGTCAGTCAGTCAG... in transcript Space. More on this below

gffread merged.gtf \
  -g GRCh38.primary_assembly.genome.fa\
  -w merged.transcripts.fa

## 4. Run TransDecoder: needs FASTA (merged.transcripts.fa)
# The thing is: ShortStop works in "Genome Space", and TransDecoder in "Transcript Space"
# It produces a FASTA file of predicted ORFs (transcript space) and a GFF3 file with
# the transcript id and positions within the transcript, NOT THE GENOME

TransDecoder.LongOrfs -t merged.transcripts.fa

TransDecoder.Predict -t merged.transcripts.fa

	# That produces:
	#	-> merged.transcripts.fa.transdecoder.pep (predicted protein sequences)
	#	-> merged.transcripts.fa.transdecoder.gff3 (ORF coordinates in transcripts)

## 5. Filter ORFs to small ORFs (10–150 aa) using python script. 
# This lets us focus on  real smORFs. The resulting file is a subset of TransDecoder ORFs,
# but it is still transcript-based, which won't work with ShortStop. Still, this is a nice 
# filtering step to make the file executable and run the script:

chmod +x filter_smorf_pep.py

./filter_smorf_pep.py \
    merged.transcripts.fa.transdecoder.pep \
    --min_len 10 \
    --max_len 150 \
    --out_fasta smorfs.fa \
    --out_ids smorf_ids.txt 

grep -F -f smorf_ids.txt merged.transcripts.fa.transdecoder.gff3 > smorfs.gff3

## 6. Run smorfs_transcript_to_genome_gtf.py to converte smorfs.gtf to a genome-based GTF

chmod +x smorfs_transcript_to_genome_gtf.py

./smorfs_transcript_to_genome_gtf.py \
  --merged_gtf merged.gtf \
  --smorfs_gff3 smorfs.gff3 \
  --out_gtf smorfs_shortstop.gtf

# Thus, in this stage, we are giving a FASTA with chr1, chr2...
# a smORF GTF with genome coordinates for the transcripts and whose seqnames are chromosomes (just like what ShortStop wants)

## 7. Running the demo first downloads the GTF and FASTA formated specifically how it likes it
shortstop demo 

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$CONDA_PREFIX/lib64"

# 8. Finally, run ShortStop feature extraction:

shortstop feature_extract \
  --genome GRCh38.primary_assembly.genome.fa \
  --putative_smorfs_gtf smorfs_shortstop.gtf \
  --outdir shortstop_output \
  --threads 8

        #   --positive_gtf gencode.v43.primary_assembly.basic.annotation.gtf \
        #   do not use this line, its better to run the demo first

# Then, predict SAM-PRISM scores for the smORFs:

shortstop predict \
  --genome GRCh38.primary_assembly.genome.fa \
  --putative_smorfs_gtf smorfs_shortstop.gtf \
  --outdir shortstop_output \
  --threads 8

###########

python gtf_patcher.py \
  smorfs_shortstop.gtf \
  OB3.gtf \
  smorfs_patched.gtf
