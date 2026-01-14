```text
          ____            _   _     _                      
 _   _   |  __ \         | | (_)   | |                     
| | | |  | |__) |__ _ __ | |_ _  __| | ___  _ __ ___   ___ 
| | | |  |  ___/ _ \ '_ \| __| |/ _` |/ _ \| '_ ` _ \ / _ \
| |_| |  | |  |  __/ |_) | |_| | (_| | (_) | | | | | |  __/
| ___/   |_|   \___| .__/ \__|_|\__,_|\___/|_| |_| |_|\___|
| |                | |                                       
|_|                |_|                  

```

## Do you have data from transcriptomic experiments?

This workflow processes **RNA-seq** data to identify small open reading frames (smORFs) encoding putative microproteins (10-150 amino acids).
The pipeline performs transcriptome assembly, ORF prediction, size-based filtering, and coordinate transformation to generate genome-based annotations compatible with **ShortStop** for functional prediction. 

---


## Procedure

> Following quality control and alignment to the reference genome using STAR, StringTie performs reference-guided transcriptome assembly on sorted BAM files. For each sample, StringTie generates a GTF file containing both reference-annotated transcripts and novel isoforms discovered in the data:
> ```bash
> stringtie sample.bam \
>  -G gencode.v38.primary_assembly.annotation.gtf \
>  -o sample.gtf \
>  -p 8
> ```


- The *De_novo_transcripts.py* script quantifies novel versus annotated transcripts by parsing the ref_gene_id attribute in the merged GTF. Transcripts lacking this attribute are classified as novel:
  ```bash
chmod +x De_novo_transcripts.py
./De_novo_transcripts.py
  ```
This provides assembly quality metrics, reporting total transcripts, annotated transcripts, novel transcripts, and displaying examples of novel transcript coordinates.

- Transcript sequences are extracted from the merged GTF using *gffread*, which splices exons according to the annotation and generates a multi-FASTA file where each entry represents a complete spliced transcript:
  ```bash
gffread merged.gtf \
  -g GRCh38.primary_assembly.genome.fa \
  -w merged.transcripts.fa
  ```

- *TransDecoder* predicts ORFs in transcript coordinates through a two-step process. First, candidate ORFs meeting minimum length criteria are identified:
  ```bash
TransDecoder.LongOrfs -t merged.transcripts.fa
  ```
Then, predictions are refined based on coding potential scores and homology evidence:
``bash
TransDecoder.Predict -t merged.transcripts.fa
  ```
This generates:

 -- *merged.transcripts.fa.transdecoder.pep*: Predicted protein sequences
 -- *merged.transcripts.fa.transdecoder.gff3*: ORF coordinates in transcript space




## License and Contributions

This project is licensed for **non-commercial academic research use only**.  
See [LICENSE.md](./LICENSE.md) for full terms.

By contributing to this repository, you agree to the [Contributor License Agreement (CLA)](./CLA.md).  


By downloading or using this tool, you agree to the terms in LICENSE.md and CLA.md.
