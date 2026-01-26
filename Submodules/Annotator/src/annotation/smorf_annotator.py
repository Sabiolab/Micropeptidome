import re
import argparse
from ..pipeline import PipelineStructure
from collections import defaultdict

# === Precompiled Regex for GTF/GFF attributes ===
ATTR_RE = re.compile(r'(\S+) "([^"]+)"')

def parse_attributes(attr_str):
    return dict(ATTR_RE.findall(attr_str))

class smORFAnnotator(PipelineStructure):
    def __init__(self, args):
        super().__init__(args)
        self.intersect_file = args.intersect_output
        self.non_intersect_file = args.non_intersect_output
        self.output_file = args.output_file

    def process_gtf_files(self):
        priority_order = ['psORF', 'uoORF','doORF','oORF', 'dORF', 'uORF', 'lncRNA', 'riORF', 'aORF', 'eORF']

        def get_priority(annot):
            return priority_order.index(annot) if annot in priority_order else float('inf')

        gene_data = defaultdict(lambda: ('UA', 'Unknown'))

        PSEUDO_BIOTYPES = {
            'processed_pseudogene', 'unprocessed_pseudogene', 'translated_unprocessed_pseudogene',
            'translated_processed_pseudogene', 'transcribed_processed_pseudogene',
            'transcribed_unprocessed_pseudogene', 'unitary_pseudogene', 'polymorphic_pseudogene'
        }

        NONCODING_BIOTYPES = {'lncRNA', 'lincRNA', 'antisense', 'sense_intronic', 'sense_overlapping'}

        with open(self.intersect_file, 'r') as file:
            for line in file:
                if not line.startswith('chr'):
                    continue
                parts = line.strip().split('\t')
                if parts[2] != 'CDS':
                    continue

                attrs = parse_attributes(parts[8])
                gene_info = parts[17] if len(parts) > 17 else parts[8]
                gene_attrs = parse_attributes(gene_info)

                gene_id = attrs.get('gene_id', 'Unknown')
                gene_name = gene_attrs.get('gene_name', 'Unnamed')
                gene_biotype = gene_attrs.get('gene_biotype', 'Unnamed')
                transcript_biotype = gene_attrs.get('transcript_biotype', 'Unknown')
                feature_type = parts[11] if len(parts) > 11 else 'exon'

                # === Annotation Logic ===
                if feature_type == 'three_prime_utr':
                    annotation = 'dORF'
                elif feature_type == 'five_prime_utr':
                    annotation = 'uORF'
                elif feature_type == 'CDS':
                    annotation = 'oORF'
                elif feature_type == 'retrotransposed':
                    annotation = 'psORF'
                elif transcript_biotype in PSEUDO_BIOTYPES:
                    annotation = 'psORF'
                elif any(x in transcript_biotype for x in [
                    'non_stop_decay', 'nonsense_mediated_decay',
                    'ambiguous_orf', 'protein_coding_CDS_not_defined']):
                    annotation = 'aORF'
                elif 'retained_intron' in transcript_biotype:
                    annotation = 'riORF'
                elif gene_biotype in NONCODING_BIOTYPES:
                    annotation = 'lncRNA'
                elif feature_type == 'exon':
                    annotation = 'eORF'
                else:
                    annotation = 'UA'

                # === Conflict Resolution ===
                existing_annotation, _ = gene_data[gene_id]
                
                if gene_id not in gene_data:
                    gene_data[gene_id] = (annotation, gene_name)
                else:
                    existing_annotation, _ = gene_data[gene_id]

                    # Handle combined upstream + overlapping
                    if {existing_annotation, annotation} == {'uORF', 'oORF'}:
                        gene_data[gene_id] = ('uoORF', gene_name)

                    # Handle combined downstream + overlapping
                    elif {existing_annotation, annotation} == {'dORF', 'oORF'}:
                        gene_data[gene_id] = ('doORF', gene_name)

                    # Handle upstream + downstream (rare)
                    elif {existing_annotation, annotation} == {'uORF', 'dORF'}:
                        gene_data[gene_id] = ('udORF', gene_name)

                    # Keep whichever has higher priority (lower index)
                    elif get_priority(annotation) < get_priority(existing_annotation):
                        gene_data[gene_id] = (annotation, gene_name)
                
                if get_priority(annotation) < get_priority(existing_annotation):
                    gene_data[gene_id] = (annotation, gene_name)

        # === Add Intergenic Genes ===
        with open(self.non_intersect_file, 'r') as file:
            for line in file:
                if not line.startswith('chr'):
                    continue
                parts = line.strip().split('\t')
                attrs = parse_attributes(parts[8])
                gene_id = attrs.get('gene_id', 'Unknown')
                if gene_id not in gene_data:
                    gene_data[gene_id] = ('Intergenic', 'Intergenic')

        # === Write to Output File ===
        with open(self.output_file, 'w') as out:
            for gene_id, (annotation, gene_name) in gene_data.items():
                out.write(f'{gene_id}\t{annotation}\t{gene_name}\n')

        print(f"[✓] Output written to {self.output_file}")

        # === Optional: Print summary ===
        summary = defaultdict(int)
        for annot, _ in gene_data.values():
            summary[annot] += 1

        print("\n[📊] Annotation summary:")
        for annot in sorted(summary.keys()):
            print(f"  {annot:12} : {summary[annot]}")