from Bio import SeqIO
from collections import defaultdict
import pandas as pd

possible_nucleotides = ['A', 'T', 'C', 'G', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '-']

def parse_sequence_file(file_path):
    with open(file_path, "r") as handle:
        sequences = list(SeqIO.parse(handle, "fasta"))
    return sequences

def find_tbg_specific_snps(sequences):
    snp_positions = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for seq_record in sequences:
        description_fields = seq_record.description.split('|')
        hits_field = None
        for field in description_fields:
            if "Hits:" in field:
                hits_field = field
                break
        if hits_field:
            hits = [hit.strip() for hit in hits_field.split(':')[1].split('&')]
            for i, nucleotide in enumerate(str(seq_record.seq)):
                for hit in hits:
                    hit_parts = hit.split('_')
                    if len(hit_parts) >= 2:
                        species, strain = hit_parts[:2]
                    else:
                        species = hit_parts[0]
                        strain = 'unknown'
                    snp_positions[i+1][(species, strain)][nucleotide] += 1

    return snp_positions

def write_to_file(tbg_specific_snps, output_file):
    data = []
    for position, species_strains in tbg_specific_snps.items():
        for species_strain, snps in species_strains.items():
            for nucleotide, count in snps.items():
                data.append([position, nucleotide, '_'.join(species_strain), count])

    df = pd.DataFrame(data, columns=['Position', 'Nucleotide', 'Species_Strain', 'Count'])
    df = df.pivot_table(index=['Position', 'Nucleotide'], columns='Species_Strain', values='Count', fill_value=0)
    df.to_csv(output_file)

sequences = parse_sequence_file('GCGCAGTT_unique_sequences_175_180bp.fasta')
tbg_specific_snps = find_tbg_specific_snps(sequences)
write_to_file(tbg_specific_snps, 'SNP_counts.csv')
