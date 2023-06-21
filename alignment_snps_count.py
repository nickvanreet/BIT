from Bio import SeqIO
from collections import defaultdict
import os

def parse_sequence_file(file_path):
    with open(file_path, "r") as handle:
        sequences = list(SeqIO.parse(handle, "fasta"))
    return sequences

def find_tbg_specific_snps(sequences):
    snp_positions = defaultdict(lambda: defaultdict(int))
    tbg_specific_snps = defaultdict(list)
    possible_nucleotides = ['A', 'T', 'C', 'G', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '-']

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
                    snp_positions[i+1][hit.split('_')[0] + '_' + nucleotide] += 1

    for position, snps in snp_positions.items():
        for possible_nucleotide in possible_nucleotides:
            tbg_count = sum(count for snp, count in snps.items() if 'Tbg_' + possible_nucleotide in snp)
            tbb_tbr_count = sum(count for snp, count in snps.items() if 'Tbb_' + possible_nucleotide in snp or 'Tbr_' + possible_nucleotide in snp)
            tbg_specific_snps[position].append((possible_nucleotide, tbg_count, tbb_tbr_count))

    return tbg_specific_snps

def write_to_file(tbg_specific_snps, output_file):
    with open(output_file, 'w') as file:
        file.write('Position\tSNP in Tbg\tTotal Tbg Clones\tSNP in Tbb/Tbr\tTotal Tbb/Tbr Clones\n')
        for position, snps in tbg_specific_snps.items():
            for snp in snps:
                file.write(f'{position}\tTbg_{snp[0]}\t{snp[1]}\tTbb/Tbr_{snp[0]}\t{snp[2]}\n')

sequences = parse_sequence_file('unique_sequences_175_180bp_GCGCAGTT.fa')
tbg_specific_snps = find_tbg_specific_snps(sequences)
write_to_file(tbg_specific_snps, 'SNP_counts.txt')
