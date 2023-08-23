import argparse
from Bio import SeqIO
from collections import defaultdict
import re
import pandas as pd

def parse_sequence_file(file_path):
    with open(file_path, "r") as handle:
        sequences = []
        for record in SeqIO.parse(handle, "fasta"):
            match = re.search(r"(\w+)_U(\d+)", record.id)
            if match is None:
                continue
            strain, unique_index = match.groups()
            sequences.append((record, strain, int(unique_index)))
        return sequences

def process_record(seq_data):
    seq_record, strain, unique_index = seq_data
    snp_positions = defaultdict(lambda: defaultdict(int))
    for i, nucleotide in enumerate(str(seq_record.seq)):
        snp_positions[i+1][nucleotide] += 1
    return strain, snp_positions

def find_strain_specific_snps(sequences):
    strain_specific_snps = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for seq_data in sequences:
        strain, snp_positions = process_record(seq_data)
        for position, nucleotides in snp_positions.items():
            for nucleotide, count in nucleotides.items():
                strain_specific_snps[strain][position][nucleotide] += count
    return strain_specific_snps

def save_to_csv(strain_specific_snps, output_file):
    data = {(strain, position, nucleotide): count 
            for strain, positions in strain_specific_snps.items() 
            for position, nucleotides in positions.items() 
            for nucleotide, count in nucleotides.items()}
    df = pd.Series(data).unstack(level=0).fillna(0).astype(int)  # Convert to int here
    df.columns.name = None
    df.to_csv(output_file)


# Argument parsing and main
parser = argparse.ArgumentParser(description='Count SNP')
parser.add_argument('-i', '--input', help='Input file name', default='filtered')
parser.add_argument('-o', '--output', help='Output CSV file name', required=True)
args = parser.parse_args()

sequences = parse_sequence_file(args.input)
strain_specific_snps = find_strain_specific_snps(sequences)
save_to_csv(strain_specific_snps, args.output)
