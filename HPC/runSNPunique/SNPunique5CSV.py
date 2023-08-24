#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from collections import defaultdict
import csv

CHUNK_SIZE = 10000

def read_strain_lookup(file_path):
    strain_lookup = {}
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            identifier, strain, copy_number = line.strip().split(',')
            strain_lookup[(identifier, strain)] = int(copy_number)
    return strain_lookup

def read_subspecies_lookup(file_path):
    subspecies_lookup = {}
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            data = line.strip().split(';')
            strain = data[0]
            subspecies = data[2]
            subspecies_lookup[strain] = subspecies
    return subspecies_lookup

def process_record(seq_record, strain_lookup, subspecies_lookup, agg_dict):
    for strain, subspecies in subspecies_lookup.items():
        identifier = "U" + seq_record.description.split()[0][1:]
        if (identifier, strain) in strain_lookup:
            copy_number = strain_lookup[(identifier, strain)]
            for i, nucleotide in enumerate(str(seq_record.seq)):
                if nucleotide in "agtcnx-":
                    key = (i + 1, nucleotide, subspecies, strain)
                    agg_dict[key] += copy_number

def iterative_processing(input_file, strain_lookup, subspecies_lookup, aggregated_path, strain_path):
    strains_count = defaultdict(int)
    strains_seen = defaultdict(set)
    agg_dict = defaultdict(int)
    
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            process_record(record, strain_lookup, subspecies_lookup, agg_dict)

    with open(aggregated_path, 'w', newline='') as agg_handle:
        csv_writer = csv.writer(agg_handle)
        csv_writer.writerow(['Position', 'Nucleotide', 'Subspecies', 'Strain', 'Count'])
        for (position, nucleotide, subspecies, strain), count in agg_dict.items():
            csv_writer.writerow([position, nucleotide, subspecies, strain, count])
            if strain not in strains_seen[subspecies]:
                strains_count[subspecies] += 1
                strains_seen[subspecies].add(strain)

    with open(strain_path, 'w', newline='') as sp_handle:
        csv_writer = csv.writer(sp_handle)
        csv_writer.writerow(['Subspecies', 'Number of Strains'])
        for subspecies, count in strains_count.items():
            csv_writer.writerow([subspecies, count])
            
def main():
    parser = argparse.ArgumentParser(description='Count SNP')
    parser.add_argument('-i', '--input', help='Input fasta file name', required=True)
    parser.add_argument('-o', '--output', help='Output CSV file name for aggregated data by subspecies', required=True)
    parser.add_argument('-a', '--aggregated', help='Output CSV file name for aggregated data by subspecies and strain', required=True)
    parser.add_argument('-l', '--lookup', help='Lookup file for strains and copy number', required=True)
    parser.add_argument('-s', '--subspecies_lookup', help='Lookup file for subspecies', required=True)
    parser.add_argument('-sp', '--strain_path', help='Output CSV file name for strains per subspecies', required=True)
    args = parser.parse_args()

    strain_lookup = read_strain_lookup(args.lookup)
    subspecies_lookup = read_subspecies_lookup(args.subspecies_lookup)

    iterative_processing(args.input, strain_lookup, subspecies_lookup, args.aggregated, args.strain_path)

if __name__ == "__main__":
    main()
