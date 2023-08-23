#!/usr/bin/env python
import os
import argparse
import csv
import glob
import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def search_and_reorganize_sequences(args):
    fasta_path, base_output_dir, cn_threshold = args
    # Create a new output directory based on the copy number
    output_dir = os.path.join(base_output_dir, f'CN_{cn_threshold}')
    strain_name = os.path.basename(os.path.dirname(fasta_path))
    print(f"Processing strain: {strain_name}")

    selected_sequences = {}  # Store sequences >= 170 bp and CN > threshold for each strain

    sequence_dict = defaultdict(int)
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequence = str(record.seq)
        if len(sequence) >= 170:
            sequence = sequence.replace("X", "")  # Remove the nucleotide "X" from the sequence
            sequence_dict[sequence] += 1
 
    for sequence, count in sequence_dict.items():
        if count > cn_threshold:  # Use the provided threshold
            selected_sequences[sequence] = count
  
    # Return the strain name along with its selected sequences
    return strain_name, selected_sequences

def identify_unique_sequences(all_strain_sequences):
    combined_sequences = defaultdict(list)
    for strain, sequences in all_strain_sequences:
        for seq, count in sequences.items():
            combined_sequences[seq].append((strain, count))

    unique_sequences = {seq: strains for seq, strains in combined_sequences.items()}
    return unique_sequences


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate stats from FASTA TBR reads.')
    parser.add_argument('-o', '--output', help='The base output directory', required=True)
    parser.add_argument('-d', '--directory', help='The directory with fasta files', required=True)
    parser.add_argument('-cn', '--copynumber', type=int, help='The copy number threshold', required=True)
    args = parser.parse_args()  # Capture parsed arguments
    output_dir = os.path.join(args.output, f'CN_{args.copynumber}')
    os.makedirs(output_dir, exist_ok=True)

    fasta_files = glob.glob(os.path.join(args.directory, '*', '*.TBR.bam.fasta'))

    with mp.Pool(mp.cpu_count()) as pool:
        all_strain_sequences = pool.map(search_and_reorganize_sequences, [(fasta_file, args.output, args.copynumber) for fasta_file in fasta_files])
   
    unique_sequences = identify_unique_sequences(all_strain_sequences)

    # Write the output files
    with open(os.path.join(output_dir, f'unique_sequences_CN{args.copynumber}.fasta'), 'w') as f:
        for idx, (sequence, strains) in enumerate(unique_sequences.items(), 1):
            f.write(f">U{idx}\n{sequence}\n")

    with open(os.path.join(output_dir, f'copy_number_lookup_CN{args.copynumber}.csv'), 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Unique_ID', 'Strain', 'CopyNumber'])
        for idx, (sequence, strains) in enumerate(unique_sequences.items(), 1):
            for strain, count in strains:
                writer.writerow([f"U{idx}", strain, count])
