#!/usr/bin/env python
import os
import argparse
import csv
import glob
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def remove_indels(sequence):
    return sequence.replace('-', '')

def search_and_reorganize_sequences(alignment_fasta, target_sequence, output_dir):
    print(f"Processing alignment file: {alignment_fasta}")

    spacer = 'N'*27  # Define spacer
    reconstructed_sequences = []
    sequences_not_found = []

    fasta_entries = SeqIO.parse(alignment_fasta, 'fasta')

    for entry in fasta_entries:
        header = entry.id
        sequence = str(entry.seq)

        mismatch_position = -1
        match_found = False
        sequence_length = len(sequence)
        spacer_length = len(spacer)

        circular_sequence = sequence + spacer + sequence  # Circularize the sequence by concatenating it with itself

        for i in range(sequence_length):
            window = circular_sequence[i:i + len(target_sequence)]
            mismatches = sum([1 for a, b in zip(target_sequence, window) if a != b])

            if mismatches <= 1:
                if mismatches == 1:
                    mismatch_position = i + [j for j, (a, b) in enumerate(zip(target_sequence, window)) if a != b][0]

                match_found = True
                reconstructed_sequence = circular_sequence[i:i + sequence_length + spacer_length]
                reconstructed_sequence = remove_indels(reconstructed_sequence)  # Remove indels after reconstruction
                reconstructed_sequences.append((header, reconstructed_sequence, i, mismatch_position))
                break

        if not match_found:
            sequences_not_found.append((header, sequence))

    output_file = os.path.join(output_dir, f'reorganized_{target_sequence}.fasta')
    not_found_file = os.path.join(output_dir, f'sequences_not_found_{target_sequence}.fasta')

    # Write reconstructed sequences to output file
    with open(output_file, 'w') as f:
        for header, seq, _, _ in reconstructed_sequences:
            f.write(f'>{header}\n{seq}\n')

    # Write sequences not found to a separate file
    with open(not_found_file, 'w') as f:
        for header, seq in sequences_not_found:
            f.write(f'>{header}\n{seq}\n')

    return output_dir


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-o', '--output', help='The output directory', required=True)
    parser.add_argument('-f', '--fasta', help='The alignment fasta file', required=True)
    parser.add_argument('-t', '--target', help='The target sequence', required=True)

    args = parser.parse_args()

    search_and_reorganize_sequences(args.fasta, args.target, args.output)
