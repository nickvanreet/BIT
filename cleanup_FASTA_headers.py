import os
import argparse
from Bio import SeqIO

def shorten_fasta_headers(input_dir):
    fasta_files = [filename for filename in os.listdir(input_dir) if filename.endswith(".fasta")]

    for filename in fasta_files:
        fasta_file = os.path.join(input_dir, filename)
        
        # Read the FASTA file and modify the headers
        sequences = SeqIO.parse(fasta_file, "fasta")
        shortened_sequences = []
        for seq in sequences:
            seq.id = seq.id.split()[0]  # Modify this line based on your needs
            seq.description = ''  # Clear the description to prevent it from being added to the header
            shortened_sequences.append(seq)
        
        # Overwrite the original file with the sequences having shortened headers
        SeqIO.write(shortened_sequences, fasta_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Shorten FASTA headers.')
    parser.add_argument('input_dir', help='Path to the directory containing the FASTA files.')
    args = parser.parse_args()

    shorten_fasta_headers(args.input_dir)
