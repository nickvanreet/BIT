import os
import subprocess
import re

def run_trf(sequence_file):
    trf_path = "/usr/local/bin/trf"  # Path to the TRF executable

    # Run TRF on the sequence file and capture the output
    command = f"{trf_path} {sequence_file} 2 5 7 80 10 50 2000 -d -h"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    return result.stdout


def extract_repeated_sequences(trf_output):
    repeated_sequences = []

    # Extract repeated sequences using regular expressions
    pattern = r"^\d+\s\d+\s(\d+)\s(.+)$"
    matches = re.findall(pattern, trf_output, re.MULTILINE)
    for match in matches:
        sequence_length = match[0]
        sequence = match[1]
        repeated_sequences.append(sequence)

    return repeated_sequences

def analyze_fasta_sequence(sequence_file):
    trf_output = run_trf(sequence_file)
    repeated_sequences = extract_repeated_sequences(trf_output)
    return repeated_sequences

# Example usage
fasta_file = "gamb1801e06p1k.fasta"

repeated_sequences = analyze_fasta_sequence(fasta_file)

# Print the extracted repeated sequences
for sequence in repeated_sequences:
    print(sequence)
