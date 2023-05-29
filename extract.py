import re

def extract_sequences(trf_output_file):
    sequences = []

    with open(trf_output_file, 'r') as file:
        # Read the TRF output file
        lines = file.readlines()

    # Process each line in the TRF output
    for line in lines:
        if not line.startswith('@'):
            fields = line.strip().split()
            if len(fields) >= 14:
                sequence = fields[13]
                length = len(sequence)
                if length == 177 or length == 176:
                    sequences.append(sequence)

    return sequences

trf_output_file = "test.fasta.2.5.7.80.10.50.2000.dat"

repeated_sequences = extract_sequences(trf_output_file)

# Print the extracted sequences
for sequence in repeated_sequences:
    print(sequence)
