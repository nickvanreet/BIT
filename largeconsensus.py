fasta_file = "output._gamb2028d07.p1k.fasta.fasta"  # Replace with the actual filename

sequences = {}  # Dictionary to store sequence sizes

with open(fasta_file, "r") as file:
    sequence_name = ""
    sequence = ""
    for line in file:
        if line.startswith(">"):
            if sequence_name and sequence:
                sequences[sequence_name] = sequence
            sequence_name = line[1:].strip()
            sequence = ""
        else:
            sequence += line.strip()
    # Add the last sequence to the dictionary
    if sequence_name and sequence:
        sequences[sequence_name] = sequence

# Find the sequence with the largest size
largest_sequence_name = max(sequences, key=lambda x: len(sequences[x]))

# Print the extracted sequence
print(largest_sequence_name)
print(sequences[largest_sequence_name])