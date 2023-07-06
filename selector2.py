import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


# Load data from CSV file into a DataFrame
df = pd.read_csv('identifier_to_subspecies.csv')

# Get list of all subspecies (ignoring 'Identifier' and 'unknown')
subspecies = [col for col in df.columns if col not in ('Identifier', 'unknown')]

# For each subspecies, get the identifiers that are unique to it and their counts
unique_identifiers = {}
unique_counts = {}
for sp in subspecies:
    # Check if all other subspecies (excluding 'unknown') are 0
    other_subspecies = [col for col in subspecies if col != sp]
    condition = (df[sp] > 0) & (df[other_subspecies].sum(axis=1) == 0)
    unique_identifiers[sp] = df.loc[condition, 'Identifier'].tolist()
    unique_counts[sp] = df.loc[condition, sp].tolist()

# This list will store all sequences from every subspecies for a final multi-fasta file
all_filtered_sequences = []

# Filter and write sequences for each subspecies
for subspecies, identifiers in unique_identifiers.items():
    if identifiers:  # Check if the list is not empty
        identifiers_set = set(identifiers)  # Convert to set for faster lookup
        output_file = f"filtered_{subspecies.replace(' ', '_')}_sequences.fasta"
        sequences = SeqIO.parse("all_strains_unique_RECO20.fasta", "fasta")

        # Filter sequences and modify sequence headers
        filtered_sequences = []
        for seq in sequences:
            if seq.id in identifiers_set:
                index = identifiers.index(seq.id)
                count = unique_counts[subspecies][index]
                # Remove spaces in the new_id
                new_id = f"{subspecies.replace(' ', '')}_{seq.id}_CN{count}"
                seq.id = new_id
                seq.description = ''
                # Remove 'X' characters from the sequence
                seq.seq = Seq(str(seq.seq).replace('X', ''))
                filtered_sequences.append(seq)
                all_filtered_sequences.append(seq)  # Add sequence to the multi-fasta sequence list

        # Write sequences to file
        SeqIO.write(filtered_sequences, output_file, "fasta")
        print(f"Filtered sequences for {subspecies} saved to {output_file}.")

# Finally, write the multi-fasta file with all filtered sequences
SeqIO.write(all_filtered_sequences, "filtered_Trypanozoon_sequences.fasta", "fasta")
print("Filtered sequences for all subspecies saved to filtered_Trypanozoon_sequences.fasta.")
