import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Load data from CSV file into a DataFrame
df = pd.read_csv('identifier_to_strain.csv')

# Get list of all subspecies_strains (ignoring 'Identifier')
subspecies_strains = [col for col in df.columns if col != 'Identifier']

# Identify abbreviations of subspecies
subspecies_abbreviations = set([s.split('_')[0] for s in subspecies_strains if s != 'Uk'])

# For each subspecies, get the identifiers that are unique to it
unique_identifiers = {}
for sp in subspecies_abbreviations:
    # Generate list of strains for this subspecies
    strains_for_this_sp = [s for s in subspecies_strains if s.split('_')[0] == sp]

    # Check if all other subspecies (excluding 'Uk') are 0
    other_subspecies_strains = [s for s in subspecies_strains if s.split('_')[0] != sp and s != 'Uk']
    condition = (df[strains_for_this_sp].sum(axis=1) > 0) & (df[other_subspecies_strains].sum(axis=1) == 0)
    unique_identifiers[sp] = df.loc[condition, 'Identifier'].tolist()

# This list will store all sequences from every subspecies_strain for a final multi-fasta file
all_filtered_sequences = []

# Filter and write sequences for each subspecies_strain
for subspecies_strain in subspecies_strains:
    identifiers = unique_identifiers.get(subspecies_strain.split('_')[0], [])  # Check if the list is not empty
    if identifiers:
        identifiers_set = set(identifiers)  # Convert to set for faster lookup
        output_file = f"filtered_{subspecies_strain.replace(' ', '_')}_sequences.fasta"
        sequences = SeqIO.parse("all_strains_unique_RECO20.fasta", "fasta")

        # Filter sequences and modify sequence headers
        filtered_sequences = []
        for seq in sequences:
            if seq.id in identifiers_set:
                # Get count for this subspecies_strain
                count = df.loc[df['Identifier'] == seq.id, subspecies_strain].values[0]
                # Remove spaces in the new_id
                new_id = f"{subspecies_strain.replace(' ', '_')}_{seq.id}_CN{count}"
                seq.id = new_id
                seq.description = ''
                # Remove 'X' characters from the sequence
                seq.seq = Seq(str(seq.seq).replace('X', ''))
                filtered_sequences.append(seq)
                all_filtered_sequences.append(seq)  # Add sequence to the multi-fasta sequence list

        # Write sequences to file
        SeqIO.write(filtered_sequences, output_file, "fasta")
        print(f"Filtered sequences for {subspecies_strain} saved to {output_file}.")

# Finally, write the multi-fasta file with all filtered sequences
SeqIO.write(all_filtered_sequences, "filtered_all_sequences.fasta", "fasta")
print("Filtered sequences for all subspecies saved to filtered_all_sequences.fasta.")
