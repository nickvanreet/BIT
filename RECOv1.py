import os
import glob
import argparse
from collections import defaultdict
from Bio import SeqIO

def create_multifasta(output_dir, output_filename, prefix, target_sequence):
    files = glob.glob(os.path.join(output_dir, f'{prefix}*_{target_sequence}_*.fasta'))
    print(f"Creating multifasta file: {output_filename}")
    print(f"Files to be merged: {files}")

    # Construct the correct filenames
    files = [filename for filename in files if prefix in filename]
    print(f"Files to be merged after filtering: {files}")

    with open(output_filename, 'w') as outfile:
        for filename in files:
            with open(filename, 'r') as readfile:
                outfile.write(readfile.read())
            os.remove(filename)  # remove the file after its content has been written to the output file
    return len(files)


def search_and_reorganize_sequences(variant_output_filename, target_sequence, output_dir):
    print(f"Processing file: {variant_output_filename}")
    fasta_entries = SeqIO.parse(variant_output_filename, 'fasta')
    reconstructed_sequences = []
    sequences_not_found = []
    fasta_filename = os.path.basename(variant_output_filename)

    for entry in fasta_entries:
        header = entry.id
        sequence = str(entry.seq)

        mismatch_position = -1
        match_found = False
        sequence_length = len(sequence)
        circular_sequence = sequence + sequence  # Circularize the sequence by concatenating it with itself

        for i in range(sequence_length):
            window = circular_sequence[i:i + len(target_sequence)]
            mismatches = sum([1 for a, b in zip(target_sequence, window) if a != b])

            if mismatches <= 1:
                if mismatches == 1:
                    mismatch_position = i + [j for j, (a, b) in enumerate(zip(target_sequence, window)) if a != b][0]

                match_found = True
                reconstructed_sequence = circular_sequence[i:i + sequence_length]
                reconstructed_sequences.append((header, reconstructed_sequence, i, mismatch_position))
                break

        if not match_found:
            sequences_not_found.append((header, sequence))

    output_file = os.path.join(output_dir, f'reorganized_{target_sequence}_{fasta_filename}')
    not_found_file = os.path.join(output_dir, f'sequences_not_found_{target_sequence}_{fasta_filename}')

    # Save reconstructed sequences to a new FASTA file
    with open(output_file, 'w') as f:
        for header, sequence, start_position, mismatch_position in reconstructed_sequences:
            reconstructed_sequence_length = len(sequence)
            if mismatch_position == -1:
                f.write(f'>{header}|Matched at: {start_position}|Original length: {len(sequence)}|New length: {reconstructed_sequence_length}\n{sequence}\n')
            elif mismatch_position != -1:
                f.write(f'>{header}|Mismatched at: {mismatch_position}|Original length: {len(sequence)}|New length: {reconstructed_sequence_length}\n{sequence}\n')

    # Save sequences not found to a separate file
    with open(not_found_file, 'w') as f:
        for header, sequence in sequences_not_found:
            f.write(f'>{header}|Target sequence not found|Original length: {sequence_length}\n{sequence}\n')

    reconstructed_count = len(reconstructed_sequences)
    not_found_count = len(sequences_not_found)
    print(f'Total number of variants reconstructed: {reconstructed_count}')
    print(f'Total number of variants not found: {not_found_count}')


def write_sequences(output_file, sequences, strain):
    sequence_counts = defaultdict(int)

    for sequence in sequences:
        sequence_counts[str(sequence.seq)] += 1

    sequences_with_counts = [(count, seq) for seq, count in sequence_counts.items()]
    sequences_with_counts.sort(reverse=True)  # Sort sequences based on count in descending order

    with open(output_file, 'w') as file:
        for count, sequence in sequences_with_counts:
            header = f">{strain} | Count: {count}"
            file.write(f"{header}\n")
            file.write(f"{sequence}\n")


def process_sequences(input_dir, output_dir, target_sequence):
    # Check if the output directory exists and create it if necessary
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    multifasta_files = glob.glob(os.path.join(input_dir, '*_multi_variants.fasta'))

    total_file_count = 0  # Initialize the variable here

    for multifasta_file in multifasta_files:
        search_and_reorganize_sequences(multifasta_file, target_sequence, output_dir)

    # Create the multifasta file for reorganized files and remove individual FASTA files
    output_reorganized_file = os.path.join(output_dir, f'reorg_multi_{target_sequence}.fasta')
    reorganized_file_count = create_multifasta(output_dir, output_reorganized_file, 'reorganized', target_sequence)
    total_file_count += reorganized_file_count

    # Create the multifasta file for sequences not found files and remove individual FASTA files
    output_sequences_not_found_file = os.path.join(output_dir, f'not_found_multi_{target_sequence}.fasta')
    sequences_not_found_file_count = create_multifasta(output_dir, output_sequences_not_found_file, 'sequences_not_found', target_sequence)
    total_file_count += sequences_not_found_file_count

    # Create a multifasta file for unique and non-unique sequences
    unique_sequences = defaultdict(list)
    non_unique_sequences = []

    # Read sequences from the file
    with open(output_reorganized_file, "r") as file:
        sequences = SeqIO.parse(file, "fasta")

        # Iterate through each sequence
        for sequence in sequences:
            # Determine the strain based on the filename
            filename = sequence.description.split("|")[0].strip().split("_")
            strain = filename[0]

            # Add the sequence to the corresponding strain dictionary
            unique_sequences[strain].append(sequence)

    # Write unique sequences to multifasta files for each strain
    for strain, sequences in unique_sequences.items():
        unique_output_file = os.path.join(output_dir, f"unique_sequences_{strain}.fasta")
        write_sequences(unique_output_file, sequences, strain)
        print(f"Number of unique sequences in {strain}: {len(sequences)}")
        print(f"Unique sequences saved in: {unique_output_file}")

    # Write non-unique sequences to a multifasta file
    non_unique_output_file = os.path.join(output_dir, f"non_unique_sequences_{target_sequence}.fasta")
    write_sequences(non_unique_output_file, non_unique_sequences, target_sequence)
    print(f"Number of non-unique sequences: {len(non_unique_sequences)}")
    print(f"Non-unique sequences saved in: {non_unique_output_file}")

    # Create an overview of unique sequences found in each strain
    overview_file = os.path.join(output_dir, "unique_sequences_overview.txt")
    with open(overview_file, "w") as file:
        file.write("Strain\tSequence ID\n")
        for strain, sequences in unique_sequences.items():
            for sequence in sequences:
                sequence_id = sequence.id
                file.write(f"{strain}\t{sequence_id}\n")

    print(f"Unique sequences overview saved in: {overview_file}")
    print(f"Total number of processed files: {total_file_count}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process multi-FASTA files to search and reorganize sequences.')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing multi-FASTA files')
    parser.add_argument('-o', '--output', required=True, help='Output directory for results')
    parser.add_argument('-t', '--target', required=True, help='Target sequence to search in the multi-FASTA files')

    args = parser.parse_args()

    process_sequences(args.input, args.output, args.target)
