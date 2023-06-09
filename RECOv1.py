import os
import glob
import argparse
import csv
from collections import defaultdict
from Bio import SeqIO

# USAGE: python script.py -i your/input/directory -o your/output/directory -t your_target_sequence

def create_multifasta(output_dir, output_filename, prefix, target_sequence):
    files = glob.glob(os.path.join(output_dir, f'{prefix}*_{target_sequence}_*.fasta'))
    print(f"Creating multifasta file: {output_filename}")
    #print(f"Files to be merged: {files}")

    # Construct the correct filenames
    files = [filename for filename in files if prefix in filename]
    #print(f"Files to be merged after filtering: {files}")

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

    return len(reconstructed_sequences), len(sequences_not_found)

def process_sequences(input_dir, output_dir, target_sequence):
    # Check if the output directory exists and create it if necessary
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    multifasta_files = glob.glob(os.path.join(input_dir, '*_multi_variants.fasta'))
    total_file_count = 0  # Initialize the variable here
    results = []

    for multifasta_file in multifasta_files:
        filename = os.path.basename(multifasta_file)
        
        # Use the returned values from search_and_reorganize_sequences function
        number_of_variants_reconstructed, number_of_variants_not_found = search_and_reorganize_sequences(multifasta_file, target_sequence, output_dir)

        # Create a dictionary and a list to hold unique and non-unique sequences
        # Move these lines inside the loop so they are reset for each file


        # Update unique_sequences and non_unique_sequences based on the data in the file
        with open(multifasta_file, "r") as file:
            sequences = SeqIO.parse(file, "fasta")
            
            unique_sequences = defaultdict(list)
            non_unique_sequences = []
            
            for sequence in sequences:
                sequence_str = str(sequence.seq)

                if sequence_str not in unique_sequences:
                    unique_sequences[sequence_str].append(sequence)
                else:
                    non_unique_sequences.append(sequence)

            # Calculate the number of unique and non-unique sequences
            number_of_unique_sequences = len(unique_sequences)
            number_of_unique_sequences_175_180bp = sum([1 for seq in unique_sequences if 175 <= len(seq) <= 180])
            number_of_non_unique_sequences = len(non_unique_sequences)
            
            results.append([filename, target_sequence, number_of_variants_reconstructed, number_of_variants_not_found,
                            number_of_unique_sequences, number_of_unique_sequences_175_180bp, number_of_non_unique_sequences])

    # Write the results to a CSV file
    # with open(os.path.join(output_dir, 'results.csv'), 'w', newline='') as csvfile:
    #    writer = csv.writer(csvfile)
    #    writer.writerow(["Filename", "Target Sequence", "Number of Variants Reconstructed", "Number of Variants Not Found",
    #                     "Number of Unique Sequences", "Number of Unique Sequences between 175 and 180 bp", "Number of Non-Unique Sequences"])
    #    writer.writerows(results)
         
    # Create the multifasta file for reorganized files and remove individual FASTA files
    output_reorganized_file = os.path.join(output_dir, f'{target_sequence}_reorg_multi.fasta')
    reorganized_file_count = create_multifasta(output_dir, output_reorganized_file, 'reorganized', target_sequence)
    total_file_count += reorganized_file_count
    print(f"Reorganized multifasta file created: {output_reorganized_file}")
    print(f"Number of reorganized files merged: {reorganized_file_count}")

    # Create the multifasta file for sequences not found files and remove individual FASTA files
    output_sequences_not_found_file = os.path.join(output_dir, f'{target_sequence}_not_found_multi.fasta')
    sequences_not_found_file_count = create_multifasta(output_dir, output_sequences_not_found_file, 'sequences_not_found', target_sequence)
    total_file_count += sequences_not_found_file_count
    print(f"Sequences not found multifasta file created: {output_sequences_not_found_file}")
    print(f"Number of sequences not found files merged: {sequences_not_found_file_count}")

    # Create a multifasta file for unique and non-unique sequences
    unique_sequences = defaultdict(list)
    non_unique_sequences = []

    # Read sequences from the file
    with open(output_reorganized_file, "r") as file:
        sequences = SeqIO.parse(file, "fasta")

        # Iterate through each sequence
        for sequence in sequences:
            # Check if the sequence is unique
            sequence_str = str(sequence.seq)
            if sequence_str not in unique_sequences:
                unique_sequences[sequence_str] = [(sequence, sequence.description)]
            else:
                unique_sequences[sequence_str].append((sequence, sequence.description))
                non_unique_sequences.append(sequence)

    # Write unique sequences to a multifasta file and count the number of sequences
    unique_output_file = os.path.join(output_dir, f"{target_sequence}_unique_sequences.fasta")
    with open(unique_output_file, "w") as file:
        unique_sequences_list = list(unique_sequences.items())

        # Sort the unique sequences based on count (in descending order)
        sorted_unique_sequences = sorted(unique_sequences_list, key=lambda x: len(x[1]), reverse=True)

        count = 0
        for i, (sequence_str, sequences) in enumerate(sorted_unique_sequences):
            count += len(sequences)
            original_headers = [seq.description.split("|")[0].strip().lstrip(">/").split("/")[-1] for seq, _ in sequences]
            header = f"Unique_{i} | Count: {len(sequences)} | Hits: {' & '.join(original_headers)}"
            file.write(f">{header}\n{sequence_str}\n")

    print(f"Number of unique sequences: {count}")
    print(f"Unique sequences saved in: {unique_output_file}")

    # Write unique sequences 175 - 180 bp to a multifasta file and count the number of sequences
    unique_output_file_175_180bp = os.path.join(output_dir, f"{target_sequence}_unique_sequences_175_180bp.fasta")
    with open(unique_output_file_175_180bp, "w") as file:
        count = 0
        for i, (sequence_str, sequences) in enumerate(sorted_unique_sequences):
            if 175 <= len(sequence_str) <= 180:  # Only consider sequences between 175 and 180 bp
                count += len(sequences)
                original_headers = [seq.description.split("|")[0].strip().lstrip(">/").split("/")[-1] for seq, _ in sequences]
                header = f"Unique_{i} | Count: {len(sequences)} | Hits: {' & '.join(original_headers)}"
                file.write(f">{header}\n{sequence_str}\n")

    print(f"Number of unique sequences between 175 and 180 bp: {count}")
    print(f"Unique sequences between 175 and 180 bp saved in: {unique_output_file_175_180bp}")

    # Write non-unique sequences to a multifasta file and count the number of sequences
    non_unique_output_file = os.path.join(output_dir, f"{target_sequence}_non_unique_sequences.fasta")
    with open(non_unique_output_file, "w") as file:
        count = 0
        for i, sequence in enumerate(non_unique_sequences):
            count += 1
            header = f"{target_sequence} | Length: {len(sequence.seq)}"
            file.write(f">{header}\n{sequence.seq}\n")
        print(f"Number of non-unique sequences: {count}")
        print(f"Non-unique sequences saved in: {non_unique_output_file}")

    # Write the results to a CSV file
    with open(os.path.join(output_dir, f'{target_sequence}_counts.csv'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Filename", "Target Sequence", "Number of Variants Reconstructed", "Number of Variants Not Found",
                         "Number of Unique Sequences", "Number of Unique Sequences between 175 and 180 bp", "Number of Non-Unique Sequences"])
        writer.writerows(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process multi-FASTA files to search and reorganize sequences.')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing multi-FASTA files')
    parser.add_argument('-o', '--output', required=True, help='Output directory for results')
    parser.add_argument('-t', '--target', required=True, help='Target sequence to search in the multi-FASTA files')

    args = parser.parse_args()

    process_sequences(args.input, args.output, args.target)
