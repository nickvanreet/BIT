import os
import sys
import re
import subprocess
from Bio import SeqIO

def run_trf(fasta_file, trf_path, output_file):
    cmd = [trf_path, fasta_file, '2', '7', '7', '80', '10', '50', '500', '-f', '-h']
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    name_for_trf = f'{fasta_file}.2.7.7.80.10.50.500.dat'

    with open(name_for_trf) as fasta, open(output_file, 'a') as out:
        tr_number = 0
        seq_id = []
        a = [str(line) for line in fasta]
     
        for lin in a:
            # identification of start of the sequence
            # name of the sequence
            if lin[0] == 'S':
                seq_id = lin.split(' ')[1].rstrip()

            # identification of the sequence
            ## if line starts with an integer then it corresponds to the TR line
            check = lin.split(' ')[0]
            if check.isdigit():
                split_lin = lin.split(' ')
                tr_number += 1
                seq_seq = split_lin[13]
                s = str(tr_number)
                Rn = split_lin[3]

                out.write('>' + seq_id + '|' + s + '|' + str(Rn) + '|' + str(len(seq_seq)) + '\n' + seq_seq + '\n')

    print('\nTotal number of repeats found by TRF:', str(tr_number))


def extract_repeat_sequences(trf_output_file, output_file):
    repeat_sequences = []

    with open(trf_output_file) as trf_file:
        for line in trf_file:
            if line.startswith('>'):
                seq_id = line.strip().lstrip('>')
            elif line.startswith('Parameters:'):
                # Skip the parameters line
                continue
            else:
                sequence = line.strip()
                repeat_length = len(sequence)

                if repeat_length >= 150 and repeat_length <= 200:
                    repeat_sequences.append((seq_id, sequence))

    if not repeat_sequences:
        print("No repeat sequences with a length between 150 and 200 bp were found.")
    else:
        with open(output_file, 'w') as consensus:
            rep_number = 0
            for seq_id, sequence in repeat_sequences:
                consensus.write(f'>{seq_id}\n{sequence}\n')
                rep_number += 1

        print("Total number of consensus sequences found:", str(rep_number))


def create_blast_db(fasta_file, db_name):
    cmd = f'makeblastdb -in {fasta_file} -dbtype nucl -parse_seqids -out {db_name} > /dev/null 2>&1'
    os.system(cmd)

    # print(f"BLAST database '{db_name}' created.")


def run_blast(query_file, db_name, output_file):
    cmd = f'blastn -query {query_file} -db {db_name} -out {output_file} -outfmt "7 qstart qend sstart send qseq sseq" > /dev/null 2>&1'
    #print("Executing command:", cmd)  # Added for debugging
    os.system(cmd)

    


def extract_hits_from_alignment_results(alignment_file):
    hits = []
    num_sequences = 0  # Counter for the number of sequences

    with open(alignment_file, 'r') as file:
        alignment_results = file.read()

    pattern = r'>.*?Length=(\d+).*?\n(.*?)\n\n'
    matches = re.findall(pattern, alignment_results, re.DOTALL)
    for match in matches:
        length = int(match[0])
        subject_lines = match[1].split('\n')[1:-1]
        subject_sequence = ''.join([line.split()[-1] for line in subject_lines])
        start = int(re.search(r'\bstart=(\d+)', subject_lines[0]).group(1))
        end = int(re.search(r'\bend=(\d+)', subject_lines[-1]).group(1))
        hits.append({'start': start, 'end': end, 'sequence': subject_sequence})
        num_sequences += 1
    return hits, num_sequences

def save_alignment_results_as_fasta(alignment_file, variant_output_filename, fasta_filename):
    num_variants = 0  # Counter for the number of variants
    
    with open(alignment_file, "r") as blast_file, open(variant_output_filename, "w") as variant_file:
        for line in blast_file:
            if line.startswith("#"):
                continue  # Skip header lines

            fields = line.strip().split("\t")
            subject_start = int(fields[2])
            subject_end = int(fields[3])
            subject_sequence = fields[5]

            # Generate the FASTA header with start and end information
            fasta_header = f">{fasta_filename}_{subject_start}_{subject_end}"

            # Write the FASTA header and sequence to the output file
            variant_file.write(fasta_header + "\n")
            variant_file.write(subject_sequence + "\n")

            num_variants += 1

    print("Total number of variant files extracted:", num_variants)       

def search_and_reorganize_sequences(variant_output_filename, target_sequence, output_file, not_found_file):
    fasta_entries = SeqIO.parse(variant_output_filename, 'fasta')
    reconstructed_sequences = []
    sequences_not_found = []

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


def process_fasta_files(input_dir, target_sequence):
    # Create output directory based on the target sequence name
    output_dir = os.path.join(os.path.dirname(target_sequence), os.path.basename(target_sequence).split('.')[0] + "_output")
    os.makedirs(output_dir, exist_ok=True)
    ascii_art = r'''
      ____ __ __ ______  ___  ___ ___  ____ ______ ____   __      ____    ___ ____   ___  ____ ______      _____ ____ ____  ___     ___ ____  
     /    |  |  |      |/   \|   |   |/    |      |    | /  ]    |    \  /  _]    \ /  _]/    |      |    |     |    |    \|   \   /  _]    \ 
    |  o  |  |  |      |     | _   _ |  o  |      ||  | /  /     |  D  )/  [_|  o  )  [_|  o  |      |    |   __||  ||  _  |    \ /  [_|  D  )
    |     |  |  |_|  |_|  O  |  \_/  |     |_|  |_||  |/  /      |    /|    _]   _/    _]     |_|  |_|    |  |_  |  ||  |  |  D  |    _]    / 
    |  _  |  :  | |  | |     |   |   |  _  | |  |  |  /   \_     |    \|   [_|  | |   [_|  _  | |  |      |   _] |  ||  |  |     |   [_|    \ 
    |  |  |     | |  | |     |   |   |  |  | |  |  |  \     |    |  .  \     |  | |     |  |  | |  |      |  |   |  ||  |  |     |     |  .  \
    |__|__|\__,_| |__|  \___/|___|___|__|__| |__| |____\____|    |__|\_|_____|__| |_____|__|__| |__|      |__|  |____|__|__|_____|_____|__|\_|
    '''
    print(ascii_art)
    fasta_files = [filename for filename in os.listdir(input_dir) if filename.endswith(".fasta")]
    num_files = len(fasta_files)

    print(f"{num_files} FASTA file(s) detected in {input_dir}.")
    print(f"Target sequence requested: {target_sequence}\n")

    # Process each FASTA file in the input directory
    for index, filename in enumerate(fasta_files, start=1):
        fasta_file = os.path.join(input_dir, filename)
        fasta_filename = os.path.basename(fasta_file)

        print(f"Processing FASTA {index}/{num_files}: {fasta_filename}")

        # Read the FASTA file and get the sequence length
        with open(fasta_file, 'r') as file:
            fasta_data = file.read().strip().split('\n', 1)
            sequence_length = len(fasta_data[1].replace('\n', ''))

        print(f"Sequence length: {sequence_length} bases")

        # Run TRF and process the results
        trf_output_file = f'output_{fasta_filename}.fasta'
        run_trf(fasta_file, '/usr/local/bin/trf', trf_output_file)

        # Extract and save repeat sequences
        output_file = os.path.join(output_dir, f'consensus_{fasta_filename}')
        extract_repeat_sequences(trf_output_file, output_file)

        # Create BLAST database
        db_name = fasta_filename.replace('.fasta', '')
        create_blast_db(fasta_file, db_name)

        # Run BLAST search
        blast_output_file = os.path.join(input_dir, f'alignment_results_{fasta_filename}.txt')
        run_blast(output_file, db_name, blast_output_file)

        # Save alignment results as a multi-FASTA file
        variant_output_filename = os.path.join(output_dir, f"variants_{fasta_filename}")
        save_alignment_results_as_fasta(blast_output_file, variant_output_filename, fasta_file)

        # Reconstruct the variants multi-FASTA file
        output_file = os.path.join(output_dir, f'reorganized_{target_sequence}_{fasta_filename}')
        not_found_file = os.path.join(output_dir, f'sequences_not_found_{target_sequence}_{fasta_filename}')
        search_and_reorganize_sequences(variant_output_filename, target_sequence, output_file, not_found_file)

        # Remove intermediate files
        os.remove(f'{fasta_file}.2.7.7.80.10.50.500.dat')
        os.remove(trf_output_file)
        os.remove(blast_output_file)
        os.remove(f"{db_name}.nhr")
        os.remove(f"{db_name}.nin")
        os.remove(f"{db_name}.nsq")
        os.remove(f"{db_name}.ndb")
        os.remove(f"{db_name}.nog")
        os.remove(f"{db_name}.nos")
        os.remove(f"{db_name}.not")
        os.remove(f"{db_name}.ntf")
        os.remove(f"{db_name}.nto")

        print("\n")

    print("Processing of FASTA files in the directory is complete.\nOutput is saved in:", input_dir, output_dir)
   

# Check if the input directory argument is provided
if len(sys.argv) < 3:
    print("Please provide the input directory path and the target sequence as arguments.")
    sys.exit(1)

# Get the input directory and target sequence from the command-line arguments
input_directory = sys.argv[1]
target_sequence = sys.argv[2]

# Process the FASTA files in the input directory
process_fasta_files(input_directory, target_sequence)