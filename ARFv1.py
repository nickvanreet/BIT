import os
import sys
import re

def run_trf(fasta_file, trf_path, output_file):
    cmd = f'{trf_path} {fasta_file} 2 7 7 80 10 50 500 -f -h'
    os.system(cmd)

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

    print('Total number of sequences found by TRF:', str(tr_number))


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
            for seq_id, sequence in repeat_sequences:
                consensus.write(f'>{seq_id}\n{sequence}\n')

        print("Repeat sequences with a length between 150 and 200 bp were saved to the consensus file.")


def create_blast_db(fasta_file, db_name):
    cmd = f'makeblastdb -in {fasta_file} -dbtype nucl -parse_seqids -out {db_name}'
    os.system(cmd)

    print(f"BLAST database '{db_name}' created.")


def run_blast(query_file, db_name, output_file):
    cmd = f'blastn -query {query_file} -db {db_name} -out {output_file} -outfmt "7 qstart qend sstart send qseq sseq"'
    print("Executing command:", cmd)  # Added for debugging
    os.system(cmd)

    print("BLAST search completed. Results saved to alignment_results.txt")


def extract_hits_from_alignment_results(alignment_file):
    hits = []
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

    return hits

def save_alignment_results_as_fasta(alignment_file, output_filename):
    with open(alignment_file, "r") as blast_file, open(output_filename, "w") as fasta_file:
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
            fasta_file.write(fasta_header + "\n")
            fasta_file.write(subject_sequence + "\n")

    print("FASTA file generated successfully!")

# Check if the input file argument is provided
if len(sys.argv) < 2:
    print("Please provide the input FASTA file as an argument.")
    sys.exit(1)

# Get the input FASTA file from the command-line argument
fasta_file = sys.argv[1]
fasta_filename = os.path.basename(fasta_file)

# Run TRF and process the results
trf_output_file = f'output_{fasta_filename}.fasta'
run_trf(fasta_file, '/usr/local/bin/trf', trf_output_file)

# Extract and save repeat sequences
output_file = f'consensus_{fasta_filename}'
extract_repeat_sequences(trf_output_file, output_file)

# Create BLAST database
db_name = fasta_filename.replace('.fasta', '')
create_blast_db(fasta_file, db_name)

# Run BLAST search
blast_output_file = f'alignment_results_{fasta_filename}.txt'
run_blast(output_file, db_name, blast_output_file)

# Save alignment results as a multi-FASTA file
output_filename = f"variants_{fasta_filename}"
save_alignment_results_as_fasta(blast_output_file, output_filename)