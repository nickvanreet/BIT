import os
import sys
import re
from Bio import SeqIO

def run_trf(fasta_file, trf_path, output_file):
    cmd = f'{trf_path} {fasta_file} 2 7 7 80 10 50 500 -f -h > /dev/null 2>&1'
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

    pattern = r'>.*\n.*'
    matches = re.findall(pattern, alignment_results, re.MULTILINE)

    for match in matches:
        lines = match.strip().split('\n')
        header = lines[0][1:]
        seq = lines[1]
        num_sequences += 1

        if seq != 'N/A':
            hits.append((header, seq))

    print("Total number of hits found:", num_sequences)

    return hits


def main(fasta_file, trf_path):
    # Create a directory for the target sequence
    target_sequence = os.path.splitext(os.path.basename(fasta_file))[0]
    output_dir = target_sequence + '_output'
    os.makedirs(output_dir, exist_ok=True)

    # Run TRF and process the results
    trf_output_file = os.path.join(output_dir, f'output_{target_sequence}.fasta')
    run_trf(fasta_file, trf_path, trf_output_file)

    # Extract consensus sequences
    consensus_output_file = os.path.join(output_dir, f'consensus_{target_sequence}.fasta')
    extract_repeat_sequences(trf_output_file, consensus_output_file)

    # Create a BLAST database
    blast_db_name = os.path.join(output_dir, target_sequence)
    create_blast_db(fasta_file, blast_db_name)

    # Run BLAST and process the results for variants
    blast_variants_output_file = os.path.join(output_dir, f'blast_variants_{target_sequence}.fasta')
    run_blast(consensus_output_file, blast_db_name, blast_variants_output_file)
    variants = extract_hits_from_alignment_results(blast_variants_output_file)

    # Run BLAST and process the results for recognized
    blast_recognized_output_file = os.path.join(output_dir, f'blast_recognized_{target_sequence}.fasta')
    run_blast(fasta_file, blast_db_name, blast_recognized_output_file)
    recognized = extract_hits_from_alignment_results(blast_recognized_output_file)

    # Run BLAST and process the results for not found
    blast_not_found_output_file = os.path.join(output_dir, f'blast_not_found_{target_sequence}.fasta')
    run_blast(trf_output_file, blast_db_name, blast_not_found_output_file)
    not_found = extract_hits_from_alignment_results(blast_not_found_output_file)

    print("\nVariants:")
    for header, seq in variants:
        print(f">{header}\n{seq}")
    
    print("\nRecognized:")
    for header, seq in recognized:
        print(f">{header}\n{seq}")
    
    print("\nNot Found:")
    for header, seq in not_found:
        print(f">{header}\n{seq}")

    print("\nProcess completed successfully!")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py fasta_file trf_path")
        sys.exit(1)

    fasta_file = sys.argv[1]
    trf_path = sys.argv[2]
    main(fasta_file, trf_path)
