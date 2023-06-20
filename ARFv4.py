import argparse
import os
import csv
import re
import subprocess
import shutil
from collections import defaultdict

from Bio import SeqIO

def run_trf(fasta_file, trf_path, name_for_trf, trf_output_file):
    cmd = f'{trf_path} {fasta_file} 2 7 7 80 10 50 500 -d -h'
    subprocess.run(cmd, shell=True)

    counts = {'<50': 0, '50-100': 0, '100-250': 0, '250-500': 0, '>500': 0}
    Rn_values = {}
    rep_number = 0

    # Modify this line to include the output directory in the TRF output file path
    with open(name_for_trf, 'r') as fasta, open(trf_output_file, 'a') as out:
        tr_number = 0
        seq_id = None
        seq_seq = ''
        a = [str(line) for line in fasta]

        for lin in a:
            if lin[0] == 'S':
                seq_id = lin.split(' ')[1].rstrip()
            elif lin[0].isdigit():
                split_lin = lin.split(' ')
                tr_number += 1
                seq_seq = split_lin[13]
                s = str(tr_number)
                Rn = split_lin[3]
                Rn_values[seq_id] = Rn

                out.write('>' + seq_id + '|' + s + '|' + str(Rn) + '|' + str(len(seq_seq)) + '\n' + seq_seq + '\n')

            repeat_length = len(seq_seq)

            if repeat_length < 50:
                counts['<50'] += 1
            elif 50 <= repeat_length < 100:
                counts['50-100'] += 1
            elif 100 <= repeat_length < 250:
                counts['100-250'] += 1
            elif 250 <= repeat_length < 500:
                counts['250-500'] += 1
            else:
                counts['>500'] += 1

    print('\nTotal number of repeats found by TRF:', str(tr_number))
    shutil.move(name_for_trf, os.path.dirname(trf_output_file))
    
    return tr_number, counts, Rn_values

def extract_repeat_sequences(trf_output_file, output_file_consensus):
    repeat_sequences = []

    with open(trf_output_file) as trf_file:
        seq_id = None
        for line in trf_file:
            if line.startswith('>'):
                seq_id = line.strip().lstrip('>')
            elif line.startswith('Parameters:'):
                continue
            else:
                sequence = line.strip()
                repeat_length = len(sequence)

                if 100 <= repeat_length <= 250:
                    repeat_sequences.append((seq_id, sequence))

    if not repeat_sequences:
        print("No repeat sequences with a length between 150 and 200 bp were found.")
        return {"rep_number": 0, "repeat_sequences": []}

    else:
        with open(output_file_consensus, 'w') as consensus:
            rep_number = 0
            for seq_id, sequence in repeat_sequences:
                consensus.write(f'>{seq_id}\n{sequence}\n')
                rep_number += 1

        print("Total number of consensus sequences found:", str(rep_number))
        return {"rep_number": rep_number, "repeat_sequences": repeat_sequences}

def create_blast_db(fasta_file, db_name):
    cmd = f'makeblastdb -in {fasta_file} -dbtype nucl -parse_seqids -out {db_name}'
    subprocess.run(cmd, shell=True)

def run_blast(query_file, db_name, output_file_consensus):
    cmd = f'blastn -query {query_file} -db {db_name} -out {output_file_consensus} -outfmt "7 qstart qend sstart send qseq sseq"'
    subprocess.run(cmd, shell=True)

def extract_hits_from_alignment_results(alignment_file):
    hits = []
    num_sequences = 0

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
    num_variants = 0
    variant_headers = []

    with open(alignment_file, "r") as blast_file, open(variant_output_filename, "w") as variant_file:
        for line in blast_file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            subject_start = int(fields[2])
            subject_end = int(fields[3])
            subject_sequence = fields[5]

            fasta_header = f">{fasta_filename}_{subject_start}_{subject_end}"
            variant_headers.append(fasta_header)

            variant_file.write(fasta_header + "\n")
            variant_file.write(subject_sequence + "\n")

            num_variants += 1

    print("Total number of variant files extracted:", num_variants)
    return variant_headers, num_variants


def process_fasta_files(input_dir, trf_path, output_dir, min_sequence_length=200):
    ascii_art = r'''
    
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    \/////// \\//////// \\///////// \
    /// \\/// \/// \\/// \/// \\\\\\\
    ///////// \//////// \\/////// \\\
    /// \\/// \/// \\/// \/// \\\\\\\
    /// \\/// \/// \\/// \/// \\\\\\\
    \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    '''
    print(ascii_art)

    # Check if the output directory exists, and if not, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # process the FASTA files
    fasta_files = [filename for filename in os.listdir(input_dir) if filename.endswith(".fasta")]
    num_files = len(fasta_files)

    results = []
    print(f"{num_files} FASTA file(s) detected in {input_dir}.\n")
    print(f"ARF will save all output to {output_dir}.\n ")

    for index, filename in enumerate(fasta_files, start=1):
        fasta_file = os.path.join(input_dir, filename)

        print(f"Processing FASTA {index}/{num_files}: {filename}")

        with open(fasta_file, 'r') as file:
            fasta_data = file.read().strip().split('\n', 1)
            sequence_length = len(fasta_data[1].replace('\n', ''))

        print(f"Sequence length: {sequence_length} bases")

        tr_number = 0

        if sequence_length < min_sequence_length:
            print(f"Sequence is shorter than {min_sequence_length} bases. Skipping to the next sequence.")
            skipped_filename = f"variant_{filename}"
            skipped_file_path = os.path.join(output_dir, skipped_filename)
            with open(skipped_file_path, 'w') as skipped_file:
                skipped_file.write('\n'.join(fasta_data)) # fasta_data contains both header and sequence

            results.append({'fasta_file': filename, 'sequence_length': sequence_length, 'tr_number': 'Skipped', 'rep_number': 1, 'num_variants': 0})
            continue

        else:
            name_for_trf = os.path.join(input_dir, f'{filename}.2.7.7.80.10.50.500.dat')
            trf_output_file = os.path.join(output_dir, f'output_{filename}.fasta')
            tr_number, counts, Rn_values = run_trf(fasta_file, trf_path, name_for_trf, trf_output_file)

            output_file_consensus = os.path.join(output_dir, f'consensus_{filename}')
            consensus_data = extract_repeat_sequences(trf_output_file, output_file_consensus)
            num_consensus_sequences = consensus_data['rep_number']

            if num_consensus_sequences == 0:
                print(f"No sequences between 150-200 bp found in {filename}. Saving the original sequence as a variant.")
                variant_output_filename = os.path.join(output_dir, f"variants_{filename}")
                with open(fasta_file, "r") as original, open(variant_output_filename, "w") as variant:
                    variant.write(original.read())
                    results.append({'fasta_file': filename, 'sequence_length': sequence_length, 'tr_number': tr_number, 'rep_number': 0, 'num_variants': 1})
                continue
            
            extract_repeat_sequences(trf_output_file, output_file_consensus)

            db_name = os.path.join(output_dir, filename)
            create_blast_db(fasta_file, db_name)

            blast_output_file = os.path.join(output_dir, f'alignment_results_{filename}.txt')
            run_blast(output_file_consensus, db_name, blast_output_file)

            variant_output_filename = os.path.join(output_dir, f"variants_{filename}")
            variant_headers, num_variants = save_alignment_results_as_fasta(blast_output_file, variant_output_filename, filename)

            hits, num_sequences = extract_hits_from_alignment_results(blast_output_file)

            result = {'fasta_file': filename, 'sequence_length': sequence_length, 'tr_number': tr_number, 'Rn': Rn_values, 'rep_number': num_consensus_sequences, 'num_sequences': num_sequences, 'num_variants': num_variants}
            result.update(counts)
            results.append(result)

    with open(os.path.join(output_dir, 'results.csv'), 'w', newline='') as csvfile:
        fieldnames = ['fasta_file', 'sequence_length', 'tr_number', '<50', '50-100', '100-250', '250-500', '>500', 'Rn', 'rep_number', 'num_sequences', 'num_variants']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in results:
            writer.writerow(result)

parser = argparse.ArgumentParser(description="Process some FASTA files.")
parser.add_argument('InputDirectory', type=str, help='The path to the directory that contains FASTA files')
parser.add_argument('TRFPath', type=str, help='The path to the Tandem Repeats Finder executable')
parser.add_argument('OutputDirectory', type=str, help='The path to the directory to store output files')
args = parser.parse_args()

process_fasta_files(args.InputDirectory, args.TRFPath, args.OutputDirectory)
