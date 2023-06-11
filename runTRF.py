## modified by Nick Van Reet from TRF_run_parse.py
import os

def run_trf(fasta_file, trf_path, output_file):
    cmd = f'{trf_path} {fasta_file} 2 7 7 80 10 50 500 -f -h '
    os.system(cmd)

    name_for_trf = f'{fasta_file}.2.7.7.80.10.50.500.dat'

    with open(name_for_trf) as fasta, open(output_file, 'a') as out:
        tr_number = 0
        seq_id = []
        q = 0
        a = [str(line) for line in fasta]

        for lin in a:
            # identification of start of the sequence
            # name of the sequence
            if lin[0] == 'S':
                seq_id = lin.split(' ')[1].rstrip()

            # identification of the sequence
            ## if line starts with an integer then it corresponds to the TR line
            check = lin.split(' ')[0]
            if check.isdigit() == True:
                split_lin = lin.split(' ')
                tr_number += 1
                seq_seq = split_lin[13]
                s = str(tr_number)
                Rn = split_lin[3]

                out.write('>' + seq_id + '|' + s + '|' + str(Rn) + '|' + str(len(seq_seq)) + '\n' + seq_seq + '\n')

    print('Total number of sequences found by TRF:', str(tr_number))


# Get all *.fasta files in the directory
fasta_files = [filename for filename in os.listdir('.') if filename.endswith('.fasta')]

for fasta_file in fasta_files:
    # Run TRF and process the results
    run_trf(fasta_file, '/usr/local/bin/trf', f'output.{fasta_file}.fasta')
