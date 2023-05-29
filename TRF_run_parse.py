import os

class TRF_run_parse():
    def __init__(self, fasta_file, trf_path, output_file):
        self.fasta_file = fasta_file
        self.trf_path = trf_path
        self.output_file = output_file

    def run(self):
        cmd = f'{self.trf_path} {self.fasta_file} 2 7 7 80 10 50 5000 -f -h -m'
        os.system(cmd)

        name_for_trf = f'{self.fasta_file}.2.7.7.80.10.50.5000.dat'
        out_fasta = self.output_file

        with open(name_for_trf) as fasta, open(out_fasta, 'w') as out:
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
                ## if line starts with integer then it corresponds to TR line
                check = lin.split(' ')[0]
                if check.isdigit() == True:
                    split_lin = lin.split(' ')
                    tr_number += 1
                    seq_seq = split_lin[13]
                    s = str(tr_number)
                    Rn = split_lin[3]

                    out.write('>' + seq_id + '|' + s + '|' + str(Rn) + '|' + str(len(seq_seq)) + '\n' + seq_seq + '\n')

        print('Total number of sequences found by TRF:', str(tr_number))
        print('Number of sequences in the output file:', str(q))


# Create an instance of TRF_run_parse
trf_parser = TRF_run_parse('test.fasta', '/usr/local/bin/trf', 'output.fasta')

# Run TRF and process the results
trf_parser.run()
