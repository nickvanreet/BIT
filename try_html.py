import os
from bs4 import BeautifulSoup

class TRF_run_parse():
    def __init__(self, fasta_file, trf_path, output_file):
        self.fasta_file = fasta_file
        self.trf_path = trf_path
        self.output_file = output_file

    def run(self):
        cmd = f'{self.trf_path} {self.fasta_file} 2 7 7 80 10 50 5000'
        os.system(cmd)

        name_for_trf = f'{self.fasta_file}.2.7.7.80.10.50.5000.1.txt.html'
        out_fasta = self.output_file

        with open(name_for_trf) as html_file, open(out_fasta, 'w') as out:
            soup = BeautifulSoup(html_file, 'html.parser')
            sequences = soup.find_all('pre')

            tr_number = 0
            q = 0

            for seq in sequences:
                if len(seq.contents) >= 3:  # Check if seq has at least 3 contents
                    seq_id = seq.contents[0].strip()
                    seq_seq = seq.contents[2].strip()

                    tr_number += 1
                    s = str(tr_number)

                    out.write('>' + seq_id + '|' + s + '\n' + seq_seq + '\n')

                    q += 1

        print('Total number of sequences found by TRF:', str(tr_number))
        print('Number of sequences in the output file:', str(q))


# Create an instance of TRF_run_parse
trf_parser = TRF_run_parse('test.fasta', '/usr/local/bin/trf', 'output.fasta')

# Run TRF and process the results
trf_parser.run()
