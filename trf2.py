import os
import subprocess

def run_trf(sequence_file):
    trf_path = "/usr/local/bin/trf"  # Path to the TRF executable

    # Generate the output file name
    output_file = f"{sequence_file}.trf_output.dat"

    # Run TRF on the sequence file and generate the output file
    command = f"{trf_path} {sequence_file} 2 5 7 80 10 50 2000 -d -h > {output_file}"
    subprocess.run(command, shell=True)

    return output_file

# Example usage
fasta_file = "test.fasta"
trf_output_file = run_trf(fasta_file)
print(f"TRF output file generated: {trf_output_file}")