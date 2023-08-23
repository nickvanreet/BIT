#!/usr/bin/env python
import os
import argparse
import subprocess

def run_mafft(input_file, output_file, num_threads):
    # Run MAFFT command
    mafft_cmd = f"mafft --auto --leavegappyregion --thread {num_threads} {input_file} > {output_file}"
    subprocess.run(mafft_cmd, shell=True, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align a _unique.fasta file using MAFFT")
    parser.add_argument("-i", "--input-file", help="Path to the input fasta file", required=True)
    parser.add_argument("-o", "--output-file", help="Path to the output alignment file", required=True)
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use for MAFFT")
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    num_threads = args.threads

    print(f"Aligning {input_file} file.")
    run_mafft(input_file, output_file, num_threads)
    print("MAFFT alignment completed.")
