#!/usr/bin/env python
import os
import argparse
import csv
import glob
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp

def search_and_reorganize_sequences(args):
    fasta_path, _ = args  # No need for output_dir since we aren't using individual CSVs
    strain_name = os.path.basename(os.path.dirname(fasta_path))
    print(f"Processing strain: {strain_name}")

    sequence_dict = defaultdict(int)
    total_reads = 0
    reads_below_170 = 0
    reads_above_170 = 0
    unique_reads_over_100 = 0
    unique_reads_over_50 = 0
    unique_reads_over_20 = 0
    unique_reads_over_5 = 0
    single_occurrence = 0  # To count sequences that occur only once

    for record in SeqIO.parse(fasta_path, "fasta"):
        sequence = str(record.seq)
        sequence_dict[sequence] += 1
        total_reads += 1

    for sequence, count in sequence_dict.items():
        if len(sequence) < 170:
            reads_below_170 += 1
        elif len(sequence) >= 170:
            reads_above_170 += 1
            if count > 100:
                unique_reads_over_100 += 1
            elif count > 50:
                unique_reads_over_50 += 1
            elif count > 20:
                unique_reads_over_20 += 1
            elif count > 5:
                unique_reads_over_5 += 1
            elif count == 1:
                single_occurrence += 1
            
    return [strain_name, total_reads, reads_below_170, reads_above_170, single_occurrence, unique_reads_over_5, unique_reads_over_20, unique_reads_over_50, unique_reads_over_100]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate stats from FASTA TBR reads.')
    parser.add_argument('-o', '--output', help='The output directory', required=True)
    parser.add_argument('-d', '--directory', help='The directory with fasta files', required=True)
    
    args = parser.parse_args()
    print(f"Looking for fasta files in directory: {args.directory}")

    fasta_files = glob.glob(os.path.join(args.directory, '*', '*.TBR.bam.fasta'))

    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(search_and_reorganize_sequences, [(fasta_file, args.output) for fasta_file in fasta_files])

    # Write aggregated results to a single CSV file
    aggregate_stats_file = os.path.join(args.output, 'aggregate_stats.csv')
    with open(aggregate_stats_file, 'w') as f:
        print(f"Writing aggregated stats to file: {aggregate_stats_file}")
        writer = csv.writer(f)
        writer.writerow(['Strain', 'Total Reads', 'Total < 170', 'Total > 170', 'CN = 1', 'CN > 5', 'CN > 20', 'CN > 50', 'CN > 100'])
        writer.writerows(results)
