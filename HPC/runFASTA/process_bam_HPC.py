#!/usr/bin/env python
import os
import pysam
from Bio.Seq import Seq
import glob
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="Directory containing BAM files")
parser.add_argument("output_dir", help="Directory to save output FASTA files")
args = parser.parse_args()

# Input and output directories
input_dir = args.input_dir
output_dir = args.output_dir

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Input validation: Check if input directory exists
if not os.path.exists(input_dir):
    print(f"Error: Input directory {input_dir} does not exist.")
    exit(1)

# Get list of BAM files
bam_files = glob.glob(os.path.join(input_dir, '*.TBR.bam'))

# Error handling: Check if bam_files list is empty
if not bam_files:
    print(f"No BAM files found in directory {input_dir}")
    exit(1)

for bam_file in bam_files:
    print(f"Processing file {bam_file}")
    
    # Define output file path
    output_file = os.path.join(output_dir, os.path.basename(bam_file) + '.fasta')
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as infile, open(output_file, "w") as outfile:
            # Iterate over each read in the BAM file
            for read in infile:
                # Get the mapped sequence from the read
                mapped_sequence = read.query_sequence[read.query_alignment_start:read.query_alignment_end]

                # Always treat the sequence as if it were aligned to the forward strand
                sequence = mapped_sequence

                # Padding the sequence
                sequence = sequence + 'X' * 27
                outfile.write(">" + read.query_name + "\n" + sequence + "\n")
    except IOError as e:
        print(f"Error opening/processing file {bam_file}: {e}")
        continue

    print(f"Finished processing file {bam_file}")
