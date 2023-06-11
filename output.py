output_filename = "output.fasta"

with open("alignment_results.txt", "r") as blast_file, open(output_filename, "w") as fasta_file:
    for line in blast_file:
        if line.startswith("#"):
            continue  # Skip header lines
            
        fields = line.strip().split("\t")
        subject_start = int(fields[2])
        subject_end = int(fields[3])
        subject_sequence = fields[5]

        # Generate the FASTA header with start and end information
        fasta_header = f">{subject_start}_end{subject_end}"

        # Write the FASTA header and sequence to the output file
        fasta_file.write(fasta_header + "\n")
        fasta_file.write(subject_sequence + "\n")

print("FASTA files generated successfully!")