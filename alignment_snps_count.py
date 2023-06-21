from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def parse_sequence_file(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequence_info = record.description.split(' | ')
        sequence = str(record.seq)
        sequences.append((sequence_info, sequence))
    return sequences

def find_unique_tbg_hits(sequences):
    unique_tbg_hits = {}
    for sequence_info, sequence in sequences:
        hits = sequence_info[2].replace('Hits: ', '').split(' & ')
        for hit in hits:
            if hit.startswith('Tbb'):
                if hit in unique_tbg_hits:
                    unique_tbg_hits[hit] = unique_tbg_hits.get(hit, 0) + 1
                else:
                    unique_tbg_hits[hit] = 1
    return unique_tbg_hits

def find_snps(sequences):
    snps = []
    for sequence_info, sequence in sequences:
        hits = sequence_info[2].replace('Hits: ', '').split(' & ')
        for hit in hits:
            seq = sequences[hits.index(hit)][1]
            alignments = pairwise2.align.globalxx(sequence, seq)
            for alignment in alignments:
                aligned_seq1, aligned_seq2, _, _, _ = alignment
                snp_positions = [i for i, (a, b) in enumerate(zip(aligned_seq1, aligned_seq2)) if a != b]
                for position in snp_positions:
                    snp = f"{hit} - {aligned_seq1[position]}-{aligned_seq2[position]} at {position}"
                    snps.append(snp)
    return snps

# Path to your sequence file
file_path = 'unique_sequences_GCGCAGTT.fasta'

# Parse the sequence file
sequences = parse_sequence_file(file_path)

# Find unique TBG hits and their counts
unique_tbg_hits = find_unique_tbg_hits(sequences)

# Find SNPs in aligned sequences
snps = find_snps(sequences)

# Print the unique TBG hits and their counts
print("Unique TBG Hits:")
for hit, count in unique_tbg_hits.items():
    print(f"{hit}: {count}")

# Print the SNPs
print("\nSNPs:")
for snp in snps:
    print(snp)
