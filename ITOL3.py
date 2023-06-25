from collections import defaultdict
from Bio import SeqIO
import random

# Define your labels and colors
subspecies_labels = ['Tbb', 'Tbg', 'Tbr']
subspecies_colors = ['#ff0000', '#00ff00', '#0000ff']

# Define a dictionary to store your data
data = defaultdict(lambda: defaultdict(int))

# Parse your fasta file
for record in SeqIO.parse('GCGCAGTT_unique_sequences_175_180bp.fasta', 'fasta'):
    # Get the unique sequence ID, count, and hits from the line
    fields = record.description.split(' | ')
    unique_id = fields[0]
    count = int(fields[1][7:])
    all_hits = fields[2][6:].split(' & ')

    for hit_string in all_hits:
        hit_parts = hit_string.split('_')
        subspecies = hit_parts[0]
        strain = hit_parts[1]
        clone = "_".join(hit_parts[2:])
        data[unique_id][subspecies] += count

# Write to output file
with open('itol_multibar.txt', 'w') as f:
    # Write the header
    f.write('DATASET_MULTIBAR\nSEPARATOR TAB\nDATASET_LABEL\tmy_dataset\nCOLOR\t#ff0000\n')
    f.write('FIELD_LABELS\t' + '\t'.join(subspecies_labels) + '\n')
    f.write('FIELD_COLORS\t' + '\t'.join(subspecies_colors) + '\n')
    f.write('DATA\n')
    
    # Write the data
    for unique_id, counts in data.items():
        f.write(unique_id + '\t' + '\t'.join(str(counts[label]) for label in subspecies_labels) + '\n')
