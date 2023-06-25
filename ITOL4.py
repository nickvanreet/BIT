from collections import defaultdict
from Bio import SeqIO
import random
import math

# Define a dictionary to store your data
data = defaultdict(lambda: defaultdict(int))

# Parse your fasta file
for record in SeqIO.parse('short.fasta', 'fasta'):
    # Get the unique sequence ID, count, and hits from the line
    fields = record.description.split(' | ')
    unique_id = fields[0][1:]
    count = int(fields[1][7:])
    all_hits = fields[2][6:].split(' & ')

    for hit_string in all_hits:
        hit_parts = hit_string.split('_')
        subspecies = hit_parts[0]
        strain = hit_parts[1]
        clone = "_".join(hit_parts[2:])
        label = f"{subspecies}_{strain}"
        data[unique_id][label] += 1

# Define labels and colors
labels = set()
colors = []

# Add combined subspecies and strain labels and colors
for unique_hits in data.values():
    for label in unique_hits.keys():
        labels.add(label)

labels = sorted(labels)
colors = ['#%06x' % random.randint(0, 0xFFFFFF) for _ in range(len(labels))]

# Write to output file
with open('itol_multibar.txt', 'w') as f:
    # Write the header
    f.write('DATASET_MULTIBAR\nSEPARATOR TAB\nDATASET_LABEL\tmy_dataset\n')
    f.write('FIELD_LABELS\t' + '\t'.join(labels) + '\n')
    f.write('FIELD_COLORS\t' + '\t'.join(colors) + '\n')
    f.write('DATA\n')
    
    # Write the data
    for unique_id, counts in data.items():
        modified_unique_id = 'Unique_' + unique_id.split('_')[1]
        values = []
        for label in labels:
            values.append(str(counts.get(label, 0)))
        f.write(modified_unique_id + '\t' + '\t'.join(values) + '\n')
