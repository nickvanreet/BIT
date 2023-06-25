from collections import defaultdict
from Bio import SeqIO
import random
import os

# Define the directory where you want to save the output files
output_dir = '/home/nick/ITOL'

# If the directory doesn't exist, create it
if not os.path.exists(os.path.expanduser(output_dir)):
    os.makedirs(os.path.expanduser(output_dir))

# Define a dictionary to store your data
data = defaultdict(lambda: defaultdict(int))

# Parse your fasta file
file_path = '/home/nick/RECO/GCGCAGTT_unique_sequences_175_180bp.fasta'
for record in SeqIO.parse(file_path, 'fasta'):
    # Get the unique sequence ID, count, and hits from the line
    fields = record.description.split(' | ')
    unique_id = fields[0][1:]
    count = int(fields[1][7:])
    
    # Get sequence length
    seq_length = len(record.seq)
    label = f'{seq_length}bp'
    data[unique_id][label] += 1

# Define labels and colors for multibar annotation
labels = sorted(set(label for unique_hits in data.values() for label in unique_hits.keys()))
colors = ['#%06x' % random.randint(0, 0xFFFFFF) for _ in range(len(labels))]

# Write multibar annotation to output file
with open(os.path.join(os.path.expanduser(output_dir),'itol_multibar_sequence_lengths.txt'), 'w') as f:
    # Write the header
    f.write('DATASET_MULTIBAR\nSEPARATOR TAB\nDATASET_LABEL\tmy_dataset\n')
    f.write('FIELD_LABELS\t' + '\t'.join(labels) + '\n')
    f.write('FIELD_COLORS\t' + '\t'.join(colors) + '\n')
    f.write('LEGEND_TITLE\tLegend\n')

    # Define lists for legend shapes, colors, and labels
    legend_shapes = ['1' for _ in labels]
    legend_colors = colors
    legend_labels = labels

    # Write legend information to the file
    f.write('LEGEND_SHAPES\t' + '\t'.join(legend_shapes) + '\n')
    f.write('LEGEND_COLORS\t' + '\t'.join(legend_colors) + '\n')
    f.write('LEGEND_LABELS\t' + '\t'.join(legend_labels) + '\n')

    f.write('DATA\n')

    # Write the data
    for unique_id, counts in data.items():
        modified_unique_id = 'Unique_' + unique_id.split('_')[1]
        values = []
        for label in labels:
            values.append(str(counts.get(label, 0)))
        f.write(modified_unique_id + '\t' + '\t'.join(values) + '\n')
