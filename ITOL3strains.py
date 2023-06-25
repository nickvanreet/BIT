from collections import defaultdict
from Bio import SeqIO
import random
import math

def get_random_color():
    r = lambda: random.randint(0,255)
    return('#%02X%02X%02X' % (r(),r(),r()))

# Define a dictionary to store your data
data = defaultdict(lambda: defaultdict(int))

# Strain names provided
strain_names = ['AnTat2.2', 'AnTat5.2', 'AnTat11.17', 'J10', 'Lister427', 'TREU927', 'TSW196', '348BT', 'AnTat9.1', 'AnTat17.1', 'AnTat22.2', 'DAL972', 'KOBIR', 'MM01', 'OUSOU', 'RUMPHI']

# Assign colors to strains
strain_colors = {name: get_random_color() for name in strain_names}

# Parse your fasta file
for record in SeqIO.parse('GCGCAGTT_unique_sequences_175_180bp.fasta', 'fasta'):
    # Get the unique ID from the record ID
    unique_id = record.id.strip()

    # Get the count from the record description
    count = int(record.description.split(' | ')[1][7:])

    # Get the hits from the record description
    all_hits = record.description.split(' | ')[2][6:].split(' & ')

    for hit_string in all_hits:
        hit_parts = hit_string.split('_')
        strain = hit_parts[1]
        if strain in strain_names:
            data[unique_id][strain] += math.log(count, 2)

# Write to output file
with open('itol_strain_multibar.txt', 'w') as f:
    # Write the header
    f.write('DATASET_MULTIBAR\nSEPARATOR TAB\nDATASET_LABEL\tmy_dataset\nCOLOR\t#ff0000\n')
    f.write('FIELD_LABELS\t' + '\t'.join(strain_names) + '\n')
    f.write('FIELD_COLORS\t' + '\t'.join(strain_colors[name] for name in strain_names) + '\n')
    f.write('DATA\n')
    
    # Write the data
    for unique_id, counts in data.items():
        f.write(unique_id + '\t' + '\t'.join(str(counts[name]) for name in strain_names) + '\n')
