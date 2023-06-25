import random
from Bio import SeqIO
from collections import defaultdict

def get_random_color():
    r = lambda: random.randint(0,255)
    return('#%02X%02X%02X' % (r(),r(),r()))

# Predefined colors for subspecies
subspecies_colors = {"Tbb": "#FF0000", "Tbg": "#00FF00", "Tbr": "#0000FF"}
strain_colors = defaultdict(get_random_color)

# Dictionaries to store unique color combinations for each unique ID
subspecies_colors_combinations = defaultdict(set)
strain_colors_combinations = defaultdict(set)

# Parse fasta file
for record in SeqIO.parse('GCGCAGTT_unique_sequences_175_180bp.fasta', 'fasta'):
    # Split header by ' | '
    fields = record.description.split(' | ')
    
    # Get unique id
    unique_id = fields[0]
    
    # Get hits and split by ' & ' to get a list of all hits
    all_hits = fields[2][6:].split(' & ')

    for hits in all_hits:
        # Split the hits into subspecies, strain, and clone
        hits_parts = hits.split('_')
        subspecies = hits_parts[0]
        strain = hits_parts[1]
        clone = "_".join(hits_parts[2:])

        # Store the color combinations
        subspecies_colors_combinations[unique_id].add(subspecies_colors[subspecies])
        strain_colors_combinations[unique_id].add(strain_colors[strain])

# Write to subspecies and strain files
with open('itol_subspecies.txt', 'w') as outfile:
    outfile.write('TREE_COLORS\nSEPARATOR TAB\nDATA\n')
    for unique_id, colors in subspecies_colors_combinations.items():
        for color in colors:
            outfile.write(f'{unique_id}\trange\t{color}\t1\tsubspecies\n')

with open('itol_strain.txt', 'w') as outfile:
    outfile.write('TREE_COLORS\nSEPARATOR TAB\nDATA\n')
    for unique_id, colors in strain_colors_combinations.items():
        for color in colors:
            outfile.write(f'{unique_id}\trange\t{color}\t1\tstrain\n')
