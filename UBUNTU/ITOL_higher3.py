from collections import defaultdict
from Bio import SeqIO
import os

# Define the directory where you want to save the output files
output_dir = '/home/nick/ITOL'

# If the directory doesn't exist, create it
if not os.path.exists(os.path.expanduser(output_dir)):
    os.makedirs(os.path.expanduser(output_dir))

# Parse your fasta file
file_path = '/home/nick/RECO/GCGCAGTT_unique_sequences_175_180bp.fasta'
with open(os.path.join(os.path.expanduser(output_dir),'itol_text_labels.txt'), 'w') as f:
    # Write the header
    f.write('DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\ttest_labels\nCOLOR\t#000000\nDATA\n')
    
    for record in SeqIO.parse(file_path, 'fasta'):
        # Get the unique sequence ID and count from the line
        fields = record.description.split(' | ')
        unique_id = fields[0][1:]
        modified_unique_id = 'Unique_' + unique_id.split('_')[1]
        
        # Write the label for this sequence
        f.write(f'{modified_unique_id}\ttest\n')
