from collections import defaultdict
from Bio import SeqIO

# Define a dictionary to store your data
data = defaultdict(int)

# Parse your fasta file
for record in SeqIO.parse('short.fasta', 'fasta'):
    # Get the unique sequence ID, count, and hits from the line
    fields = record.description.split(' | ')
    unique_id = fields[0][1:]
    
    # Extract the count value
    count = 0
    count_field = next((field for field in fields if field.startswith('Count:')), None)
    if count_field:
        try:
            count = int(count_field.split(':')[1])
        except ValueError:
            pass

    all_hits = fields[2][6:].split(' & ')

    # Check if all hits are Tbg
    if all(hit.startswith('Tbg') for hit in all_hits):
        data[unique_id] = count

# Write to output file
with open('itol_multibar.txt', 'w') as f:
    # Write the header
    f.write('DATASET_MULTIBAR\nSEPARATOR TAB\nDATASET_LABEL\tmy_dataset\n')
    f.write('FIELD_LABELS\tTbg\n')
    f.write('FIELD_COLORS\t#ff0000\n')
    f.write('DATA\n')
    
    # Write the data
    for unique_id, count in data.items():
        f.write(f'{unique_id}\t{count}\n')
