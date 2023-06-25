from Bio import SeqIO

# Open output file
with open('itol_metadata.txt', 'w') as outfile:
    # Write header
    outfile.write('DATASET_SIMPLEBAR\nSEPARATOR TAB\nDATASET_LABEL\tmy_dataset\nDATA\n')
    
    # Parse fasta file
    for record in SeqIO.parse('GCGCAGTT_unique_sequences_175_180bp.fasta', 'fasta'):
        # Split header by ' | '
        fields = record.description.split(' | ')
        
        # Get unique id
        unique_id = fields[0]
        
        # Get count (removing leading ' Count: ')
        count = fields[1][7:]
        
        # Get hits and split by ' & ' to get a list of all hits
        all_hits = fields[2][6:].split(' & ')

        for hits in all_hits:
            # Split the hits into subspecies, strain, and clone
            hits_parts = hits.split('_')
            subspecies = hits_parts[0]
            strain = hits_parts[1]
            clone = "_".join(hits_parts[2:])
            
            # Write to output file
            outfile.write(f'{unique_id}\t{count}\t{subspecies}\t{strain}\n')
