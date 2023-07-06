from Bio import AlignIO, Align
from Bio.Phylo.TreeConstruction import DistanceCalculator

import pandas as pd
import numpy as np

# load the alignment file
alignment = AlignIO.read('GCGCAGTT_aligned.fa', 'fasta')

# calculate the distance matrix
calculator = DistanceCalculator('identity')  # uses the simple identity model to calculate Hamming distance
dm = calculator.get_distance(alignment)

# Create a DataFrame from the distance matrix, including labels
df = pd.DataFrame(dm.matrix, columns = dm.names, index = dm.names)

# Save the DataFrame to a CSV file
df.to_csv("distance_matrix.csv")

