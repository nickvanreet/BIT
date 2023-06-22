import pandas as pd
import os
import glob

# List all the target_sequence_counts files
files = glob.glob('*_counts.csv')

# Sort files to get the first file and remaining files
files.sort()
first_file = files[0]
remaining_files = files[1:]

# Read in the first file
df = pd.read_csv(first_file)

# Read in the remaining files (without headers)
for file in remaining_files:
    temp_df = pd.read_csv(file, header=None, skiprows=1)
    df = pd.concat([df, temp_df])

# Save the concatenated DataFrame to a new csv file
df.to_csv("combined_counts.csv", index=False)
