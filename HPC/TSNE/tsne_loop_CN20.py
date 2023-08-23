import pandas as pd
from sklearn.manifold import TSNE

# Load the data
data = pd.read_csv('copy_number_lookup_CN20.csv')

# Normalize the copy numbers
strain_totals = data.groupby('Strain')['CopyNumber'].sum()
data['NormalizedCopyNumber'] = data.apply(lambda row: row['CopyNumber'] / strain_totals[row['Strain']], axis=1)

# Pivot the normalized data to wide format
data_wide = data.pivot(index='Strain', columns='Unique_ID', values='NormalizedCopyNumber').fillna(0)

# Check data shape
print(f"Data shape: {data_wide.shape}")

# List of perplexities and learning rates to test
perplexities = [10, 30, 50, 70]
learning_rates = [100, 200, 500, 1000]

# Iterate over different t-SNE parameters and compute results
for perplexity in perplexities:
    for learning_rate in learning_rates:
        # Inform about starting t-SNE
        print(f"Starting t-SNE computation with perplexity={perplexity} and learning_rate={learning_rate}...")

        # Run t-SNE
        tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity, n_iter=1000, learning_rate=learning_rate)
        tsne_results = tsne.fit_transform(data_wide)

        # Inform about finished t-SNE computation
        print("Finished t-SNE computation.")

        # Save the t-SNE results
        filename = f'tsne_results_CN5_perplexity_{perplexity}_lr_{learning_rate}.csv'
        df_tsne = pd.DataFrame(data={'Dim1': tsne_results[:,0], 'Dim2': tsne_results[:,1]})
        df_tsne['Strain'] = data_wide.index
        df_tsne.to_csv(filename, index=False)

        print(f"Results saved to {filename}")
