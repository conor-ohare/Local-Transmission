import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from sklearn.metrics import pairwise_distances

# Create 'plots' folder if it doesn't exist
plots_folder = "plots"
if not os.path.exists(plots_folder):
    os.makedirs(plots_folder)
    print(f"'{plots_folder}' folder created for saving plots.")

# Load metadata with columns: AccessionID, PostalCode, Serotype, Genotype, MajorLineage
metadata_file = "data/NEA_cleaned.metadata.csv"
metadata_df = pd.read_csv(metadata_file)

# Filter metadata to exclude entries with PostalCode equal to 0
metadata_df = metadata_df[metadata_df['PostalCode'] != 0]

# Load Singapore postal code data for coordinates
postal_data_file = "data/SG_postal.csv"
postal_df = pd.read_csv(postal_data_file)

# Ensure both PostalCode columns are strings to prevent merge issues
metadata_df['PostalCode'] = metadata_df['PostalCode'].astype(str)
postal_df['postal_code'] = postal_df['postal_code'].astype(str)

# Merge metadata with postal data to get coordinates for each postal code
metadata_df = metadata_df.merge(
    postal_df[['postal_code', 'lat', 'lon']],
    left_on='PostalCode', right_on='postal_code',
    how='inner'
)

# Remove rows where coordinates are missing after the merge
metadata_df.dropna(subset=['lat', 'lon'], inplace=True)
metadata_df['Coordinates'] = metadata_df[['lat', 'lon']].apply(tuple, axis=1)
print("Coordinates mapped from SG_postal.csv.")

# Create a new column merging Genotype and MajorLineage
metadata_df['Genotype_Lineage'] = metadata_df['Genotype'] + "_" + metadata_df['MajorLineage']

# --- Genotype-Lineage-Based Plotting and Distance Analysis ---

# Set seaborn style
sns.set(style="whitegrid", palette="colorblind")

# Define a color palette that is colorblind friendly with more contrast
palette = sns.color_palette("tab10")

# Plotting the coordinates colored by Genotype_Lineage for each serotype
serotypes = metadata_df['Serotype'].unique()
for serotype in serotypes:
    serotype_metadata = metadata_df[metadata_df['Serotype'] == serotype]
    genotype_lineages = serotype_metadata['Genotype_Lineage'].unique()
    
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
    
    for genotype_lineage in genotype_lineages:
        subset = serotype_metadata[serotype_metadata['Genotype_Lineage'] == genotype_lineage]
        
        # Calculate pairwise distances within this genotype-lineage for homotypic distance
        homotypic_coordinates = subset[['lat', 'lon']].values
        homotypic_dist_matrix = pairwise_distances(homotypic_coordinates, metric='euclidean')
        homotypic_mean = np.mean(homotypic_dist_matrix[np.triu_indices(len(homotypic_dist_matrix), k=1)])

        # Calculate pairwise distances between this genotype-lineage and others for heterotypic distance
        other_subset = serotype_metadata[serotype_metadata['Genotype_Lineage'] != genotype_lineage]
        if not other_subset.empty:
            other_coordinates = other_subset[['lat', 'lon']].values
            heterotypic_dist_matrix = pairwise_distances(homotypic_coordinates, other_coordinates, metric='euclidean')
            heterotypic_mean = np.mean(heterotypic_dist_matrix)
        else:
            heterotypic_mean = np.nan  # No other genotype-lineages available for comparison

        # Plot points for this genotype-lineage
        sns.scatterplot(
            x='lon', y='lat', data=subset, label=f"{genotype_lineage}",
            palette=palette, ax=ax, s=50, edgecolor='black', linewidth=0.5
        )
        
        # Display homotypic and heterotypic distance statistics for this genotype-lineage on the plot
        plt.text(
            0.05, 0.95 - 0.05 * np.where(genotype_lineages == genotype_lineage)[0][0],
            f'{genotype_lineage} - Homotypic: {homotypic_mean:.2f}, Heterotypic: {heterotypic_mean:.2f}',
            transform=ax.transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(facecolor='white', edgecolor='black', alpha=0.8)  # Opaque background with a border
        )

    plt.xlabel('Longitude', fontsize=16, weight='bold', labelpad=10)
    plt.ylabel('Latitude', fontsize=16, weight='bold', labelpad=10)
    plt.title(f'Genotype-Lineage Clustering for Serotype {serotype}', fontsize=18, weight='bold', pad=15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(title="Genotype_Lineage")
    plt.tight_layout()

    # Save the plot
    plot_path = os.path.join(plots_folder, f"genotype_lineage_clustering_map_serotype_{serotype}.png")
    plt.savefig(plot_path, format='png')
    plt.show()
    print(f"Plot saved as '{plot_path}'")
