import pandas as pd
from geopy.distance import geodesic
import dendropy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Create 'plots' folder if it doesn't exist
plots_folder = "plots"
if not os.path.exists(plots_folder):
    os.makedirs(plots_folder)
    print(f"'{plots_folder}' folder created for saving plots.")

# File to save the results
results_file = "pairwise_evolutionary_spatial_distances.csv"

# Load metadata with columns: AccessionID, PostalCode, Serotype
metadata_file = "NEA_cleaned.metadata.csv"
metadata_df = pd.read_csv(metadata_file)

# Filter metadata to exclude entries with PostalCode equal to 0
metadata_df = metadata_df[metadata_df['PostalCode'] != 0]

# Load Singapore postal code data for coordinates
postal_data_file = "SG_postal.csv"
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

def calculate_tmrca(tree, sample1, sample2):
    try:
        # Find the MRCA of the two samples
        mrca = tree.mrca(taxon_labels=[sample1, sample2])
        if mrca is None:
            print(f"No MRCA found for {sample1} and {sample2}. Skipping this pair.")
            return None

        # Find the nodes for the two samples
        node1 = tree.find_node_with_taxon_label(sample1)
        node2 = tree.find_node_with_taxon_label(sample2)

        # Function to calculate distance from a node to an ancestor
        def distance_to_ancestor(node, ancestor):
            distance = 0.0
            current = node
            while current is not None and current != ancestor:
                if current.edge_length is not None:
                    distance += current.edge_length
                current = current.parent_node
            return distance

        # Calculate the evolutionary distance from each sample to the MRCA
        distance_to_mrca_sample1 = distance_to_ancestor(node1, mrca)
        distance_to_mrca_sample2 = distance_to_ancestor(node2, mrca)

        # Sum the distances from each sample to the MRCA
        evolutionary_time = distance_to_mrca_sample1 + distance_to_mrca_sample2

        return evolutionary_time
    except KeyError:
        print(f"One or both of the samples ({sample1}, {sample2}) are not found in the tree. Skipping this pair.")
        return None

# Check if results file exists
if os.path.exists(results_file):
    print(f"Loading results from '{results_file}'")
    results_df = pd.read_csv(results_file)
else:
    # Dictionary to store tree files by serotype
    tree_files = {
        "D1": "data/NEA_D1.tree",
        "D2": "data/NEA_D2.tree",
        "D3": "data/NEA_D3.tree",
        "D4": "data/NEA_D4.tree"
    }

    # Store results across all serotypes
    all_results = []

    # Process each serotype independently, adding to the combined dataset
    for serotype, tree_file in tree_files.items():
        print(f"Processing serotype {serotype} with tree {tree_file}...")

        # Load the phylogenetic tree for the specific serotype using DendroPy
        tree = dendropy.Tree.get(
            path=tree_file,
            schema="nexus",
            preserve_underscores=True
        )

        # Filter metadata for the current serotype
        serotype_metadata = metadata_df[metadata_df['Serotype'] == serotype]
        num_pairs = (len(serotype_metadata) * (len(serotype_metadata) - 1)) // 2
        print(f"Calculating pairwise distances and MRCAs for {len(serotype_metadata)} samples ({num_pairs} pairs).")

        # Calculate pairwise distances and tMRCA for each pair within the serotype
        pair_count = 0
        for i, row1 in serotype_metadata.iterrows():
            for j, row2 in serotype_metadata.iterrows():
                if i < j:
                    sample1, sample2 = row1['AccessionID'], row2['AccessionID']
                    coords1, coords2 = row1['Coordinates'], row2['Coordinates']

                    # Calculate spatial distance
                    spatial_distance = geodesic(coords1, coords2).kilometers

                    # Calculate tMRCA
                    tmrca = calculate_tmrca(tree, sample1, sample2)
                    if tmrca is None:
                        continue

                    # Append results to the main list
                    all_results.append({
                        'Sample1': sample1,
                        'Sample2': sample2,
                        'EvolutionaryTime': tmrca,
                        'SpatialDistance': spatial_distance
                    })

                    pair_count += 1
                    if pair_count % 100 == 0:
                        print(f"Processed {pair_count} / {num_pairs} pairs for serotype {serotype}")

        print(f"Completed serotype {serotype}.")

    # Convert results to DataFrame
    if all_results:
        results_df = pd.DataFrame(all_results)
        print("Pairwise calculations complete.")
        results_df.to_csv(results_file, index=False)
        print(f"Results saved as '{results_file}'")
    else:
        print("No results were collected. Please check for errors in the tree or metadata files.")

# --- Plotting Section ---

if 'EvolutionaryTime' not in results_df.columns or 'SpatialDistance' not in results_df.columns:
    print("No evolutionary time or spatial distance data available.")
else:
    # Set seaborn style
    sns.set(style="whitegrid", palette="colorblind")

    # --- Line Plot with Confidence Interval using Seaborn ---

    # Define finer bins for evolutionary time intervals (e.g., 0.1-year intervals)
    bin_width = 0.1
    bins = np.arange(0, results_df['EvolutionaryTime'].max() + bin_width, bin_width)
    results_df['TimeBin'] = pd.cut(results_df['EvolutionaryTime'], bins)

    # Filter out bins with no entries
    non_empty_bins = results_df.dropna(subset=['EvolutionaryTime', 'SpatialDistance'])
    grouped = non_empty_bins.groupby('TimeBin')

    # Calculate median and confidence intervals for each non-empty bin
    agg_data = grouped.agg(
        MedianDistance=('SpatialDistance', 'median'),
        CI_Lower=('SpatialDistance', lambda x: np.percentile(x, 2.5) if len(x) > 0 else np.nan),
        CI_Upper=('SpatialDistance', lambda x: np.percentile(x, 97.5) if len(x) > 0 else np.nan)
    ).reset_index()

    # Calculate the midpoint of each bin
    agg_data['TimeBinCenter'] = agg_data['TimeBin'].apply(lambda x: x.mid if pd.notnull(x) else np.nan)
    agg_data.dropna(subset=['TimeBinCenter'], inplace=True)

    # Set up the plot with seaborn
    plt.figure(figsize=(10, 8), dpi=300)
    sns.lineplot(
        x='TimeBinCenter', y='MedianDistance',
        data=agg_data, marker='o', label='Median Spatial Distance'
    )

    plt.fill_between(agg_data['TimeBinCenter'], agg_data['CI_Lower'], agg_data['CI_Upper'], color='lightblue', alpha=0.3, label='95% Confidence Interval')
    plt.xlim(0, 10)
    plt.xlabel('Evolutionary Time (years)', fontsize=16, weight='bold', labelpad=10)
    plt.ylabel('Spatial Distance (km)', fontsize=16, weight='bold', labelpad=10)
    plt.title('Spatial Distance Over Evolutionary Time (Finer Bins)', fontsize=18, weight='bold', pad=15)
    plt.legend(fontsize=14, loc='upper right', frameon=False)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()

    plot_path_line = os.path.join(plots_folder, "finer_evolutionary_time_vs_spatial_distance_seaborn.png")
    plt.savefig(plot_path_line, format='png')
    plt.show()
    print(f"Plot saved as '{plot_path_line}'")
