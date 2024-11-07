import os  # Import os module for checking file existence
import pandas as pd
import numpy as np
from geopy.distance import geodesic
import dendropy
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import timedelta
from scipy.stats import linregress, sem


# Function to calculate branch length from MRCA to recent sample manually
def calculate_branch_length(mrca, recent_node):
    branch_length = 0.0
    current = recent_node
    while current != mrca:
        if current.edge_length is not None:
            branch_length += current.edge_length
        current = current.parent_node
    return branch_length

# Function to determine if MRCA is within six months of the earlier sample
def is_mrca_within_six_months(tree, early_sample, recent_sample, early_date, recent_date):
    try:
        # Find the MRCA of the two samples
        mrca = tree.mrca(taxon_labels=[early_sample, recent_sample])
    except KeyError:
        # Skip if either sample is not in the tree
        print(f"Skipping comparison: {early_sample} or {recent_sample} not found in tree.")
        return False, None

    if mrca is None:
        return False, None

    # Find the recent sample node and calculate the branch length to MRCA
    try:
        recent_node = tree.find_node_with_taxon_label(recent_sample)
    except KeyError:
        print(f"Skipping comparison: {recent_sample} not found in tree.")
        return False, None

    branch_length_to_recent = calculate_branch_length(mrca, recent_node)

    # Calculate the inferred MRCA date relative to the recent sample
    mrca_date = recent_date - timedelta(days=branch_length_to_recent * 365.25)  # Convert branch length (years) to days

    # Check if the MRCA date is within six months of the early sample's collection date
    is_within_six_months = abs((early_date - mrca_date).days) <= 183  # 183 days ~ 6 months
    return is_within_six_months, branch_length_to_recent

# Function for pretty printing
def pretty_print(message, icon="✨"):
    print(f"{icon} {message}")

# Load metadata and filter out PostalCode = 0
metadata_file = "NEA_cleaned.metadata.csv"  # Update with your metadata file path
pretty_print(f"Loading metadata from '{metadata_file}'...", "📂")
metadata_df = pd.read_csv(metadata_file)
metadata_df['CollectionDate'] = pd.to_datetime(metadata_df['CollectionDate'])
metadata_df = metadata_df[metadata_df['PostalCode'] != 0]
pretty_print(f"Metadata loaded with {len(metadata_df)} entries (PostalCode 0 excluded).", "✅")

# Load postal data
postal_data_file = "SG_postal.csv"  # Update with the path to your postal data file
pretty_print(f"Loading postal data from '{postal_data_file}'...", "📍")
postal_df = pd.read_csv(postal_data_file)
pretty_print("Postal data loaded successfully.", "✅")

# Merge postal data to add coordinates
pretty_print("Mapping coordinates from postal data...", "🗺️")
metadata_df['PostalCode'] = metadata_df['PostalCode'].astype(str)
postal_df['postal_code'] = postal_df['postal_code'].astype(str)
metadata_df = metadata_df.merge(postal_df[['postal_code', 'lat', 'lon']], left_on='PostalCode', right_on='postal_code', how='inner')
metadata_df.dropna(subset=['lat', 'lon'], inplace=True)
metadata_df['Coordinates'] = metadata_df[['lat', 'lon']].apply(tuple, axis=1)
pretty_print("Coordinates mapped successfully.", "📌")

# Define tree files by serotype
tree_files = {
    "D1": "data/NEA_D1.tree",
    "D2": "data/NEA_D2.tree",
    "D3": "data/NEA_D3.tree",
    "D4": "data/NEA_D4.tree"
}

# Prepare smaller distance bins (e.g., every 0.1 km up to 5.4 km)
distance_bins = np.arange(0, 5.4, 0.1)
bin_centers = (distance_bins[:-1] + distance_bins[1:]) / 2
same_chain_counts = [0] * (len(distance_bins) - 1)
total_counts = [0] * (len(distance_bins) - 1)

# Define path for pairwise comparison results CSV
pairwise_results_csv = "pairwise_comparisons_results.csv"

# Check if the CSV already exists
if not os.path.exists(pairwise_results_csv):
    pretty_print("No existing pairwise comparison results found.", "🔍")
    pretty_print("Running pairwise comparison analysis...", "⚙️")

    # List to store pairwise comparison results
    pairwise_results = []

    # Process each serotype independently
    for serotype, tree_file in tree_files.items():
        pretty_print(f"Processing serotype {serotype}...", "🧬")

        # Load the phylogenetic tree for the specific serotype
        tree = dendropy.Tree.get(path=tree_file, schema="nexus", preserve_underscores=True)
        pretty_print(f"Tree for {serotype} loaded successfully.", "🌳")

        # Filter metadata for the current serotype
        serotype_metadata = metadata_df[metadata_df['Serotype'] == serotype]
        if serotype_metadata.empty:
            pretty_print(f"No samples found for serotype {serotype}. Skipping...", "⚠️")
            continue  # Skip this serotype if there are no records

        pretty_print(f"{len(serotype_metadata)} samples loaded for serotype {serotype}.", "📊")

        # Initialize a counter for pairwise comparisons
        pairwise_comparison_count = 0

        # Calculate pairwise distances and MRCA for each pair within the serotype
        for i, row1 in serotype_metadata.iterrows():
            for j, row2 in serotype_metadata.iterrows():
                if i < j:
                    pairwise_comparison_count += 1

                    # Calculate spatial distance
                    spatial_distance = geodesic(row1['Coordinates'], row2['Coordinates']).kilometers

                    # Bin the pair by distance
                    bin_index = np.digitize(spatial_distance, distance_bins) - 1
                    if 0 <= bin_index < len(total_counts):
                        total_counts[bin_index] += 1  # Count all pairs in the bin

                        # Check if the cases are within 6 months of each other
                        time_diff_days = abs((row1['CollectionDate'] - row2['CollectionDate']).days)
                        time_diff_months = time_diff_days / 30.44
                        if time_diff_months > 6:
                            continue  # Skip if they aren't within 6 months of each other

                        # Determine the earlier and recent sample based on CollectionDate
                        if row1['CollectionDate'] <= row2['CollectionDate']:
                            early_sample, recent_sample = row1['AccessionID'], row2['AccessionID']
                            early_date, recent_date = row1['CollectionDate'], row2['CollectionDate']
                        else:
                            early_sample, recent_sample = row2['AccessionID'], row1['AccessionID']
                            early_date, recent_date = row2['CollectionDate'], row1['CollectionDate']

                        # Check if MRCA date is within six months of the earlier sample's date
                        mrca_within_six_months, branch_length_to_recent = is_mrca_within_six_months(
                            tree, early_sample, recent_sample, early_date, recent_date
                        )
                        if mrca_within_six_months:
                            same_chain_counts[bin_index] += 1  # Only count pairs that satisfy the same chain criteria

                        # Save results to pairwise_results list
                        pairwise_results.append({
                            'Serotype': serotype,
                            'Sample1': early_sample,
                            'Sample2': recent_sample,
                            'Distance': spatial_distance,
                            'SameChain': mrca_within_six_months,
                            'BinIndex': bin_index
                        })

                    # Debug statement to show progress of adding data
                    if pairwise_comparison_count % 500 == 0:
                        pretty_print(f"{pairwise_comparison_count} pairwise comparisons done for serotype {serotype}.", "🔄")

        pretty_print(f"Completed serotype {serotype}. Processed {pairwise_comparison_count} total pairs.", "✅")

    # Save pairwise results to a CSV file
    if len(pairwise_results) > 0:
        pairwise_results_df = pd.DataFrame(pairwise_results)
        pairwise_results_df.to_csv(pairwise_results_csv, index=False)
        pretty_print(f"Pairwise comparison results saved as '{pairwise_results_csv}'", "💾")
    else:
        pretty_print("No pairwise comparisons met the criteria for any serotype. CSV will not be saved.", "⚠️")

else:
    pretty_print(f"Loading pairwise comparison results from '{pairwise_results_csv}'...", "📄")
    pairwise_results_df = pd.read_csv(pairwise_results_csv)
    pretty_print(f"Pairwise comparison results loaded successfully from '{pairwise_results_csv}'", "✅")

    # Recalculate same_chain_counts and total_counts from the loaded data
    for _, row in pairwise_results_df.iterrows():
        bin_index = row['BinIndex']
        if 0 <= bin_index < len(total_counts):
            total_counts[bin_index] += 1
            if row['SameChain']:
                same_chain_counts[bin_index] += 1

# Calculate proportions and confidence intervals for each distance bin
pretty_print("Calculating proportions and confidence intervals for each distance bin...", "📊")
proportions = []
confidence_intervals = []
effective_chains = []

for same_count, total_count in zip(same_chain_counts, total_counts):
    if total_count > 0:
        proportion = same_count / total_count
        effective_chain = 1 / proportion
        ci = 1.96 * sem([effective_chain] * total_count)  # Using SEM and scaling for 95% CI
        proportions.append(proportion)
        effective_chains.append(effective_chain)
        confidence_intervals.append(ci)
    else:
        proportions.append(np.nan)
        effective_chains.append(np.nan)
        confidence_intervals.append(np.nan)

pretty_print("Proportions and confidence intervals calculated successfully.", "✅")

# Create 'plots' folder if it doesn't exist
plots_folder = "plots"
if not os.path.exists(plots_folder):
    os.makedirs(plots_folder)
    pretty_print(f"'{plots_folder}' folder created for saving plots.", "📁")

# Plotting with confidence intervals and enhanced aesthetics
pretty_print("Generating plot: Proportion of Same-Chain Pairs vs Distance...", "📊")
sns.set(style="whitegrid")
plt.figure(figsize=(10, 8), dpi=300)  # Larger figure size and higher DPI for better quality

# Plot with square markers, thicker lines, and confidence intervals with shading
plt.plot(bin_centers, proportions, marker="s", color="#1f77b4", markersize=8, linewidth=2.5, label="Proportion same chain")
plt.fill_between(
    bin_centers, 
    np.array(proportions) - np.array(confidence_intervals), 
    np.array(proportions) + np.array(confidence_intervals), 
    color="#1f77b4", 
    alpha=0.2,
    label="95% CI"
)

# Customize the plot with professional aesthetics
plt.xlabel("Distance (km)", fontsize=16, weight="bold", labelpad=10)
plt.ylabel("Proportion of Same-Chain Pairs", fontsize=16, weight="bold", labelpad=10)
plt.ylim(0, 0.8)  # Set y-axis limit similar to example image
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Proportion of Same-Chain Pairs vs Distance", fontsize=18, weight="bold", pad=15)

# Add subtle grid lines and remove unnecessary spines
plt.grid(visible=True, linestyle='--', linewidth=0.5, alpha=0.7)
sns.despine(trim=True)

# Add a legend
plt.legend(fontsize=14, loc='upper right', frameon=False)

# Final adjustments
plt.tight_layout()

# Save the plot as a vector graphic (PDF) and PNG
plot_path_pdf = os.path.join(plots_folder, "same_chain_proportion_with_confidence_intervals_smaller_bins.pdf")
plot_path_png = os.path.join(plots_folder, "same_chain_proportion_with_confidence_intervals_smaller_bins.png")
plt.savefig(plot_path_pdf, format='pdf')
plt.savefig(plot_path_png, format='png')
plt.show()

pretty_print(f"Plot saved as '{plot_path_pdf}' and '{plot_path_png}'", "📊")

# Calculate the effective number of transmission chains
pretty_print("Calculating effective number of transmission chains...", "🔢")
effective_chains = []
for proportion in proportions:
    if proportion > 0:
        effective_chain = 1 / proportion
        effective_chains.append(effective_chain)
    else:
        effective_chains.append(np.nan)

pretty_print("Effective number of transmission chains calculated successfully.", "✅")

# Average population density of Singapore (people per square kilometer)
population_density = 8383.25  # people/km²

# Debugging print to verify confidence intervals
for i, (chain, ci) in enumerate(zip(effective_chains, confidence_intervals)):
    print(f"Bin {i}: Effective Chains = {chain}, 95% CI = ±{ci}")

# Calculate population size for each distance bin (assuming the distance is the radius of a circle)
areas = np.pi * (bin_centers ** 2)  # Area of circle with radius equal to distance
population_sizes = areas * population_density  # Calculate population size for each bin

# Log-transform population sizes for linear fitting
log_population_sizes = np.log10(population_sizes)
effective_chains = np.array(effective_chains)
confidence_intervals = np.array(confidence_intervals)

# Remove NaN values for regression calculation
valid_indices = ~np.isnan(effective_chains) & ~np.isinf(log_population_sizes)
log_population_sizes = log_population_sizes[valid_indices]
effective_chains = effective_chains[valid_indices]
confidence_intervals = confidence_intervals[valid_indices]

# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(log_population_sizes, effective_chains)

# Generate predicted values for the regression line
regression_line = slope * log_population_sizes + intercept

# Plotting the effective number of transmission chains vs population size
plt.figure(figsize=(10, 8), dpi=300)  # Larger figure size and higher DPI for better quality

# Plot effective number of chains vs population size with confidence intervals
plt.errorbar(population_sizes, effective_chains, yerr=confidence_intervals, fmt='o', color="#ff7f0e", markersize=6, 
             ecolor='orange', capsize=5, elinewidth=2, markeredgewidth=1, label="Effective number of chains")

# Plot the regression line
plt.plot(10**log_population_sizes, regression_line, color='blue', linestyle='--', linewidth=2, label="Line of best fit")

# Set log scale for x-axis (population size)
plt.xscale('log')

# Set labels and title
plt.xlabel("Population Size (log scale)", fontsize=16, weight="bold", labelpad=10)
plt.ylabel("Effective Number of Transmission Chains", fontsize=16, weight="bold", labelpad=10)
plt.ylim(0, max(effective_chains + confidence_intervals) + 1)  # Set y-axis limit to fit data
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Effective Number of Transmission Chains vs Population Size", fontsize=18, weight="bold", pad=15)

# Add subtle grid lines and remove unnecessary spines
plt.grid(visible=True, linestyle='--', linewidth=0.5, alpha=0.7)
sns.despine(trim=True)

# Add a secondary x-axis at the top for distance in km
ax = plt.gca()
secax = ax.secondary_xaxis('top')
secax.set_xscale('log')
secax.set_xlim(ax.get_xlim())
secax.set_xlabel('Distance (km)', fontsize=14, weight='bold', labelpad=10)

# Define the ticks for the top axis to match relevant distance values in kilometers
distance_ticks = [0.2, 0.5, 1, 2, 5]  # Selected relevant distance ticks for better readability
secax.set_xticks([np.pi * (d ** 2) * population_density for d in distance_ticks])
secax.set_xticklabels([str(d) for d in distance_ticks], fontsize=12)

# Add a legend
plt.legend(fontsize=14, loc='upper left', frameon=False)

# Final adjustments
plt.tight_layout()

# Create 'plots' folder if it doesn't exist
plots_folder = "plots"
if not os.path.exists(plots_folder):
    os.makedirs(plots_folder)

# Save the plot as a vector graphic (PDF) and PNG
plot_path_pdf = os.path.join(plots_folder, "effective_number_of_chains_vs_population_size.pdf")
plot_path_png = os.path.join(plots_folder, "effective_number_of_chains_vs_population_size.png")
plt.savefig(plot_path_pdf, format='pdf')
plt.savefig(plot_path_png, format='png')
plt.show()

pretty_print(f"Plot saved as '{plot_path_pdf}' and '{plot_path_png}'", "📈")