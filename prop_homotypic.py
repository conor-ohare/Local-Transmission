import pandas as pd
import numpy as np
from geopy.distance import geodesic
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import timedelta
import os

# Function for pretty printing messages with an icon
def pretty_print(message, icon="✨"):
    print(f"{icon} {message}")

# Create 'plots' folder if it doesn't exist to store plot files
plots_folder = "plots"
if not os.path.exists(plots_folder):
    os.makedirs(plots_folder)
    pretty_print(f"'{plots_folder}' folder created for saving plots.", icon="📂")

# Path for saving analysis results
results_file = "homotypic_proportion_results.csv"

# Function to calculate and save homotypic proportion analysis results
def calculate_homotypic_proportion(metadata_df, distance_bins):
    pretty_print("Starting homotypic proportion calculation...", icon="🔍")
    
    # Initialize bin centers and counters for homotypic and total pairs per bin
    bin_centers = (distance_bins[:-1] + distance_bins[1:]) / 2
    homotypic_counts = [0] * (len(distance_bins) - 1)
    total_counts = [0] * (len(distance_bins) - 1)

    # Iterate over metadata rows to perform pairwise comparisons
    for i, row1 in metadata_df.iterrows():
        for j, row2 in metadata_df.iterrows():
            if i < j:
                # Calculate spatial distance between two cases
                spatial_distance = geodesic(row1['Coordinates'], row2['Coordinates']).kilometers

                # Identify bin for the spatial distance
                bin_index = np.digitize(spatial_distance, distance_bins) - 1
                if 0 <= bin_index < len(total_counts):
                    # Calculate time difference (in months) between cases
                    time_diff_days = abs((row1['CollectionDate'] - row2['CollectionDate']).days)
                    time_diff_months = time_diff_days / 30.44
                    if time_diff_months > 6:
                        continue  # Skip pairs that are not within 6 months

                    # Increment total pair count in the distance bin
                    total_counts[bin_index] += 1

                    # Check if the cases are of the same serotype (homotypic)
                    if row1['Serotype'] == row2['Serotype']:
                        homotypic_counts[bin_index] += 1

    # Calculate homotypic proportions and confidence intervals for each distance bin
    proportions = []
    confidence_intervals = []
    for homotypic_count, total_count in zip(homotypic_counts, total_counts):
        if total_count > 0:
            proportion = homotypic_count / total_count
            proportions.append(proportion)
            
            # Calculate 95% confidence interval for the proportion
            ci = 1.96 * np.sqrt((proportion * (1 - proportion)) / total_count)
            confidence_intervals.append(ci)
        else:
            proportions.append(np.nan)
            confidence_intervals.append(np.nan)
    
    # Save results to a DataFrame and write to CSV
    results_df = pd.DataFrame({
        "DistanceBinCenter": bin_centers,
        "ProportionHomotypic": proportions,
        "ConfidenceInterval": confidence_intervals
    })
    results_df.to_csv(results_file, index=False)
    pretty_print(f"Results saved to '{results_file}'", icon="💾")
    return results_df

# Check if results file exists, otherwise, calculate and save the results
if os.path.exists(results_file):
    pretty_print(f"Loading results from '{results_file}'", icon="📄")
    results_df = pd.read_csv(results_file)
else:
    pretty_print("Calculating homotypic proportions...", icon="⚙️")

    # Load metadata and filter out entries with PostalCode = 0
    metadata_file = "NEA_cleaned.metadata.csv"  # Update with your metadata file path
    metadata_df = pd.read_csv(metadata_file)
    metadata_df['CollectionDate'] = pd.to_datetime(metadata_df['CollectionDate'])
    metadata_df = metadata_df[metadata_df['PostalCode'] != 0]
    pretty_print(f"Metadata loaded with {len(metadata_df)} entries (PostalCode 0 excluded).", icon="📊")

    # Load postal code data and merge to add geographical coordinates
    postal_data_file = "SG_postal.csv"  # Update with the path to your postal data file
    postal_df = pd.read_csv(postal_data_file)
    pretty_print("Postal data loaded.", icon="🌍")

    # Merge postal data to add lat/lon coordinates to metadata
    metadata_df['PostalCode'] = metadata_df['PostalCode'].astype(str)
    postal_df['postal_code'] = postal_df['postal_code'].astype(str)
    metadata_df = metadata_df.merge(postal_df[['postal_code', 'lat', 'lon']], left_on='PostalCode', right_on='postal_code', how='inner')
    metadata_df.dropna(subset=['lat', 'lon'], inplace=True)
    metadata_df['Coordinates'] = metadata_df[['lat', 'lon']].apply(tuple, axis=1)
    pretty_print("Coordinates mapped from postal data.", icon="📍")

    # Define distance bins (e.g., every 0.2 km up to 5.4 km)
    distance_bins = np.arange(0, 5.4, 0.2)
    results_df = calculate_homotypic_proportion(metadata_df, distance_bins)

# Plotting results with confidence intervals and enhanced aesthetics
sns.set(style="whitegrid")
plt.figure(figsize=(10, 8), dpi=300)

# Plot homotypic proportions with confidence intervals shaded
plt.plot(results_df["DistanceBinCenter"], results_df["ProportionHomotypic"], marker="s", color="#1f77b4", markersize=8, linewidth=2.5, label="Proportion homotypic")
plt.fill_between(
    results_df["DistanceBinCenter"], 
    results_df["ProportionHomotypic"] - results_df["ConfidenceInterval"], 
    results_df["ProportionHomotypic"] + results_df["ConfidenceInterval"], 
    color="#1f77b4", 
    alpha=0.2,
    label="95% CI"
)

# Reference line at y = 0.5
plt.axhline(y=0.5, color='gray', linestyle='--', linewidth=1.5, label="y = 0.5")

# Axis labels and title for publication quality
plt.xlabel("Distance (km)", fontsize=16, weight="bold", labelpad=10)
plt.ylabel("Proportion of Homotypic Pairs", fontsize=16, weight="bold", labelpad=10)
plt.ylim(0, 1)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Proportion of Homotypic Case-Pairs vs Distance", fontsize=18, weight="bold", pad=15)

# Add grid, legend, and finalize layout
plt.grid(visible=True, linestyle='--', linewidth=0.5, alpha=0.7)
sns.despine(trim=True)
plt.legend(fontsize=14, loc='upper right', frameon=False)
plt.tight_layout()

# Save the plot in PDF and PNG formats
plot_path_pdf = os.path.join(plots_folder, "homotypic_proportion_vs_distance.pdf")
plot_path_png = os.path.join(plots_folder, "homotypic_proportion_vs_distance.png")
plt.savefig(plot_path_pdf, format='pdf')
plt.savefig(plot_path_png, format='png')
plt.show()

pretty_print(f"Plot saved as '{plot_path_pdf}' and '{plot_path_png}'", icon="📊")
