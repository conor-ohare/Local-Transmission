# Local Structure Analysis of Dengue Transmission

## Summary

This project is inspired by *"Dengue diversity across spatial and temporal scales: local structure and the impact of host population size"* (Salje et al., Science, 2017). It aims to analyze local dengue transmission dynamics using geolocated case data, serotype classification, and phylogenetic methods. This project explores the structure of transmission chains, their spatial and temporal distribution, and the effect of population density on dengue spread.

### Key Features

- **Local Transmission Analysis**: Determines the proportion of dengue cases in close proximity that belong to the same transmission chain.
- **Population Density and Chain Saturation**: Evaluates the impact of population density on transmission chain diversity and saturation.
- **Temporal and Spatial Dispersion**: Examines how transmission chains change over time and distance.

---

## Scripts Description

### 1. `prop_homotypic.py`
   - **Purpose**: Analyzes the proportion of case-pairs from the same serotype (homotypic) as a function of distance.
   - **Description**: This script calculates the proportion of case-pairs that are homotypic within varying distance ranges, highlighting the likelihood of close-proximity transmission for the same serotype.

### 2. `prop_same_chain.py`
   - **Purpose**: Calculates the proportion of case-pairs that belong to the same transmission chain over different distances.
   - **Description**: This script evaluates pairs of cases to determine if they are part of the same transmission chain, showing the relationship between distance and localized dengue transmission.

### 3. `time_vs_distance.py`
   - **Purpose**: Explores the relationship between spatial distance and evolutionary time.
   - **Description**: This script analyzes spatial distance as a function of evolutionary time, providing insights into how far the virus spreads over different evolutionary time scales.

---

## Requirements

- **Python 3.x**
- **Packages**:
  - `pandas`: Data manipulation and analysis
  - `numpy`: Numerical operations
  - `geopy`: Calculating distances between geographic coordinates
  - `dendropy`: Phylogenetic tree analysis
  - `seaborn`: Data visualization
  - `matplotlib`: Plotting and visualization
  - `scipy`: Statistical calculations
  - `datetime`: Date and time manipulation
  - `os`: File and directory handling