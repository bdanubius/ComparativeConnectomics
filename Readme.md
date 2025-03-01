# Physical Network Constraints Define the Lognormal Architecture of the Brain's Connectome

Thank you for your interest in our project. 

This repository contains the code and data access to reproduce the results of [our paper](https://www.biorxiv.org/content/10.1101/2025.02.27.640551v1), organized into several directories for ease of use. Below is a description of the folder structure and contents.

## Directory Structure

### Datasets/

This directory contains subfolders for each system (e.g., Larva, FlyWire). Within each system’s subdirectory, there are two folders:

#### Skeletons/ 
– Raw skeleton data (.swc files) for each dataset.

#### Synapses/ 
– Raw synapse data for each dataset.

### Processing/

This directory contains six subfolders, corresponding to different stages of data processing:

1. SWC Processing – Converts .swc files into a set of connected chains or segments by splitting the skeleton tree at nodes with degree ≥ 3.

2. Synapse Processing – Matches each synapse point to the nearest position in the corresponding skeleton.

3. Geometry Processing – Computes geometric properties (length, curvature, torsion) from processed .swc data and saves aggregated data in .mat format.

4. Degree Processing – Uses synapse data to count the degree and total synapses for each neuron.

5. Data Aggregation – Compiles processed data into spreadsheets, which are directly used to generate figures.

6. Filtering – Applies dataset-specific filtering for Hemibrain, MANC, and FlyWire as specified in the supplementary information.

### Processed_Data/

This directory contains data at various stages of processing. Running the scripts in Processing/ will populate this directory with intermediate and final processed data.

### Figure_Generation/

This directory contains scripts for generating figures using the processed data:

Figure 2 Generation – Uses aggregated data spreadsheets to reproduce Figure 2.

Helper Functions – Contains additional functions used for plotting.

More figure scripts will be added in the future.

### Figures/

This directory stores the .PDF figures generated from the scripts in Figure_Generation/.

## Usage Notes

If you have downloaded the GitHub version of this project, only the Larva dataset preprocessing scripts and figure generation will be fully functional.

If you have downloaded the full dataset, all scripts should work as expected, including data and processed results.

Regardless of the version, all datasets include the final aggregated spreadsheets (From Step 5: Data Aggregation), enabling figure reproduction.

## Future Updates

We plan to expand the Figure_Generation/ directory with additional scripts to facilitate figure reproduction and analysis.

For any issues or questions, please feel free to open an issue or reach out.
