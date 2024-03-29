[![DOI](https://zenodo.org/badge/668270612.svg)](https://zenodo.org/doi/10.5281/zenodo.10666864)

# SwiFColBM_Evo NetLogo implementation
## Description

The SwiFCoIBM_evo model is an individual-based simulation model that simulates the dynamics of Classical Swine Fever (CSF) virus in a wild boar population. It incorporates explicit movement, trait evolution, and dynamic landscape processes. The model is based on previous versions published in Kramer-Schadt et al. 2009, Lange et al. 2012, Scherer et al. 2020, and Kürschner et al. 2021.

The model consists of two major components:

Wild Boar Demography Model: This component simulates the demographic processes of wild boar populations, including seasonal reproduction, herd splitting (dispersal), and mortality. It models the population dynamics of wild boars based on various parameters such as survival probabilities, reproduction probabilities, and movement distances.

CSF Virus Model: This component simulates the spread of the CSF virus within the wild boar population. It considers factors such as transmission coefficients, infection probabilities, and infectious periods. The model tracks the number of infected individuals, the strain of the virus, and the spatial distribution of infections.

# R Markdown analysis script
## Description
The R Markdown script is used for creating the figures and analysis used in the publication "Resource asynchrony and landscape homogenization as drivers of virulence evolution: the case of a directly transmitted disease in a social host" by Tobias Kürschner, Cédric Scherer, Viktoriia Radchuk, Niels Blaum and Stephanie Kramer-Schadt.
It performs various data processing and analysis tasks on the raw data generated by the IBM model SwiFColBM_Evo.

## Prerequisites
- Required R packages installed (see `source_lib.r` for package dependencies)

## Usage
1. Set up the R environment by running the setup chunk at the beginning of the script.
2. Modify the file paths in the script to point to the correct input and output directories.
3. Run each chunk sequentially to perform the desired data processing and analysis tasks.
4. The final output will be saved in the specified output folder.

## Author
- Tobias Kuerschner

## Date
- 04/05/2021
