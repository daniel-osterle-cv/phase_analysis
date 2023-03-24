# Confocal Microscopy Images Analysis

* [Project Description](#project-description)
* [Public Scripts Description](#public-scripts-description)
* [Private Scripts Description](#private-scripts-description)

## Project Description

The code in this repository is used to analyze raw data resulting from confocal microscopy images of *S. Cerevisiae* cells. In very general terms, these cells were prepared to contain two different building blocks that can assemble to form a mesoscale structure that can be observed as foci in microscopy imaging. The building blocks are based on proteins, namely a tetrameric protein with a binding domain A and a dimeric protein with a binding domain B. The two building blocks can assemble via the specific binding of domains A and B. Further, each building block is tagged with a specific fluorescent protein such that its location and abundance can be measured with microscopy. A more detailed description of these building blocks and their assembly can be found in [Nandi, Österle et al., 2022](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.129.128102) (PMID: 36179193).

This repository contains exemplary raw data from microscopy. This data is first processed using propietary scripts that have been developed in the [Emmanuel Levy Lab](https://www.weizmann.ac.il/CSB/ELevy/home). A short explanation of these scripts can be found in the section "Private Scripts Description". Running the analysis starting with raw data these private scripts are necessary. Please contact me at daniel.oesterle@mail.de to obtain them. 

## Public Scripts Description

### process_results.R
Execute this script from the command line with `./process_results.R` or run `./process_results.R --help` to see all available options that can be supplied from the command line. Without options the script will parse the directory from which this script is executed for a set of files that are required to postprocess the raw image files.
Most importantly, the script will look for the following files:
* A file of ".nd" format that contains basic information about the set-up of the (automated) microscopy screen
* A plate definition ("pdef") file of ".csv" format that contains information about the contents of each well on a 384-well plate
* Several files of ".txt" format that are the result of prior cell identification and segmentation and contain the raw imaging data, including cell size, position and fluorescence intensities for the measured wavelengths of a specific screen as described in [Heidenreich et al.,(2020)](https://rdcu.be/cE9xO) (PMID: 32661377).

Wile the files of ".nd" and ".txt" format result from the microscopy screen, the plate definition file has to be supplied by the user.

Further, this script then uses these input files and supplies them to the postprocess_result_txt() function to postprocess the raw data by mainly filtering unnecessary data with custom scripts from the [Emmanuel Levy Lab](https://www.weizmann.ac.il/CSB/ELevy/home) as described below and in this [preprint](https://doi.org/10.1101/260695).

Finally, the global environment is saved to "data_ready_to_analyze.RData". This file is sourced in "phase_analysis_main.R" for extracting the data of interest and plotting.


### phase_analysis_main.R

This script loads "phase_analysis_main.R" that contains all required functions and postprocessed data. Next, a new data frame is created by extracting columns of interest which include median measured fluorescence intensity, cellular dimensions, amount and dimensions of observed foci as well as identifying labels.
Finally, pdf files containing various plots are created by using functions from the "phase_plot_functions.R" module.

### phase_plot_functions.R

This module contains several functions to create plots with the ggplot2 R package. 

In the uploaded version of "phase_analysis_main.R" the function plot.overlap.gfp() is used which creates two x-y-scatterplots overlaid on top of each other. The bottom scatter plot uses all cells in which maximum GFP fluorescence intensity is below 100 AUF (arbitrary fluorescence units), whereas the top scatter plot uses cells where maximum GFP fluorescence intensity is above 1000 AUF. Biologically, this corresponds to low and high expression of a protein (here: a kinase) tagged with GFP, respectively. Both plots illustrate the median fluorescence intensity of protein scaffolds that are tagged with either BFP (y-axis) or RFP (x-axis).


## Private Scripts Description

### raw_data_tools.R

This module combines the functions for raw data processing based on the proprietary scripts. It calls the functions from the other private scripts and determines the order of data cleaning steps. Further, it extracts only the data that is actually required for further analysis and visualization and returns cleaned data in a data frame. It also performs conversion of fluorescence intensity values to protein concentrations.

### MicroscopeToolBox.R

Contains various functions to perform the processing of raw data including filtering background fluorescence, outliers and cells outside of a specified size range. A more detailed description can be found in this [preprint](https://doi.org/10.1101/260695).

### Multiwell.R

This module contains several function that are used to process multiwell plates such as 96 or 384 well plates and to ensure correct mapping of screening order to plate layout.

### general.R

This module contains several auxiliary functions related to color scales, matrix manipulation and memory checks.
