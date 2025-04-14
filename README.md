# Angle-Domain Multi-User Multiplexing for Non-Coherent MIMO in Rician channels

## Description

This repository contains the code used to obtain the results of the article "Angle-Domain Multi-User Multiplexing for Non-Coherent MIMO in Rician channels". The folder "Functions" contains the functions needed to
run the simulation, the folder "Development_scripts" contains the files used for the coding and development of the simulation files. 

## Getting Started

### Dependencies

* The code was done using Matlab 2024b.

### Running the simulations

* Running main.m will create a folder named after the current date with the results of the simulation with the same parameters used in the article.
* Running main_antenna_scaling.m will create a folder called "antenna_scaling" where the simulation results for increasing numbers of antennas are stored.
* Include the name of the result folder created by main.m in the folder variable in line 3 of Plot_results.m
* Running the first half of Plot_results.m shows the simulation results corresponding to Figs. 2 and 3 from the paper.
* Running the second half of Plot_results.m shows the simulation results corresponding to Figs. 4.

## Authors

Ruben de Miguel Gil (rmgil@tsc.uc3m.es)
