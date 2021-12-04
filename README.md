# 2022-BMB-FS-JRM

# Adaptive Evolution of a Species' Range

Written by Farshad Shirani (f.shirani@gatech.edu), December 3, 2021

This repository contains MATLAB codes used in the article:

"Competition, Trait Variance Dynamics, and the Evolution of a Species' Range", F. Shirani and J. R. Miller, *Bulletin of Mathematical Biology*, (2022). 

All codes are written in MATLAB R2021a.

## Description
These MATLAB codes numerically compute solutions of several systems of partaial differntial equations that simulae adaptive evoultion of the geographic range of biological species in 1- and 2-dimensional spaces.

# How to run simulations:

-	Codes for a single species as well as two species, both in 1-dimensional and 2-dimensional geographic spaces, are provided in separate folders.

-	Each folder contains a function named `initializeSimulation`. Model parameters, simulation parameters, and discretization parameters of the numerical scheme are set by the user through this function. The model parameters are defined as global parameters. Therefore, the values of model parameters should only be changed through this function.

-	Each folder contains a script named `Main_VariableVariance`, which is the main code that performs the simulations based on the parameters set in the function `initializeSimulation`. When a simulation is complete, the resulting solutions can be saved using the `save` command provided at the end of the script. The computed solutions are stored in a structure array named `populations`. The parameters are also stored in three different structures, named `modelParameters`, `simulationParameters`, and `discretizationParamaters`. 

-	The results can be plotted using the script `PlottingResults` available in each folder. 
