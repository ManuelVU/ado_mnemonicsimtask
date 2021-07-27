Adaptive Design Optimization for a Mnemonic Similatiry Task
===========================================================

This repository contains the results from the "Adaptive Design Optimization for a Mnemonic Similatiry Task" working paper. Files are divided in 4 main folders:

* __data__ : Contains experimental data and results form analysis.

* __figures__: Contains Figures used on the paper and additional graphs at the individual level.

* __models__: Contains files in '.txt' format of Bayesian graphical models used by other functions.

* __src__: Contains the main functions used in the analysis and adaptive design procedure.

# Data files

There are 5 data files in the __data__ directory:

* __StarkExp1.mat__: Contains the raw data in the experiment. This file is parsed by the functions in __data_t.R__ in the directory __src__ to return the data file used for the analysis in R.

* __memory.RData__: This file contains the data from the experiment that will be used by other functions. The R object is organized as a list with two elements:
    1. *results$rs*: is an array with 3 dimentions, trials (192), participants by age group (20) and age group (1 = Young, 2 = Elderly) 