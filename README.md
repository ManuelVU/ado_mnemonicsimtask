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
    1. *results$rs*: Array with participants responses organized in 3 dimentions, trials (192), participants by age group (20) and age group (1 = Young, 2 = Elderly). Responses are either 0 = old or 1 = new.
    2. *results$st*: Array with participants stimulus organized in 3 dimentions, trials (192), participants by age group (20) and age group (1 = Young, 2 = Elderly). Stimulus take the values 1 = old, 2 = lure 1, 3 = lure 2, 4 = lure 3, 5 = lure 4, 6 = lure 5 and 7 = new.

* __signaldt.RData__: File contains the results of the analysis using the non-contaminant model. Results are stored in a list with multiple entries with two main elements, data and optimal:
    1. *sdt$data*: Contains the results of the analysis using the original stimulus sequence in the experiment in a trial by trial basis. The results are sotered as follows.
        - *sdt$data$mean*: Array with posterior mean of the d' and k parameters. This array is divided in 4 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$vcov*: Array with the estimated posterior variance covariance matrix of the parameters in a trial by trial basis. This array is divided in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$ci*: Array with the posterior credible interval of each parameter in a trial by trial basis. This array is organized in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), quantile (2.5%,97.5%, 0.5%, 99.5%, 5%, 95%), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$ut*: Array with the KL duvergence between prior and posterior distributions on a trial by trial basis. This array is divided in 3 dimentions, trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
