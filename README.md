Adaptive Design Optimization for a Mnemonic Similarity Task
===========================================================

This repository contains the results from the "Adaptive Design Optimization for a Mnemonic Similatiry Task" working paper. In order to be able to use the files in this repository one needs to open __ado_mnemonicsimtask.Rproj__ as a project in R, this is needed to make the paths used by functions and code to be the correct ones. Files are divided in 4 main folders:

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
        - *sdt$data$ut*: Array with the KL divergence between prior and posterior distributions on a trial by trial basis. This array is divided in 3 dimentions, trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
    2. *sdt$optimal*: Contains the results of the ADO analysis using the non-contaminant model restricting the number of stimulus used to the ones in the original experiment. The results are stored as follows:
        - *sdt$optimal$mean*: Array with posterior mean of the d' and k parameters. This array is divided in 4 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$vcov*: Array with the estimated posterior variance covariance matrix of the parameters in a trial by trial basis. This array is divided in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$ci*: Array with the posterior credible interval of each parameter in a trial by trial basis. This array is organized in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), quantile (2.5%,97.5%, 0.5%, 99.5%, 5%, 95%), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$ut*: Array with the KL divergence between prior and posterior distributions on a trial by trial basis. This array is divided in 3 dimentions, trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$expdes*: Array with the order in which the trials of a participant where presented in the ADO (in trial number format). This array is organized in 3 dimentions, trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$storder*: Array with the order of the KL divergence by stimulus type, highest expected KL in first position. The array is organized in 4 dimentions, trials (192), KL value (in order as function of KL, 7), participant (20) and age group (1 = Young, 2 = Elderly).  

* __contaminant_sdt.RData__: File contains the results of the analysis using the **contaminant model**. Results are stored in a list with multiple entries with four main elements, data, optimal, pdred and group:
    1. *sdt$data*: Contains the results of the analysis using the original stimulus sequence in the experiment in a trial by trial basis. The results are sotered as follows.
        - *sdt$data$mean*: Array with posterior mean of the d' and k parameters. This array is divided in 4 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$vcov*: Array with the estimated posterior variance covariance matrix of the parameters in a trial by trial basis. This array is divided in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$ci*: Array with the posterior credible interval of each parameter in a trial by trial basis. This array is organized in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), quantile (2.5%,97.5%, 0.5%, 99.5%, 5%, 95%), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$ut*: Array with the KL divergence between prior and posterior distributions on a trial by trial basis. This array is divided in 3 dimentions, trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
    2. *sdt$optimal*: Contains the results of the ADO analysis using the non-contaminant model restricting the number of stimulus used to the ones in the original experiment. The results are stored as follows:
        - *sdt$optimal$mean*: Array with posterior mean of the d' and k parameters. This array is divided in 4 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$vcov*: Array with the estimated posterior variance covariance matrix of the parameters in a trial by trial basis. This array is divided in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$ci*: Array with the posterior credible interval of each parameter in a trial by trial basis. This array is organized in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), quantile (2.5%,97.5%, 0.5%, 99.5%, 5%, 95%), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$ut*: Array with the KL divergence between prior and posterior distributions on a trial by trial basis. This array is divided in 3 dimentions, trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$expdes*: Array with the order in which the trials of a participant where presented in the ADO (in trial number format). This array is organized in 3 dimentions, trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$storder*: Array with the order of the KL divergence by stimulus type, highest expected KL in first position. The array is organized in 4 dimentions, trials (192), KL value (in order as function of KL, 7), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$kl*: Array with the expected KL divergence by stimulus type in a trial bytrial basis. This array is divided in 4 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$group*: Array that contains the posterior mean of the contaminant indicator variable (z) in the contaminant model. This array is organized in 4 dimentions, trials (192), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
    3. *sdt$pred*: Array that contains the posterior predictive distribution of the probability of a "yes" response in a trial by trail basis for the original experiment ordering. Results are stored as follows:
        - *sdt$pred$mean*: Array with posterior mean of the probability of a "yes" response in a trial by trial basis. This array is organized in 3 dimentions, trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
        - *sdt$pred$ci*: Array with the posterior credible interval of the probability of a "yes" response by trial. The array is organized in 4 dimentions, quantile (2.5%,97.5%, 0.5%, 99.5%, 5%, 95%), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).
    4. *sdt$group*: Array with the posterior mean of the probability of contaminant in a trial by trial basis for the original experimental design. This array is organized in 4 dimentions, trials (192), trials (192), participant (20) and age group (1 = Young, 2 = Elderly).

* __simulation_multiple_designs.Rdata__: File contains the results of the simulations using the **contaminant** model. Results are stored in a list with two main elements, data and optimal:
    1. *sdt$data*: Contains the results of the analysis using simulations with 5 different experimental designs used in the original experiment for each simulated participant. Results are sotered as follows.
        - *sdt$data$mean*: Array with posterior mean of the d' and k parameters. This array is divided in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), designs (5),  participant (3) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$vcov*: Array with the estimated posterior variance covariance matrix of the parameters in a trial by trial basis. This array is divided in 6 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), designs(5), participant (3) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$ci*: Array with the posterior credible interval of each parameter in a trial by trial basis. This array is organized in 6 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), quantile (2.5%,97.5%, 0.5%, 99.5%, 5%, 95%), trials (192), design (5), participant (3) and age group (1 = Young, 2 = Elderly).
        - *sdt$data$ut*: Array with the KL divergence between prior and posterior distributions on a trial by trial basis. This array is divided in 4 dimentions, trials (192), design(5), participant (3) and age group (1 = Young, 2 = Elderly).
    2. *sdt$optimal*: Contains the results of the analysis using simulations with 3 participants using the ADO procedure. Results are sotered as follows.
        - *sdt$optimal$mean*: Array with posterior mean of the d' and k parameters. This array is divided in 4 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (3) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$vcov*: Array with the estimated posterior variance covariance matrix of the parameters in a trial by trial basis. This array is divided in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (3) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$ci*: Array with the posterior credible interval of each parameter in a trial by trial basis. This array is organized in 5 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), quantile (2.5%,97.5%, 0.5%, 99.5%, 5%, 95%), trials (192), participant (3) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$ut*: Array with the KL divergence between prior and posterior distributions on a trial by trial basis. This array is divided in 3 dimentions, trials (192), participant (3) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$expdes*: Array with the stimulus presented in the ADO simulation. This array is organized in 3 dimentions, trials (192), participant (3) and age group (1 = Young, 2 = Elderly).
        - *sdt$optimal$kl*: Array with the expected KL divergence by stimulus type in a trial bytrial basis. This array is divided in 4 dimentions, parameter (k, d'<sub>1</sub>, ..., d'<sub>new</sub>), trials (192), participant (3) and age group (1 = Young, 2 = Elderly).

# SRC files

Directory __src__ contains contains 5 types of files dedicated to data analysis, storing functions, drawing figures, writting Bayesian graphical models and data parsing. These files are:

* __Data analysis__: Files in this category contain code to implement the analysis using the non-contaminant, contaminant model and simulations. These filaes make use of functions in the __Functions.R__ file in the same directory. The files in this category are.
    1. **ado_noncontaminant.R**: Data analysis using the non-contaminant model. The analysis is carried out on the original sequence of stimulus presented to the participant and with the ADO ordering independently. Results of this analysis are stored in the file __signaldt.RData__.
    2. **ado_contaminant.R**: Data analysis using the contaminant model. The analysis is carried out on the original sequence of stimulus and on the ADO ordering independently. Results of this analysis are stored in the file __contaminant_sdt.RData__.
    3. **simulation.R**: This file simulates data usingt the non-contaminant model and does a data analysis using the contaminant model. Results from this analysis are stored on the file __simulation_multiple_designs.Rdata__.

* __Bayesian graphical models__: This files write the Bayesian graphical models used by different functions using a "*.txt*" format. The files in this category are:
    1. **write_noncontaminant.R**: Writtes non-contaminant Bayesian model used to expand the data for ADO use. The result of the code is stored on the file __noncontaminant.txt__.
    2. **write_noncontaminant_individual_pred.R**: Writes non-contaminant Bayesian model used to analyze the data at the individual level. The result of the code is stored on the file __noncontaminant_individual_pred.txt__.
    3. **write_contaminant.R**: Writtes contaminant Bayesian model used to expand the data for ADO use. The result of the code is stored on the file __contaminant.txt__.
    4. **write_contaminant_individual_pred.R**: Writes contaminant Bayesian model used to analyze the data at the individual level. The result of the code is stored on the file __contaminant_individual_pred.txt__.

* __Graphs__: Code used to generate the figures in the working paper and appendix. The files in this category are:
    1. **plots.R**: Contains the code to produce Figures 1, 4, 5, 6, 7, 8 and 9. Also contains code to generate graphs found in the appendix section.
    2. **drawfigures_1.m**: Contains code to produce Figure 2 and Figure 10. 

* __data_t.R__: The code in this file is used to parse the raw data form the experiment found i the dile __StarkExp1.mat__. Output of the code is saved in the file __memory.Rdata__. 

* __Functions.R__: This file contains the functions used for data analysis and in the ADO process along with two functions that allow to simulate behavior and calculate the KL divergence between two multivariate normal distributions. Files in the __Data analysis__ section use these functions.

* __PntoneFall2016.mat__: File contains colors used in the figures in rgb format.

* __pantoneColors.mat__: File contains colors used in figures 2 and 10 in rgb format.

# Models

This directory contains the Bayesian graphical models written for use in JAGS as ".txt" files. The files are:
* __noncontaminant.txt__: Code for the non-contaminant model that expands data in order to apply ADO method.
* __noncontaminant_individual_pred.txt__: Code for the non-contaminant model. Can be used to generate posterior predictive distribution of the probability of a "yes" response.
* __contaminant.txt__: Code for the contaminant model that expands data in order to apply ADO method.
* __contaminant_individual_pred.txt__: Code for the contaminant model. Can be used to generate posterior predictive distribution of the probability of a "yes" response.

# Figures 

This directory contains the Figures from the working paper along with the ones used as examples and in the appendix section. There are also graphs at the individual level that where not used.
