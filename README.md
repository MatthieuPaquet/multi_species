# multi_species
This folder contains R code of a two-species integrated predator-prey model in Nimble used to simulate data and fit the models to the simulated data

Data and figures can be found on the OSF project:
https://osf.io/xfa6e/

##removal
Contains a capture-removal version of the IPM (whereby individuals are only recaptured once). More details on the README file of this folder.

### Nimble_PBpred_time30ind100.R
File used to both simulate the data and fit the Integrated Predator Prey Model on them for the results presented in the main text, that is:
-100 individuals of each species are marked every year for 30 years
-Densities were centered

### Nimble_PBpred_samplesizesBG2019.R
File used to both simulate and fit the data with the same data design as in Barraquand & Gimenez 2019, that is either:
-100 individuals of each species are marked for 10 years
or
-20 individuals of each species are marked for 30 years.

### outputs_scenariosmaintext.R
File used to extract estimates and plot figures resulting from the "Nimble_BG2019_centered_T30_N100.R" file

### outputs_samplesizeBG2019_noncentered.R
File used to extract estimates and plot figures resulting from the "Nimble_BG2019_obs_error.R" file and the "Nimble_BG2019_obs_error_simul_nodd_model_dd.R" file.

### script_choose_initial_values.R
Code to choose initial values for assessing sensitivity of the results to the choice of initial values for density dependent parameters, see section "Sensitivity of parameter estimation to the choice of initial values" in the Appendix

### script_simulatedinitial_values_out_of_posterior.R
Code to perform the Sensitivity analysis to the choice of initial values presended in the Appendix.
