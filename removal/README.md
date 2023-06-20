# multi_species capture-removal Integrated population model
In a previous version of this work, we inadvertently omitted the CMR array of adults in the code, which transformed the model into a capture-removal model (i.e., individuals were re-captured only once and then removed from the population, as in hunting or fishing data). 

This folder therefore contains R codes of a two-species integrated predator-prey model in Nimble used to simulate data and fit the models to the simulated data

### Nimble_BG2019_centered_T30_N100.R
File used to both simulate the data and fit the Integrated Predator Prey Model on them for the results presented in the main text, that is:
-100 individuals of each species are marked every year for 30 years
-Densities were centered

### Nimble_BG2019_obs_error.R
File used to simulate the data with the same data design as in Barraquand & Gimenez 2019, that is either:
-100 individuals of each species are marked for 10 years
or
-20 individuals of each species are marked for 30 years.

### Nimble_BG2019_obs_error_simul_nodd_model_dd.R
File used to fit the data simulated with the previous file (Nimble_BG2019_obs_error.R). Note that here densities were not centered to be consistent with Barraquand & Gimenez for comparison.

### outputs_centered_scenariosmaintext.R
File used to extract estimates and plot figures resulting from the "Nimble_BG2019_centered_T30_N100.R" file

### outputs_nodd_noncentered_samplesizeBG2019.R
File used to extract estimates and plot figures resulting from the "Nimble_BG2019_obs_error.R" file and the "Nimble_BG2019_obs_error_simul_nodd_model_dd.R" file.

### script_initial_values.R
Code to choose initial values for assessing sensitivity of the results to the choice of initial values for density dependent parameters, see section "Sensitivity of parameter estimation to the choice of initial values" in the Appendix

### script_MCMC_simulatedinitial_values_out_of_posterior.R
Code to perform the Sensitivity analysis to the choice of initial values presended in the Appendix.
