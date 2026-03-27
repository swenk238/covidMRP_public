# covidMRP
## Replicating and investigating the issue Marnie Downes et. al. encountered in [Machalek et al. (2022)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0265858)
The issue was that the MRP estimate decreases (away from the crude value of 1.6%) when additional covariates are added to the MRP model that also estimates sensitivity/specificity, suggesting model misfit.

### File Structure

```
├── covidMRP.Rproj                            <- R project to run the code below
├── code                                      <- stan code for all Marnie's original models
│   ├── alternative_figures_realdata.Rmd      <- Rmd file to make all images for comparison of real data for models used in experiments
│   ├── alternative_figures_realdata.html     <- knitted html of Rmd
│   ├── figures_simdata.Rmd                   <- creates the figures for the experiments with simulated data
│   ├── figures_simdata.html                  <- knitted html of Rmd
│   ├── prior_predictive_checks_model.R       <- creates prior preditives for each model, used in Fig 13
│   ├── real_model_posterior.R 
│   ├── example_plots
│   │   ├── effect_sens_spec.R                <- difference between obs and true prevalence as sens and spec varied (Fig 2)
│   │   ├── exploring_priors_intercept_change.R<- Short simulation exploring impact of models changing the prior (not included in paper)
│   │   ├── exploring_priors_intercept_change.R<- Short simulation exploring impact of models changing the with default priors (not included in paper)
│   ├── realdata
│   │   ├── archived                          <- archived folder with previous code development and testing implications of SBC and R2D2M2
│   │   │   ├── data_list.R                   <- list of data for running stan code (embedded in models.R)
│   │   │   ├── R2D2M2
│   │   │   ├── SBC
│   │   │   └── test
│   │   ├── functions.R
│   │   ├── models.R                          <- runs and saves all the models using Marnie's code 
│   │   ├── sim1_noFE_noME
│   │   │   ├── models_nopsamp.R              <- removing measurement error
│   │   │   ├── og_noFE_noME.R                <- removing measurement error and fixed effects
│   │   │   └── stan                          <- stan code for sim1
│   │   ├── sim2_spec                         
│   │   │   ├── models_estsens.R              <- estimating sensitivity and specificity
│   │   │   ├── models_fixedsens.R            <- fixing sensitivity
│   │   │   ├── models_fixedspecsens.R        <- fixing specificity and sensitivity
│   │   │   └── stan                          <- stan code for sim2
│   │   ├── sim3_noFE  
│   │   │   ├── seifa_investigation.R         <- testing different coding for fixed and varying effect (not included in paper)
│   │   │   └── stan                          <- stan code for sim3
│   │   ├── stan                              <- stan code for original data 
│   ├── simdata
│   │   ├── archived                          <- archived folder (not included in paper) containing previous testing for sim2 and 3 on different prevalence, alternative priors, R2D2M2
│   │   │   ├── sim2_spec                     <- testing p = 0.1, 0.5 
│   │   │   └── sim3_noFE                     <- R2D2M2 and alt prior
│   │   ├── func                              <- functions to help read in results and fixed seed numbers
│   │   ├── sim1_noFE_noME                   
│   │   │   ├── sim1_covariates_extra.R       <- rerunning iterations that failed on cluster
│   │   │   ├── sim1_covariates.R             <- modeling covariates 
│   │   │   └── sim1_intercept_only.R         <- intercept-only models
│   │   ├── sim2_spec                         
│   │   │   ├── data                          <- generating population data
│   │   │   │   └── popn_gen_TNFP_stableN.R   <- setting p = 0.01
│   │   │   ├── sim2_estSpec.R               <- estimating specificity in stan
│   │   │   ├── sim2_noEstSpec_fixedspec.R    <- fixed sens and spec values in stan
│   │   │   ├── sim2_noEstSpec.R              <- not estimating measurement error using brm
│   │   │   └── stan                          <- stan code for sim 2
│   │   │   │   ├── estspec_noncent          
│   │   │   │   └── noEstSpec
│   │   └── sim3_fixedeffects                 <- contains RDS files
│   │   │   ├── adding_two_FE.R               <- adding fixed effects
│   │   │   └── stan                          <- stan code for sim 3
├── results
│   ├── results_list.RData                    <- real data results (not in repo due to size)
│   └── results_sim_list.RData                <- simulation results (not in repo due to size)
```

### Files from Marnie 
Files sent from Marnie including estimates for model 1-7 using two different priors

    - MODEL 1: Sampling strata
    - MODEL 2: Strata, postcode
    - MODEL 3: Strata, postcode, SEIFA
    - MODEL 4: Strata, postcode, SEIFA, sex
    - MODEL 5: Strata, postcode, SEIFA, sex, agegp
    - MODEL 6: MODEL 5 + all two-way interactions between SEIFA, sex, agegp
    - MODEL 7: MODEL 6 + two-way interactions of SEIFA,sex,agegp with sampling strata

