## replicating models 1-5 from Marnie and doing prior and posterior predictive checks
## fixing spec and sensitivity as a value
## Feb 2023
library(cmdstanr)
library(here)
source(here("code/realdata/data_list.R"), echo=F)
stan_noFE_wd <- "code/realdata/sim2_spec/stan/fixedspecsens_noFE/"
res_wd <- "results/realdata/sim2_spec/"

data_list$sens = 1
spec_vec = c(0.98, 0.99, 0.995, 1)
for(i in 1:length(spec_vec)){
  data_list$spec = spec_vec[i]
  spec = data_list$spec
  
  # model 0 -----------------------------------------------------------------
  file_m0 <- file.path(here(paste0(stan_noFE_wd, "model0_fixedspecsens.stan")))
  mod_m0 <- cmdstan_model(file_m0)
  
  ## sampling ####
  res_posterior_m0 <- mod_m0$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )

  # model 1 -----------------------------------------------------------------
  file_m1 <- file.path(here(paste0(stan_noFE_wd, "model1_fixedspecsens.stan")))
  mod_m1 <- cmdstan_model(file_m1)
  
  ## sampling ####
  res_posterior_m1 <- mod_m1$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 2
  file_m2 <- file.path(here(paste0(stan_noFE_wd, "model2_fixedspecsens.stan")))
  mod_m2 <- cmdstan_model(file_m2)
  
  ## sampling ####
  res_posterior_m2 <- mod_m2$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 3 without agegp, sex ----------------------------------------------
  file_m3 <- file.path(here(paste0(stan_noFE_wd, "model3_fixedspecsens_noFE.stan")))
  mod_m3 <- cmdstan_model(file_m3)
  
  ## sampling ####
  res_posterior_m3 <- mod_m3$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 4 - without agegp -------------------------------------------------
  file_m4 <- file.path(here(paste0(stan_noFE_wd, "model4_fixedspecsens_noFE.stan")))
  mod_m4 <- cmdstan_model(file_m4)
  
  ## sampling ####
  res_posterior_m4 <- mod_m4$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 5 - adapting Marnie's code to do prior predictive checks --------------------
  file_m5 <- file.path(here(paste0(stan_noFE_wd, "model5_fixedspecsens_noFE.stan")))
  mod_m5 <- cmdstan_model(file_m5)
  
  ## sampling ####
  res_posterior_m5 <- mod_m5$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  res_posterior_m0$save_object(file=here(paste0(res_wd, "model0spec", spec, "sens1_noFE.RDS")), compress=T)
  res_posterior_m1$save_object(file=here(paste0(res_wd, "model1spec", spec, "sens1_noFE.RDS")), compress=T)
  res_posterior_m2$save_object(file=here(paste0(res_wd, "model2spec", spec, "sens1_noFE.RDS")), compress=T)
  res_posterior_m3$save_object(file=here(paste0(res_wd, "model3spec", spec, "sens1_noFE.RDS")), compress=T)
  res_posterior_m4$save_object(file=here(paste0(res_wd, "model4spec", spec, "sens1_noFE.RDS")), compress=T)
  res_posterior_m5$save_object(file=here(paste0(res_wd, "model5spec", spec, "sens1_noFE.RDS")), compress=T)
}

res_posterior_m1$print("p_avg1", digits=5)
res_posterior_m2$print("p_avg1", digits=5)
res_posterior_m3$print("p_avg1", digits=5)
res_posterior_m4$print("p_avg1", digits=5)
res_posterior_m5$print("p_avg1", digits=5)
res_posterior_m6$print("p_avg1", digits=5)
