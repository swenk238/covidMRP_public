## replicating models 1-5 from Marnie and doing prior and posterior predictive checks
## fixing spec only as a value
## Feb 2023
library(here)
source(here("code/realdata/data_list.R"), echo=F)
stan_noFE_wd <- "code/realdata/sim2_spec/stan/fixedsens_noFE/"
res_wd <- "results/realdata/sim2_spec/"

data_list$sens = 1
data_list$prior_only = 0

total_vec <- rep(c(400,800,1200,8000),4)
FP_vec <- c(8,16,24,160, 4,8,12,80, 2,4,6,40, rep(0,4))
TN_vec <- total_vec - FP_vec # corresponding to specificity values of 0.98, 0.99, 0.995 and 1 

# no fixed effects ####
for(i in 1:length(TN_vec)){
  data_list$TN = TN_vec[i]
  data_list$FP = FP_vec[i]
  
  # model 0 -----------------------------------------------------------------
  file_m0 <- file.path(here(paste0(stan_noFE_wd, "model0_fixedsens.stan")))
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
  file_m1 <- file.path(here(paste0(stan_noFE_wd, "model1_fixedsens.stan")))
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
  file_m2 <- file.path(here(paste0(stan_noFE_wd, "model2_fixedsens.stan")))
  mod_m2 <- cmdstan_model(file_m2)
  
  ## sampling ####
  ## post pred
  data_list$prior_only = 0
  
  res_posterior_m2 <- mod_m2$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 3 without agegp, sex ----------------------------------------------
  file_m3 <- file.path(here(paste0(stan_noFE_wd, "model3_fixedsens_noFE.stan")))
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
  file_m4 <- file.path(here(paste0(stan_noFE_wd, "model4_fixedsens_noFE.stan")))
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
  file_m5 <- file.path(here(paste0(stan_noFE_wd, "model5_fixedsens_noFE.stan")))
  mod_m5 <- cmdstan_model(file_m5)
  
  ## sampling ####
  res_posterior_m5 <- mod_m5$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  res_posterior_m0$save_object(file=here(paste0(res_wd, "m0_TN", TN_vec[i], "FP", FP_vec[i],"fixedsens_noFE.RDS")), compress=T)
  res_posterior_m1$save_object(file=here(paste0(res_wd, "m1_TN", TN_vec[i], "FP", FP_vec[i],"fixedsens_noFE.RDS")), compress=T)
  res_posterior_m2$save_object(file=here(paste0(res_wd, "m2_TN", TN_vec[i], "FP", FP_vec[i],"fixedsens_noFE.RDS")), compress=T)
  res_posterior_m3$save_object(file=here(paste0(res_wd, "m3_TN", TN_vec[i], "FP", FP_vec[i],"fixedsens_noFE.RDS")), compress=T)
  res_posterior_m4$save_object(file=here(paste0(res_wd, "m4_TN", TN_vec[i], "FP", FP_vec[i],"fixedsens_noFE.RDS")), compress=T)
  res_posterior_m5$save_object(file=here(paste0(res_wd, "m5_TN", TN_vec[i], "FP", FP_vec[i],"fixedsens_noFE.RDS")), compress=T)
}
