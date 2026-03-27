## replicating models 1-5 from Marnie and doing prior and posterior predictive checks
## fixing spec only as a value
## July 2023
library(here)
source(("data_list.R"), echo=F)

options(mc.cores = 4)

## spec = 0.99 ####
data_list$spec # should be 0.99

spec_values = c(0.9, 0.92, 0.95, 0.98, 0.99, 0.995, 0.999)
for (i in 1:3){
  data_list$spec = spec_values[i]
  # model 1 
  file_m1 <- file.path(("stan/model1_fixedspec.stan"))
  mod_m1 <- cmdstan_model(file_m1)
  
  ## sampling ###
  # posterior predictive
  data_list$prior_only = 0
  
  res_posterior_m1 <- mod_m1$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 2
  file_m2 <- file.path(("stan/model2_fixedspec.stan"))
  mod_m2 <- cmdstan_model(file_m2)
  
  ## sampling 
  ## post pred
  data_list$prior_only = 0
  
  res_posterior_m2 <- mod_m2$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 3 without agegp, sex -
  file_m3 <- file.path(("stan/model3_fixedspec.stan"))
  mod_m3 <- cmdstan_model(file_m3)
  
  ## sampling 
  ## posterior pred
  data_list$prior_only = 0
  
  res_posterior_m3 <- mod_m3$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 4 - without agegp -
  file_m4 <- file.path(("stan/model4_fixedspec.stan"))
  mod_m4 <- cmdstan_model(file_m4)
  
  ## sampling
  # post pred
  data_list$prior_only = 0
  
  res_posterior_m4 <- mod_m4$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 5 - adapting Marnie's code to do prior predictive checks 
  file_m5 <- file.path(("stan/model5_fixedspec.stan"))
  mod_m5 <- cmdstan_model(file_m5)
  
  ## sampling 
  # posterior predictive checks 
  data_list$prior_only = 0
  
  res_posterior_m5 <- mod_m5$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  # model 6 - with all two-way interactions between SEIFA, sex, agegp  
  file_m6 <- file.path(("stan/model6_fixedspec.stan"))
  mod_m6 <- cmdstan_model(file_m6)
  
  ## sampling 
  # posterior predictive checks 
  data_list$prior_only = 0
  
  res_posterior_m6 <- mod_m6$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  
  res_posterior_m1$save_object(file=here(paste0("m1fixedspec", spec_values[i], ".RDS")))
  res_posterior_m2$save_object(file=here(paste0("m2fixedspec", spec_values[i], ".RDS")))
  res_posterior_m3$save_object(file=here(paste0("m3fixedspec", spec_values[i], ".RDS")))
  res_posterior_m4$save_object(file=here(paste0("m4fixedspec", spec_values[i], ".RDS")))
  res_posterior_m5$save_object(file=here(paste0("m5fixedspec", spec_values[i], ".RDS")))
  res_posterior_m6$save_object(file=here(paste0("m6fixedspec", spec_values[i], ".RDS")))
}

res_posterior_m1$print("p_avg1", digits=5)
res_posterior_m2$print("p_avg1", digits=5)
res_posterior_m3$print("p_avg1", digits=5)
res_posterior_m4$print("p_avg1", digits=5)
res_posterior_m5$print("p_avg1", digits=5)
res_posterior_m6$print("p_avg1", digits=5)


