## replicating models 1-5 from Marnie and doing prior and posterior predictive checks
## fixing spec and sensitivity as a value
## Feb 2023
library(here)

# model 1 -----------------------------------------------------------------
file_m1 <- file.path(("code/experiments/exp1/model1_fixedspecsens.stan"))
mod_m1 <- cmdstan_model(file_m1)

## sampling ####
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
file_m2 <- file.path(("code/experiments/exp1/model2_fixedspecsens.stan"))
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
file_m3 <- file.path(("code/experiments/exp1/model3_fixedspecsens.stan"))
mod_m3 <- cmdstan_model(file_m3)

## sampling ####
## posterior pred
data_list$prior_only = 0

res_posterior_m3 <- mod_m3$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 4 - without agegp -------------------------------------------------
file_m4 <- file.path(("code/experiments/exp1/model4_fixedspecsens.stan"))
mod_m4 <- cmdstan_model(file_m4)

## sampling ####
# post pred
data_list$prior_only = 0

res_posterior_m4 <- mod_m4$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 5 - adapting Marnie's code to do prior predictive checks --------------------
file_m5 <- file.path(("code/experiments/exp1/model5_fixedspecsens.stan"))
mod_m5 <- cmdstan_model(file_m5)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m5 <- mod_m5$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 6 - with all two-way interactions between SEIFA, sex, agegp  --------------------
file_m6 <- file.path(("code/experiments/exp1/model6_fixedspecsens.stan"))
mod_m6 <- cmdstan_model(file_m6)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m6 <- mod_m6$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_posterior_m1$summary("p_avg1")
res_posterior_m2$summary("p_avg1")
res_posterior_m3$summary("p_avg1")
res_posterior_m4$summary("p_avg1")
res_posterior_m5$summary("p_avg1")
res_posterior_m6$summary("p_avg1")

save(res_posterior_m1, res_posterior_m2,
     res_posterior_m3, res_posterior_m4, 
     res_posterior_m5, res_posterior_m6, file=here("results/experiments/exp3/fixedspecsens.RData"))
