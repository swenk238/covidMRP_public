## replicating models 1-5 from Marnie and doing prior and posterior predictive checks
## removing p_samp to see if the problem replicates
## Feb 2023
library(here)
source(here("code/realdata/data_list.R"), echo=F)
stan_wd <- "code/realdata/sim1_noFE_noME/stan/no_p_samp/"
res_wd <- "results/realdata/sim1_noFE_noME/no_p_samp/"

# model 0 -----------------------------------------------------------------
file_m0 <- file.path(here(paste0(stan_wd,"model0_nopsamp.stan")))
mod_m0 <- cmdstan_model(file_m0)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_m0 <- mod_m0$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 1 -----------------------------------------------------------------
file_m1 <- file.path(here(paste0(stan_wd,"model1_nopsamp.stan")))
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
file_m2 <- file.path(here(paste0(stan_wd,"model2_nopsamp.stan")))
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
file_m3 <- file.path(here(paste0(stan_wd,"model3_nopsamp.stan")))
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
file_m4 <- file.path(here(paste0(stan_wd,"model4_nopsamp.stan")))
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
file_m5 <- file.path(here(paste0(stan_wd,"model5_nopsamp.stan")))
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
file_m6 <- file.path(here(paste0(stan_wd,"model6_nopsamp.stan")))
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

res_posterior_m1$print("p_avg1", digits=5)
res_posterior_m2$print("p_avg1", digits=5)
res_posterior_m3$print("p_avg1", digits=5)
res_posterior_m4$print("p_avg1", digits=5)
res_posterior_m5$print("p_avg1", digits=5)
res_posterior_m6$print("p_avg1", digits=5)

res_posterior_m0$save_object(file=here(paste0(res_wd, "m0nopsamp.RDS")))
res_posterior_m1$save_object(file=here(paste0(res_wd, "m1nopsamp.RDS")))
res_posterior_m2$save_object(file=here(paste0(res_wd, "m2nopsamp.RDS")))
res_posterior_m3$save_object(file=here(paste0(res_wd, "m3nopsamp.RDS")))
res_posterior_m4$save_object(file=here(paste0(res_wd, "m4nopsamp.RDS")))
res_posterior_m5$save_object(file=here(paste0(res_wd, "m5nopsamp.RDS")))
res_posterior_m6$save_object(file=here(paste0(res_wd, "m6nopsamp.RDS")))
