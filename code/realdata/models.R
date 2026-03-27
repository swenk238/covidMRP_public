## replicating models 1-5 from Marnie and doing prior and posterior predictive checks
## Jan 2023
library(here)
library(cmdstanr)
source(here("code/realdata/data_list.R"), echo=F)
options(mc.cores = 8)

res_wd <- "results/realdata/original/"
stan_wd <- "code/realdata/stan/"

# model 7 - with interactions with strata --------------------
file_m7 <- file.path(here(paste0(stan_wd, "model7.stan")))
mod_m7 <- cmdstan_model(file_m7)

## sampling ####
res_posterior_m7 <- mod_m7$sample(
  data = data_list,
  seed = 2345,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# prior pred
data_list$prior_only = 1

res_prior_m7 <- mod_m7$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_prior_m7$save_object(file=here(paste0(res_wd, "m7_priorpred.RDS")), compress=T)
res_posterior_m7$save_object(file=here(paste0(res_wd, "m7_postpred.RDS")), compress=T)

# model 6 - with all two-way interactions between SEIFA, sex, agegp  --------------------
file_m6 <- file.path(here(paste0(stan_wd, "model6.stan")))
mod_m6 <- cmdstan_model(file_m6)

## sampling ####
# posterior predictive checks 
res_posterior_m6 <- mod_m6$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# prior pred
data_list$prior_only = 1

res_prior_m6 <- mod_m6$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_prior_m6$save_object(file=here(paste0(res_wd, "m6_priorpred.RDS")), compress=T)
res_posterior_m6$save_object(file = here(paste0(res_wd, "m6_postpred.RDS")), compress=T)

# model 6a - with prior for p  --------------------
file_m6a <- file.path(here(paste0(stan_wd, "model6a.stan")))
mod_m6a <- cmdstan_model(file_m6a)

# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m6a <- mod_m6a$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_posterior_m6$summary("p_avg1")
res_posterior_m6a$summary("p_avg1")

# model 5 - adapting Marnie's code to do prior predictive checks --------------------
file_m5 <- file.path(here(paste0(stan_wd, "model5.stan")))
mod_m5 <- cmdstan_model(file_m5)

## sampling ####
# prior pred
data_list$prior_only = 1

res_prior_m5 <- mod_m5$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m5 <- mod_m5$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_posterior_m5$draws("y_rep[4780]")
res_posterior_m5$summary("a_strata")

# saving
res_prior_m5$save_object(file=here(paste0(res_wd, "m5_priorpred.RDS")), compress=T)
res_posterior_m5$save_object(file=here(paste0(res_wd, "m5_postpred.RDS")), compress=T)


# model 4 - without agegp -------------------------------------------------
file_m4 <- file.path(here(paste0(stan_wd, "model4.stan")))
mod_m4 <- cmdstan_model(file_m4)

## sampling ####
## prior pred
data_list$prior_only = 1

res_prior_m4 <- mod_m4$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
  iter_sampling = 25000,
)

# post pred
data_list$prior_only = 0

res_posterior_m4 <- mod_m4$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_posterior_m4$summary("p_avg1")

# saving
res_prior_m4$save_object(file=here(paste0(res_wd, "m4_priorpred.RDS")), compress=T)
res_posterior_m4$save_object(file=here(paste0(res_wd, "m4_postpred.RDS")), compress=T)


# model 3 without agegp, sex ----------------------------------------------
file_m3 <- file.path(here(paste0(stan_wd, "model3.stan")))
mod_m3 <- cmdstan_model(file_m3)

## sampling ####
## prior pred
data_list$prior_only = 1

res_prior_m3 <- mod_m3$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
  iter_sampling = 25000,
)

## posterior pred
data_list$prior_only = 0

res_posterior_m3 <- mod_m3$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_posterior_m3$summary("b[2]")

# saving
res_prior_m3$save_object(file=here(paste0(res_wd, "m3_priorpred.RDS")), compress=T)
res_posterior_m3$save_object(file=here(paste0(res_wd, "m3_postpred.RDS")), compress=T)


# model 2 without agegp, sex, seifa ----------------------------------------------
file_m2 <- file.path(here(paste0(stan_wd, "model2.stan")))
mod_m2 <- cmdstan_model(file_m2)

## sampling ####
## prior pred
data_list$prior_only = 1

res_prior_m2 <- mod_m2$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
  iter_sampling = 25000,
)

## post pred
data_list$prior_only = 0

res_posterior_m2 <- mod_m2$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_posterior_m2$summary("p_avg1")

# saving
res_prior_m2$save_object(file=here(paste0(res_wd, "m2_priorpred.RDS")), compress=T)
res_posterior_m2$save_object(file=here(paste0(res_wd, "m2_postpred.RDS")), compress=T)


# model 1 without agegp, sex, seifa, postcode ----------------------------------------------
file_m1 <- file.path(here(paste0(stan_wd, "model1.stan")))
mod_m1 <- cmdstan_model(file_m1)

## sampling ####
## prior pred
data_list$prior_only = 1

res_prior_m1 <- mod_m1$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4
)

# posterior predictive
data_list$prior_only = 0

res_posterior_m1 <- mod_m1$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_posterior_m1$summary("p_avg1")

# saving
res_prior_m1$save_object(file=here(paste0(res_wd, "m1_priorpred.RDS")), compress=T)
res_posterior_m1$save_object(file=here(paste0(res_wd, "m1_postpred.RDS")), compress=T)

# model 0 without any covariates ----------------------------------------------
file_m0 <- file.path(here(paste0(stan_wd, "model0.stan")))
mod_m0 <- cmdstan_model(file_m0)

## sampling ####
## prior pred
data_list$prior_only = 1

res_prior_m0 <- mod_m0$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  iter_sampling = 25000,
)

# posterior predictive
data_list$prior_only = 0

res_posterior_m0 <- mod_m0$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_posterior_m0$summary("p_avg1")

# saving
res_prior_m0$save_object(file=here(paste0(res_wd, "m0_priorpred.RDS")), compress=T)
res_posterior_m0$save_object(file=here(paste0(res_wd, "m0_postpred.RDS")), compress=T)

# model A with postcode only  ----------------------------------------------
file_mA <- file.path(here("code/realdata/exp0/modelA.stan"))
mod_mA <- cmdstan_model(file_mA)

## sampling ####
## prior pred
data_list$prior_only = 1

res_prior_mA <- mod_mA$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  iter_sampling = 25000,
)

# posterior predictive
data_list$prior_only = 0

res_posterior_mA <- mod_mA$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_prior_mA$save_object(file=here(paste0(res_wd, "mA_priorpred_FP0TN10.RDS")), compress=T)
res_posterior_mA$save_object(file=here(paste0(res_wd, "mA_postpred.RDS")), compress=T)

# model Z -- without inv_logit  ----------------------------------------------
file_mZ <- file.path(here(paste0(stan_wd, "modelZ.stan")))
mod_mZ <- cmdstan_model(file_mZ)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_mZ <- mod_mZ$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_mZ$save_object(file=here("results/realdata/original/mB_postpred_FP0TN10_priorscale0.5.RDS"), compress=T)


# model B -- model 3 without seifa_cont  ----------------------------------------------
file_mB <- file.path(here(paste0(stan_wd, "modelB.stan")))
mod_mB <- cmdstan_model(file_mB)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_mB <- mod_mB$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_mB$save_object(file=here(paste0(res_wd, "mB_postpred_FP0TN10_priorscale0.5.RDS")), compress=T)

# model C -- model 3 without seifa_cont, postcode (strata and seifa) ----------------------------------------------
file_mC <- file.path(here(paste0(stan_wd, "modelC.stan")))
mod_mC <- cmdstan_model(file_mC)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_mC <- mod_mC$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_mC$save_object(file=here(paste0(res_wd, "mC_postpred_FP0TN10_priorscale0.5.RDS")), compress=T)

# model D -- model 3 without seifa_cont, postcode, strata  ----------------------------------------------
file_mD <- file.path(here(paste0(stan_wd, "modelD.stan")))
mod_mD <- cmdstan_model(file_mD)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_mD <- mod_mD$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_mD$save_object(file=here(paste0(res_wd, "mD_postpred_FP0TN10_priorscale0.5.RDS")), compress=T)
