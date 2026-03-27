## running marnie's model without fixed effects 
library(here)
library(ggplot2)
library(tidyverse)
library(arm)
library(cmdstanr)
library(brms)

source(here("code/realdata/data_list.R"), echo=F)
data_list$prior_only = 0 

options(mc.cores = 8)
stan_noFE_noME_wd <- 'code/realdata/sim1_noFE_noME/stan/'
func_wd <- 'code/realdata/'
res_noFE_wd <- 'results/realdata/sim1_noFE_noME/'

# loading functions to store results
source(here(paste0(func_wd,'functions.R')))

# model 3 without agegp, sex ----------------------------------------------
file_m3a <- file.path(here(paste0(stan_noFE_noME_wd, "model3a_noFE_noME.stan")))
mod_m3a <- cmdstan_model(file_m3a)

res_m3a <- mod_m3a$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)


# model 4 - without agegp -------------------------------------------------
file_m4a <- file.path(here(paste0(stan_noFE_noME_wd, "model4a_noFE_noME.stan")))
mod_m4a <- cmdstan_model(file_m4a)

res_m4a <- mod_m4a$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 5 ---------------------
file_m5a <- file.path(here(paste0(stan_noFE_noME_wd,"model5a_noFE_noME.stan")))
mod_m5a <- cmdstan_model(file_m5a)

res_m5a <- mod_m5a$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving results ----------------------------------------------------------
res_m3a$save_object(file = here(paste0(res_noFE_noME_wd,"res_m3a_noFE_noME.RDS")))
res_m4a$save_object(file = here(paste0(res_noFE_noME_wd,"res_m4a_noFE_noME.RDS")))
res_m5a$save_object(file = here(paste0(res_noFE_noME_wd,"res_m5a_noFE_noME.RDS")))

