## running marnie's model without fixed effects 
library(here)
library(ggplot2)
library(tidyverse)
library(arm)
library(cmdstanr)
library(brms)

source(here("code/realdata/data_list.R"), echo=F)

options(mc.cores = 8)
stan_noFE_wd <-  'code/realdata/sim3_noFE/stan/'
func_wd <- 'code/realdata/'
res_wd <- 'results/realdata/sim3_noFE/'
res_noFE_wd <- 'results/realdata/sim3_noFE/'

# loading functions to store results
source(here(paste0(func_wd,'functions.R')))


# model 3a - without agegp, sex, seifa_cont ----------------------------------------------
file_m3a <- file.path(here(paste0(stan_noFE_wd, "model3a_rmSeifaCont.stan")))
mod_m3a <- cmdstan_model(file_m3a)

res_m3a <- mod_m3a$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 3b - without agegp, sex, a_seifa ----------------------------------------------
file_m3b <- file.path(here(paste0(stan_noFE_wd, "model3b_rmASeifa.stan")))
mod_m3b <- cmdstan_model(file_m3b)

res_m3b <- mod_m3b$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 3c - without agegp, sex, factor(seifa) ----------------------------------------------
# for sample
fit3.brm <- brm(y ~ (1|strata) + (1|postcode) + as.factor(seifa), data=dat_fm, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.99, chains=1)
dat3_stan <- standata(fit3.brm)

data_list$seifa2 = dat3_stan$X[,'as.factorseifa2']
data_list$seifa3 = dat3_stan$X[,'as.factorseifa3']
data_list$seifa4 = dat3_stan$X[,'as.factorseifa4']
data_list$seifa5 = dat3_stan$X[,'as.factorseifa5']

data_list$seifa2_pop = X_pop[,'as.factorseifa_pop2']
data_list$seifa3_pop = X_pop[,'as.factorseifa_pop3']
data_list$seifa4_pop = X_pop[,'as.factorseifa_pop4']
data_list$seifa5_pop = X_pop[,'as.factorseifa_pop5']

file_m3c <- file.path(here(paste0(stan_noFE_wd, "model3c_fctSeifa.stan")))
mod_m3c <- cmdstan_model(file_m3c)

res_m3c <- mod_m3c$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_m3c$summary('p_avg1')


# model 4a - without agegp, remove seifa_cont  -------------------------------------------------
file_m4a <- file.path(here(paste0(stan_noFE_wd, "model4a_rmSeifaCont.stan")))
mod_m4a <- cmdstan_model(file_m4a)

res_m4a <- mod_m4a$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 4b - without agegp, remove a_seifa -------------------------------------------------
file_m4b <- file.path(here(paste0(stan_noFE_wd, "model4b_rmASeifa.stan")))
mod_m4b <- cmdstan_model(file_m4b)

res_m4b <- mod_m4b$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 4c - without agegp, factor(seifa) -------------------------------------------------
file_m4c <- file.path(here(paste0(stan_noFE_wd, "model4c_fctSeifa.stan")))
mod_m4c <- cmdstan_model(file_m4c)

res_m4c <- mod_m4c$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 4d - without fixed effects seifa_con, sex -------------------------------------------------
file_m4d <- file.path(here(paste0(stan_noFE_wd, "model4d_rmFE.stan")))
mod_m4d<- cmdstan_model(file_m4d)

res_m4d <- mod_m4d$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_m4d$summary('p_avg1')

# model 4e - without fixed effects seifa_con, sex -------------------------------------------------
file_m4e <- file.path(here(paste0(stan_noFE_wd, "model4e_rmSex.stan")))
mod_m4e<- cmdstan_model(file_m4e)

res_m4e <- mod_m4e$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 5a - remove seifa_cont ---------------------
file_m5a <- file.path(here(paste0(stan_noFE_wd,"model5a_rmSeifaCont.stan")))
mod_m5a <- cmdstan_model(file_m5a)

res_m5a <- mod_m5a$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 5b - rm a_seifa ---------------------
file_m5b <- file.path(here(paste0(stan_noFE_wd,"model5b_rmASeifa.stan")))
mod_m5b <- cmdstan_model(file_m5b)

res_m5b <- mod_m5b$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 5c - factor seifa ---------------------
file_m5c <- file.path(here(paste0(stan_noFE_wd,"model5c_fctSeifa.stan")))
mod_m5c <- cmdstan_model(file_m5c)

res_m5c <- mod_m5c$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 5d - without fixed effects seifa_con, sex----------------------------------------------------------------
file_m5d <- file.path(here(paste0(stan_noFE_wd, "model5d_rmFE.stan")))
mod_m5d <- cmdstan_model(file_m5d)

res_m5d <- mod_m5d$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 5e - without fixed effects for sex----------------------------------------------------------------
file_m5e <- file.path(here(paste0(stan_noFE_wd, "model5e_rmSex.stan")))
mod_m5e <- cmdstan_model(file_m5e)

res_m5e <- mod_m5e$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)


res_m3a$save_object(file=paste0(res_noFE_wd,"res_m3a_rmSeifaCont.RDS"), compress=T)
res_m3b$save_object(file=paste0(res_noFE_wd,"res_m3b_rmASeifa.RDS"))
res_m3c$save_object(file=paste0(res_noFE_wd,"res_m3c_fctSeifa.RDS"))

res_m4a$save_object(file=paste0(res_noFE_wd,"res_m4a_rmSeifaCont.RDS"))
res_m4b$save_object(file=paste0(res_noFE_wd,"res_m4b_rmASeifa.RDS"))
res_m4c$save_object(file=paste0(res_noFE_wd,"res_m4c_fctSeifa.RDS"))
res_m4d$save_object(file=paste0(res_noFE_wd,"res_m4d_rmFE.RDS"))
res_m4e$save_object(file=paste0(res_noFE_wd,"res_m4e_rmSex.RDS"))

res_m5a$save_object(file=paste0(res_noFE_wd,"res_m5a_rmSeifaCont.RDS"))
res_m5b$save_object(file=paste0(res_noFE_wd,"res_m5b_rmASeifa.RDS"))
res_m5c$save_object(file=paste0(res_noFE_wd,"res_m5c_fctSeifa.RDS"))
res_m5d$save_object(file=paste0(res_noFE_wd,"res_m5d_rmFE.RDS"))
res_m5e$save_object(file=paste0(res_noFE_wd,"res_m5e_rmSex.RDS"))

