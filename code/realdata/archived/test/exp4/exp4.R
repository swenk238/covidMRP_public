# fitting simulation model to real data 

library(here)
library(ggplot2)
library(tidyverse)
library(arm)
library(cmdstanr)
library(brms)

source(here("code/experiments/data_list.R"), echo=F)
options(mc.cores = 8)
stan_wd <- 'code/simdata/exp4/estSpec/stan/'
stan_og_wd <- 'code/experiments/stan/'
func_wd <- 'code/experiments/'
res_wd <- 'results/experiments/exp4/'

# loading functions to store results
source(here(paste0(func_wd,'functions.R')))

## renaming original data variables to match simulation x1-x5 
dat_renamed <- dat |>
  rename(x1 = agegp,
         x2 = postcode, 
         x3 = strata, 
         x4 = seifa, 
         y = positive)

## sim model brms model to extract stan data 
fit5.brm <- brm(y ~ 0 + Intercept + (1|x1) + (1|x2) + (1|x3) + (1|x4), data=dat_renamed, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.99, chains=1)
dat_stan <- standata(fit5.brm)

dat <- append(dat_stan, list(sens=1, 
                             spec = 1-(3/797),
                             TN = 797, 
                             FP = 3))


# Bayesian mixed effects model - stan code --------------------------------
# model 0 -----------------------------------------------------------------
file_m0 <- file.path(here(paste0(stan_wd,"m0_estspec.stan")))
mod_m0 <- cmdstan_model(file_m0)

res_m0 <- mod_m0$sample(
  data = dat,
  seed = 2468,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 1 -----------------------------------------------------------------
file_m1 <- file.path(here(paste0(stan_wd, "m1_estspec.stan")))
mod_m1 <- cmdstan_model(file_m1)

res_m1 <- mod_m1$sample(
  data = dat,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 2 -----------------------------------------------------------------
file_m2 <- file.path(here(paste0(stan_wd,"m2_estspec.stan")))
mod_m2 <- cmdstan_model(file_m2)

res_m2 <- mod_m2$sample(
  data = dat,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 3 -----------------------------------------------------------------
file_m3 <- file.path(here(paste0(stan_wd,"m3_estspec.stan")))
mod_m3 <- cmdstan_model(file_m3)

res_m3 <- mod_m3$sample(
  data = dat,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 4 -----------------------
file_m4 <- file.path(here(paste0(stan_wd,"m4_estspec.stan")))
mod_m4 <- cmdstan_model(file_m4)


res_m4 <- mod_m4$sample(
  data = dat,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# collating results -------------------------------------------------------
pred <- lapply(list(res_m0, res_m1, res_m2, 
                    res_m3, res_m4), res_stan_func, sample_data = dat, 
               spec = dat$spec, TN = dat$TN, FP = dat$FP)

res_tb <- do.call(rbind, pred) %>% 
  as_tibble() %>% 
  mutate(model = c(0:4))


res_tb |> 
  ggplot(aes(x=model)) +
  geom_errorbar(aes(ymin = pred_mean_q5, ymax = pred_mean_q95), col="grey") +
  geom_point(aes(y=pred_mean_q50), col='blue', alpha=0.7) + 
  geom_point(aes(y=est_int), col='orange', alpha=0.7) +
  ylim(c(0.011, 0.02)) +
  theme_bw()


save.image(file=paste0(res_wd, 'simmodel_realdata.RData'))


