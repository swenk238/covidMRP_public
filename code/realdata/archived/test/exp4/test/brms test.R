# model 1 without agegp, sex, seifa, postcode ----------------------------------------------
fit1.brm <- brm(y ~ (1|strata), data=dat_fm, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.99, chains=1)
dat1_stan <- standata(fit1.brm)
og_sublist$X_pop1 <- matrix(1, nrow=915)
dat1 <- append(dat1_stan, og_sublist)

file_m1a <- file.path(here(paste0(stan_seifaInv_brms_wd,"test/m1a.stan")))
mod_m1a <- cmdstan_model(file_m1a)

res_m1a <- mod_m1a$sample(
  data = dat1, 
  seed = 2468, 
  chains = 2, 
  parallel_chains = 4,
)

file_m1b <- file.path(here(paste0(stan_seifaInv_brms_wd,"test/m1b.stan")))
mod_m1b <- cmdstan_model(file_m1b)

res_m1b <- mod_m1b$sample(
  data = dat1, 
  seed = 2468, 
  chains = 2, 
  parallel_chains = 4,
)

file_m1c <- file.path(here(paste0(stan_seifaInv_brms_wd,"test/m1c.stan")))
mod_m1c <- cmdstan_model(file_m1c)

res_m1c <- mod_m1c$sample(
  data = dat1, 
  seed = 1357, 
  chains = 2, 
  parallel_chains = 4,
)

file_m1d <- file.path(here(paste0(stan_seifaInv_brms_wd,"test/m1d.stan")))
mod_m1d <- cmdstan_model(file_m1d)

res_m1d <- mod_m1d$sample(
  data = dat1, 
  seed = 1357, 
  chains = 2, 
  parallel_chains = 4,
)

file_m1e <- file.path(here(paste0(stan_seifaInv_brms_wd,"test/m1e.stan")))
mod_m1e <- cmdstan_model(file_m1e)

res_m1e <- mod_m1e$sample(
  data = dat1, 
  seed = 1357, 
  chains = 2, 
  parallel_chains = 4,
)

file_m1s <- file.path(here(paste0(stan_seifaInv_brms_wd,"model1_brms.stan")))
mod_m1s <- cmdstan_model(file_m1s)

res_m1s <- mod_m1s$sample(
  data = dat1, 
  seed = 2468, 
  chains = 2, 
  parallel_chains = 4,
)

res_m1a$summary('b')
res_m1b$summary('b')
res_m1c$summary('b')
res_m1d$summary('b')
res_m1e$summary('b')
res_m1s$summary('b')


res_m1a$summary('p_avg1')
res_m1b$summary('p_avg1')
res_m1c$summary('p_avg1')
res_m1d$summary('p_avg1')
res_m1e$summary('p_avg1')
res_m1s$summary('p_avg1')