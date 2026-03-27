data_list_v2 = append(data_list, list(K = 3, R2D2_mean_R2 = 0.5, R2D2_prec_R2 = 2,
                                      Kscales = 5, R2D2_cons_D2= rep(0.5, 5)))

data_list_v3 = append(data_list, list(K = 3, R2D2_mean_R2 = 0.3, R2D2_prec_R2 = 5,
                                      Kscales = 5, R2D2_cons_D2= rep(1, 5)))

stan_wd <- "code/realdata/stan/"
res_wd <- "results/realdata/original/"

# model 1 without agegp, sex, seifa, postcode ----------------------------------------------
file_m1_r2prior <- file.path(here(paste0(stan_wd, "model1_R2D2M2.stan")))
mod_m1_r2prior <- cmdstan_model(file_m1_r2prior)

res_posterior_m1 <- mod_m1_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_posterior_m1$summary("p_avg1")

res_posterior_m1$save_object(file=here(paste0(res_wd, "m1_r2prior.RDS")), compress=T)

res_m1_ori <- readRDS(file=here(paste0(res_wd, "m1_postpred.RDS")))
res_m1_ori$summary("p_avg1")

# prior predictive
data_list_v2$prior_only = 1
res_priorpred_m1_r2prior <- mod_m1_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_priorpred_m1_r2prior$summary("p_avg1")

res_priorpred_m1_r2prior$save_object(file=here(paste0(res_wd, "m1_r2prior_priorpred.RDS")), compress=T)

# model 2 without agegp, sex, seifa ----------------------------------------------
file_m2_r2prior <- file.path(here(paste0(stan_wd, "model2_R2D2M2.stan")))
mod_m2_r2prior <- cmdstan_model(file_m2_r2prior)

# posterior pred
data_list_v2$prior_only = 0
res_m2_r2prior <- mod_m2_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_m2_r2prior$save_object(file=here(paste0(res_wd, "m2_r2prior.RDS")), compress=T)

res_m2_ori <- readRDS(file=here(paste0(res_wd, "m2_postpred.RDS")))
res_m2_ori$summary("p_avg1")
res_m2_r2prior$summary("p_avg1")

# prior predictive
data_list$prior_only = 1
res_priorpred_m2_r2prior <- mod_m2_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_priorpred_m2_r2prior$save_object(file=here(paste0(res_wd, "m2_r2prior_priorpred.RDS")), compress=T)


# model 3 without agegp, sex ----------------------------------------------
file_m3_r2prior <- file.path(here(paste0(stan_wd, "model3_R2D2M2.stan")))
mod_m3_r2prior <- cmdstan_model(file_m3_r2prior)

# posterior predictive
data_list$prior_only = 0
res_m3_r2prior <- mod_m3_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_m3_r2prior$save_object(file=here(paste0(res_wd, "m3_r2prior.RDS")), compress=T)

res_m3_ori <- readRDS(file=here(paste0(res_wd, "m3_postpred.RDS")))
res_m3_ori$summary("p_avg1")
res_m3_r2prior$summary("p_avg1")

# prior predictive ####
data_list$prior_only = 1
res_priorpred_m3_r2prior <- mod_m3_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_priorpred_m3_r2prior$summary("p_avg1")
res_priorpred_m3_r2prior$save_object(file=here(paste0(res_wd, "m3_r2prior_priorpred.RDS")), compress=T)



# model 4 - without agegp -------------------------------------------------
file_m4_r2prior <- file.path(here(paste0(stan_wd, "model4_R2D2M2.stan")))
mod_m4_r2prior <- cmdstan_model(file_m4_r2prior)

# posterior predictive
data_list$prior_only = 0
res_m4_r2prior <- mod_m4_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_m4_r2prior$save_object(file=here(paste0(res_wd, "m4_r2prior.RDS")), compress=T)

res_m4_ori <- readRDS(file=here(paste0(res_wd, "m4_postpred.RDS")))
res_m4_ori$summary("p_avg1")
res_m4_r2prior$summary("p_avg1")

# prior predictive ####
data_list$prior_only = 1
res_priorpred_m4_r2prior <- mod_m4_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_priorpred_m4_r2prior$save_object(file=here(paste0(res_wd, "m4_r2prior_priorpred.RDS")), compress=T)

# model 5 - STRATA, POSTCODE, SEIFA, SEX, AGEGP --------------------
file_m5_r2prior <- file.path(here(paste0(stan_wd, "model5_R2D2M2.stan")))
mod_m5_r2prior <- cmdstan_model(file_m5_r2prior)

# posterior predictive
data_list$prior_only = 0

## sampling ####
res_m5_r2prior <- mod_m5_r2prior$sample(
  data = data_list_v2,
  seed = 2345,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_m5_r2prior$save_object(file=here(paste0(res_wd, "m5_r2prior.RDS")), compress=T)

res_m5_ori <- readRDS(file=here(paste0(res_wd, "m5_postpred.RDS")))
res_m5_ori$summary("p_avg1")
res_m5_r2prior$summary("p_avg1")

# prior predictive ####
data_list$prior_only = 1
res_priorpred_m5_r2prior <- mod_m5_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_priorpred_m5_r2prior$save_object(file=here(paste0(res_wd, "m5_r2prior_priorpred.RDS")), compress=T)


# model 6 - with all two-way interactions between SEIFA, sex, agegp  --------------------
file_m6_r2prior <- file.path(here(paste0(stan_wd, "model6_R2D2M2.stan")))
mod_m6_r2prior <- cmdstan_model(file_m6_r2prior)

# posterior predictive
data_list$prior_only = 0
res_m6_r2prior <- mod_m6_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_m6_r2prior$save_object(file=here(paste0(res_wd, "m6_r2prior.RDS")), compress=T)

res_m6_ori <- readRDS(file=here(paste0(res_wd, "m6_postpred.RDS")))
res_m6_ori$summary("p_avg1")
res_m6_r2prior$summary("p_avg1")

# prior predictive ####
data_list$prior_only = 1
res_priorpred_m6_r2prior <- mod_m6_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_priorpred_m6_r2prior$save_object(file=here(paste0(res_wd, "m6_r2prior_priorpred.RDS")), compress=T)


# model 7 - with interactions with strata --------------------
file_m7_r2prior <- file.path(here(paste0(stan_wd, "model7_R2D2M2.stan")))
mod_m7_r2prior <- cmdstan_model(file_m7_r2prior)

# posterior predictive
data_list$prior_only = 0
## sampling ####
res_posterior_m7_r2prior <- mod_m7_r2prior$sample(
  data = data_list,
  seed = 2345,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_m7_r2prior$save_object(file=here(paste0(res_wd, "m7_r2prior.RDS")), compress=T)

res_m7_ori <- readRDS(file=here(paste0(res_wd, "m7_postpred.RDS")))
res_m7_ori$summary("p_avg1")
res_m7_r2prior$summary("p_avg1")

# prior predictive
data_list$prior_only = 1
res_priorpred_m7_r2prior <- mod_m7_r2prior$sample(
  data = data_list_v2, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_priorpred_m7_r2prior$save_object(file=here(paste0(res_wd, "m7_r2prior_priorpred.RDS")), compress=T)

