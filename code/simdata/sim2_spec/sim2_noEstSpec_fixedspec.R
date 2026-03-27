## generating covariates with low prevalence 
# without estimating of sensitivity and specificity

library(ggplot2)
library(tidyverse)
library(brms)
library(here)
library(cmdstanr)

options(mc.cores=4)

on_cluster <- F
if(on_cluster){
  wd <- '/hpcfs/users/a1233233/covidMRP/' 
  setwd(wd)
  slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  iter <- as.numeric(slurm_arrayid)
  print(iter)
  func_wd <- paste0(wd, 'func/')
  popn_wd <- paste0(wd, 'data/')
  stan_wd <- paste0(wd, 'simdata/sim2_noEstSpec/stan/')
}else{
  #For testing
  iter <- 1
  func_wd <- 'code/simdata/func/'
  popn_wd <- "code/simdata/sim2_spec/data/" 
  stan_wd <- 'code/simdata/sim2_spec/stan/noEstSpec/'
}

p_set <- 0.01
samp_size <- c(400, 4000)
n_level <- c(4, 10, 20, 40)
beta1_set <- 0.3

popn_size <- 500000

# specifying sensitivity and specificity
sens <- 1
total_vec <- rep(c(400,800,1200,8000),4)
FP_vec <- c(8,16,24,160, 4,8,12,80, 2,4,6,40, rep(0,4))
TN_vec <- total_vec - FP_vec # corresponding to specificity values of 0.98, 0.99, 0.995 and 1 

spec_digits <- 3
spec_dup <- round(TN_vec / (TN_vec + FP_vec), spec_digits)  # creates duplicated spec values
spec <- unique(spec_dup) # only unique spec values

# loading functions to store results
source(here(paste0(func_wd,'functions.R')))

# empty lists
res_mat <- res_sae_mat <- tibble()

# fixed seed number
source(here(paste0(func_wd,'seed.R')))

# setting seed using array ID
set.seed(seed[iter])

popn_data <- readRDS(here(paste0(popn_wd, "popn_data_TNFP_stableN.rds")))

for(i in 1:length(samp_size)){ 
  # sample different sample sizes from population 
  samp_loc <- sample(popn_size, samp_size[i], replace=F)
  
  samp_data_ori <- popn_data[samp_loc,] 
  
  for(k in 1:length(spec)){ # spec values only
    for(j in 1:length(n_level)){ 
      
      # saving different configurations in the sample data
      samp_data <- samp_data_ori %>% 
        mutate(samp_size = samp_size[i],
               spec = format(spec[k],nsmall=spec_digits ), 
               iter = iter) %>%
        rename_with(~gsub(paste0("_cont_", n_level[j],"levels"), "", .x)) %>%   # renaming the particular level to fit the model
        rename_with(~gsub(paste0("_spec", format(spec[k], nsmall=spec_digits)), "", .x)) 
      
      dat <- list(x1 = samp_data$x1, x2 = samp_data$x2, x3 = samp_data$x3,
                  x4 = samp_data$x4, x5 = samp_data$x5, y = samp_data$y, 
                  n = samp_size[i], n_x1 = n_level[j], n_x2 = n_level[j],
                  n_x3 = n_level[j], n_x4 = n_level[j], n_x5 = n_level[j], 
                  coef_prior_scale = 0.5, prior_only = 0,
                  sens=1, spec = spec[k])
      
      # Bayesian mixed effects model - stan code
      # model 0 -----------------------------------------------------------------
      file_m0 <- file.path(here(paste0(stan_wd,"m0_noEstSpec.stan")))
      mod_m0 <- cmdstan_model(file_m0)
      
      res_m0 <- mod_m0$sample(
        data = dat,
        seed = 2468,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 1 -----------------------------------------------------------------
      file_m1 <- file.path(here(paste0(stan_wd, "m1_noEstSpec.stan")))
      mod_m1 <- cmdstan_model(file_m1)
      
      res_m1 <- mod_m1$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 2 -----------------------------------------------------------------
      file_m2 <- file.path(here(paste0(stan_wd,"m2_noEstSpec.stan")))
      mod_m2 <- cmdstan_model(file_m2)
      
      res_m2 <- mod_m2$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 3 -----------------------------------------------------------------
      file_m3 <- file.path(here(paste0(stan_wd,"m3_noEstSpec.stan")))
      mod_m3 <- cmdstan_model(file_m3)
      
      res_m3 <- mod_m3$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 4 -----------------------------------------------------------------
      file_m4 <- file.path(here(paste0(stan_wd,"m4_noEstSpec.stan")))
      mod_m4 <- cmdstan_model(file_m4)
      
      res_m4 <- mod_m4$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 5 -----------------------------------------------------------------
      file_m5 <- file.path(here(paste0(stan_wd,"m5_noEstSpec.stan")))
      mod_m5 <- cmdstan_model(file_m5)
      
      res_m5 <- mod_m5$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      pred <- lapply(list(res_m0, res_m1, res_m2, 
                          res_m3, res_m4, res_m5), res_stan_summ_func, sample_data = samp_data, 
                     iteration = iter, nl = n_level[j], sample_size = samp_size[i], 
                     spec = spec[k], TN = NA, FP = NA)
      
      res_tb <- do.call(rbind, pred) %>% 
        as_tibble() %>% 
        mutate(model = c(0:5))
      
      res_mat <- rbind(res_mat, res_tb) # appending results 
      
      # saving small area estimates 
      pred_sae <- list(res_stan_sae_noncent_func(res_m0, nl = n_level[j], sample_size = samp_size[i], spec = spec[k], TN = NA, FP = NA),
                       res_stan_sae_noncent_func(res_m1, nl = n_level[j], sample_size = samp_size[i], spec = spec[k], TN = NA, FP = NA),
                       res_stan_sae_noncent_func(res_m2, nl = n_level[j], sample_size = samp_size[i], spec = spec[k], TN = NA, FP = NA),
                       res_stan_sae_noncent_func(res_m3, nl = n_level[j], sample_size = samp_size[i], spec = spec[k], TN = NA, FP = NA),
                       res_stan_sae_noncent_func(res_m4, nl = n_level[j], sample_size = samp_size[i], spec = spec[k], TN = NA, FP = NA),
                       res_stan_sae_noncent_func(res_m5, nl = n_level[j], sample_size = samp_size[i], spec = spec[k], TN = NA, FP = NA)) 
      
      res_sae_tb <- do.call(rbind, pred_sae) %>% 
        as_tibble() 
      res_sae_mat = rbind(res_sae_mat, res_sae_tb) # appending results 
    }
  }
}

res_wd <- 'simdata/sim2_noEstSpec/'

saveRDS(res_mat, file=paste0(res_wd, "res_mat_ite", iter,".RDS"))
saveRDS(res_sae_mat, file=paste0(res_wd, "res_sae_mat_ite", iter,".RDS"))

