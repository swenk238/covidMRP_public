## generating covariates with low prevalence 
# with estimation of sensitivity and specificity
# adding fixed effects to the models 

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
  iter = as.numeric(slurm_arrayid)
  print(iter)
  func_wd <- paste0(wd, 'func/')
  popn_wd <- paste0(wd, "data/") 
  stan_wd <- paste0(wd, 'simdata/sim3_fixedeffects/stan/')
}else{
  #For testing
  iter <- 1
  func_wd <- 'code/simdata/func/'
  popn_wd <- "code/simdata/sim2_spec/data/" 
  stan_wd <- 'code/simdata/sim3_fixedeffects/stan/'
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
source(paste0(func_wd,'functions.R'))

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


# fixed seed number
source(paste0(func_wd,'seed.R'))

# empty lists
samp_data_ss400_list <- lapply(1:length(spec_dup), function(x)lapply(1:length(n_level),function(x)(list())))
samp_data_ss4000_list <- lapply(1:length(spec_dup), function(x)lapply(1:length(n_level),function(x)(list())))
res_mat <- res_sae_mat <- tibble()

# setting seed using array ID
set.seed(seed[iter])

popn_data <- readRDS(paste0(popn_wd, "popn_data_TNFP_stableN.rds"))

# generate a binary variable using x4 to mimic sex variable
popn_data$x4_bin <- ifelse(popn_data$x4_cont < 0, 0, 1)



for(i in 1:length(samp_size)){ 
  # sample different sample sizes from population 
  samp_loc <- sample(popn_size, samp_size[i], replace=F)
  
  samp_data_ori <- popn_data[samp_loc,] 
  
  for(k in 1:length(spec_dup)){ # all TN values
    for(j in 1:length(n_level)){ 
      
      # saving different configurations in the sample data
      samp_data <- samp_data_ori %>% 
        mutate(samp_size = samp_size[i],
               spec = format(spec_dup[k],nsmall=spec_digits ), 
               iter = iter) %>%
        rename_with(~gsub(paste0("_cont_", n_level[j],"levels"), "", .x)) %>%   # renaming the particular level to fit the model
        rename_with(~gsub(paste0("_spec", format(spec_dup[k], nsmall=spec_digits)), "", .x)) %>% # renaming the particular specificity value to fit the model
        mutate_at(c('x1','x2','x3','x4','x5'), as_factor) # convert to factor
      
      dat <- list(x1 = samp_data$x1, x2 = samp_data$x2, x3 = samp_data$x3,
                  x4 = samp_data$x4, x5 = samp_data$x5, y = samp_data$y, 
                  x3_con = as.numeric(samp_data$x3), x4_bin = samp_data$x4_bin, 
                  n = samp_size[i], n_x1 = n_level[j], n_x2 = n_level[j],
                  n_x3 = n_level[j], n_x4 = n_level[j], n_x5 = n_level[j], 
                  coef_prior_scale = 0.5, prior_only = 0,
                  sens=1, spec = round(spec_dup[k],spec_digits),
                  TN = TN_vec[k], FP = FP_vec[k], sd_x3_con = sd(samp_data$x3_cont),
                  K = 3, R2D2_mean_R2 = 0.5, R2D2_prec_R2 = 2,
                                Kscales = 5, R2D2_cons_D2= rep(1, 5))
      
      # model 3 with one fixed effect -----------------------------------------------------------------
      file_m3_FE_R2D2M2 <- file.path(here(paste0(stan_wd,"m3_withFE_R2D2M2.stan")))
      mod_m3_FE_R2D2M2 <- cmdstan_model(file_m3_FE_R2D2M2)
      
      res_m3_FE_R2D2M2 <- mod_m3_FE_R2D2M2$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 4 with one fixed effect -----------------------------------------------------------------
      file_m4_FE_R2D2M2 <- file.path(here(paste0(stan_wd,"m4_withFE_R2D2M2.stan")))
      mod_m4_FE_R2D2M2 <- cmdstan_model(file_m4_FE_R2D2M2)
      
      res_m4_FE_R2D2M2 <- mod_m4_FE_R2D2M2$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 5 with one fixed effect -----------------------------------------------------------------
      file_m5_FE_R2D2M2 <- file.path(here(paste0(stan_wd,"m5_withFE_R2D2M2.stan")))
      mod_m5_FE_R2D2M2 <- cmdstan_model(file_m5_FE_R2D2M2)
      
      res_m5_FE_R2D2M2 <- mod_m5_FE_R2D2M2$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 3 with two fixed effects -----------------------------------------------------------------
      file_m3_twoFE_R2D2M2 <- file.path(here(paste0(stan_wd,"m3_withtwoFE_R2D2M2.stan")))
      mod_m3_twoFE_R2D2M2 <- cmdstan_model(file_m3_twoFE_R2D2M2)
      
      res_m3_twoFE_R2D2M2 <- mod_m3_twoFE_R2D2M2$sample(
        data = dat,
        seed = 1357,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 4 with two fixed effects -----------------------------------------------------------------
      file_m4_twoFE_R2D2M2 <- file.path(here(paste0(stan_wd,"m4_withtwoFE_R2D2M2.stan")))
      mod_m4_twoFE_R2D2M2 <- cmdstan_model(file_m4_twoFE_R2D2M2)
      
      res_m4_twoFE_R2D2M2 <- mod_m4_twoFE_R2D2M2$sample(
        data = dat,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      # model 5 with two fixed effects -----------------------------------------------------------------
      file_m5_twoFE_R2D2M2 <- file.path(here(paste0(stan_wd,"m5_withtwoFE_R2D2M2.stan")))
      mod_m5_twoFE_R2D2M2 <- cmdstan_model(file_m5_twoFE_R2D2M2)
      
      res_m5_twoFE_R2D2M2 <- mod_m5_twoFE_R2D2M2$sample(
        data = dat,
        seed = 2468,
        chains = 4,
        parallel_chains = 4,
        refresh = 500, # print update every 500 iters
      )
      
      pred <- lapply(list(res_m3_FE_R2D2M2, res_m4_FE_R2D2M2, res_m5_FE_R2D2M2,
                          res_m3_twoFE_R2D2M2, res_m4_twoFE_R2D2M2, res_m5_twoFE_R2D2M2), res_stan_summ_func2, sample_data = samp_data, 
                     iteration = iter, nl = n_level[j], sample_size = samp_size[i], 
                     spec = format(spec_dup[k], nsmall=spec_digits), TN = TN_vec[k], FP = FP_vec[k])
      
      res_tb <- do.call(rbind, pred) %>% 
        as_tibble() %>% 
        mutate(model = c("3_oneFE", "4_oneFE", "5_oneFE",
                         "3_twoFE", "4_twoFE", "5_twoFE"))
      
      pred_spec <- lapply(list(res_m3_FE_R2D2M2, res_m4_FE_R2D2M2, res_m5_FE_R2D2M2,
                               res_m3_twoFE_R2D2M2, res_m4_twoFE_R2D2M2, res_m5_twoFE_R2D2M2), 
                          function(x)x$summary('spec')[,c('mean', 'q5', 'q95')] %>% 
                            as_tibble() %>% 
                            rename(spec_mean_q5 = `q5`,
                                   spec_mean = `mean`,
                                   spec_mean_q95 = `q95`) )
      
      res_spec_tb <- do.call(rbind, pred_spec) %>% 
        as_tibble() %>% 
        mutate(model = c("3_oneFE", "4_oneFE", "5_oneFE",
                         "3_twoFE", "4_twoFE", "5_twoFE"))
             
      res_join_tb <- left_join(res_tb, res_spec_tb, by = "model")
      
      res_mat <- rbind(res_mat, res_join_tb) # appending results 
               
      # saving small area estimates 
      pred_sae <- list(res_stan_sae_noncent_func(res_m3_FE_R2D2M2, nl = n_level[j], sample_size = samp_size[i], spec = format(spec_dup[k], nsmall=spec_digits), TN = TN_vec[k], FP = FP_vec[k]),
                       res_stan_sae_noncent_func(res_m4_FE_R2D2M2, nl = n_level[j], sample_size = samp_size[i], spec = format(spec_dup[k], nsmall=spec_digits), TN = TN_vec[k], FP = FP_vec[k]),
                       res_stan_sae_noncent_func(res_m5_FE_R2D2M2, nl = n_level[j], sample_size = samp_size[i], spec = format(spec_dup[k], nsmall=spec_digits), TN = TN_vec[k], FP = FP_vec[k]),
                       res_stan_sae_noncent_func(res_m3_twoFE_R2D2M2, nl = n_level[j], sample_size = samp_size[i], spec = format(spec_dup[k], nsmall=spec_digits), TN = TN_vec[k], FP = FP_vec[k]),
                       res_stan_sae_noncent_func(res_m4_twoFE_R2D2M2, nl = n_level[j], sample_size = samp_size[i], spec = format(spec_dup[k], nsmall=spec_digits), TN = TN_vec[k], FP = FP_vec[k]),
                       res_stan_sae_noncent_func(res_m5_twoFE_R2D2M2, nl = n_level[j], sample_size = samp_size[i], spec = format(spec_dup[k], nsmall=spec_digits), TN = TN_vec[k], FP = FP_vec[k])) 
      
      res_sae_tb <- do.call(rbind, pred_sae) %>% 
        as_tibble() 
      res_sae_mat = rbind(res_sae_mat, res_sae_tb) # appending results 
      
      # saving for different sample sizes and different configurations
      if(i == 1){
        samp_data_ss400_list[[k]][[j]] <- samp_data
      }
      if(i == 2){
        samp_data_ss4000_list[[k]][[j]] <- samp_data
      }
    }
  }
}

saveRDS(res_mat, file=paste0("res_mat_R2D2M2_ite", iter,".RDS"))
saveRDS(res_sae_mat, file=paste0("res_sae_mat_R2D2M2_ite", iter,".RDS"))
