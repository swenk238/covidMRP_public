## generating covariates with low prevalence 
# without estimating of sensitivity and specificity

library(ggplot2)
library(arm) # invlogit()
library(tidyverse)
library(brms)
library(here)
library(cmdstanr)

options(mc.cores=4)

on_cluster <- T
if(on_cluster){
  wd <- '~/Mona0085/skuh/covidMRP/simdata/exp4_noEstSpec_stableN' 
  setwd(wd)
  slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  iter <- as.numeric(slurm_arrayid)
  print(iter)
  func_wd <- '../func/'
  popn_wd <- "../exp4/data/" 
}else{
  #For testing
  iter <- 1
  func_wd <- 'code/simdata/func/'
  popn_wd <- "code/simdata/exp4/data/" 
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
samp_data_ss400_list <- lapply(1:length(spec_dup), function(x)lapply(1:length(n_level),function(x)(list())))
samp_data_ss4000_list <- lapply(1:length(spec_dup), function(x)lapply(1:length(n_level),function(x)(list())))
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
  
  for(k in 1:length(spec_dup)){ # all TN values
    for(j in 1:length(n_level)){ 
      
      # saving different configurations in the sample data
      samp_data <- samp_data_ori %>% 
        mutate(samp_size = samp_size[i],
               spec = format(spec_dup[k],nsmall=spec_digits ), 
               iter = iter) %>%
        rename_with(~gsub(paste0("_cont_", n_level[j],"levels"), "", .x)) %>%   # renaming the particular level to fit the model
        rename_with(~gsub(paste0("_spec", format(spec_dup[k], nsmall=spec_digits)), "", .x)) 
     
       # Bayesian mixed effects model
      fit0.brm <- brm(y ~ 0 + Intercept, data=samp_data, family = bernoulli, backend="cmdstanr", adapt_delta = 0.995)
      fit1.brm <- brm(y ~ 0 + Intercept + (1|x1), data=samp_data, family = bernoulli, backend="cmdstanr", adapt_delta = 0.995)
      fit2.brm <- brm(y ~ 0 + Intercept + (1|x1) + (1|x2), data=samp_data, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.995)
      fit3.brm <- brm(y ~ 0 + Intercept + (1|x1) + (1|x2) + (1|x3), data=samp_data, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.995)
      fit4.brm <- brm(y ~ 0 + Intercept + (1|x1) + (1|x2) + (1|x3) + (1|x4), data=samp_data, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.995)
      fit5.brm <- brm(y ~ 0 + Intercept + (1|x1) + (1|x2) + (1|x3) + (1|x4) + (1|x5), data=samp_data, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.995)
      
      # getting the predictions, then posterior distribution of the mean, then take summaries of that 
      pred_summ <- lapply(list(fit0.brm, fit1.brm, 
                               fit2.brm, fit3.brm, 
                               fit4.brm, fit5.brm), multiple_func, 
                          sample_data = samp_data, iteration = iter, nl = n_level[j], sample_size = samp_size[i], b1 = beta1_set) 
      pred <- lapply(pred_summ, function(x)mutate(x, spec = spec_dup[k], TN = TN_vec[k], FP = FP_vec[k]))
      
      res_tb <- do.call(rbind, pred) %>% 
        as_tibble() %>% 
        mutate(model = c(0:5))
      res_mat <- rbind(res_mat, res_tb) # appending results 
      
      # saving small area estimates 
      pred_sae_summ <- list(multiple_func_sae(fit0.brm, nl = n_level[j], sample_size = samp_size[i], b1 = beta1_set, iteration = iter),
                            multiple_func_sae(fit1.brm, nl = n_level[j], sample_size = samp_size[i], b1 = beta1_set, iteration = iter),
                            multiple_func_sae(fit2.brm, nl = n_level[j], sample_size = samp_size[i], b1 = beta1_set, iteration = iter),
                            multiple_func_sae(fit3.brm, nl = n_level[j], sample_size = samp_size[i], b1 = beta1_set, iteration = iter),
                            multiple_func_sae(fit4.brm, nl = n_level[j], sample_size = samp_size[i], b1 = beta1_set, iteration = iter),
                            multiple_func_sae(fit5.brm, nl = n_level[j], sample_size = samp_size[i], b1 = beta1_set, iteration = iter)) # can't use lapply because of extracting model name
      pred_sae <- lapply(pred_sae_summ, function(x)mutate(x, spec = spec_dup[k], TN = TN_vec[k], FP = FP_vec[k]))
      
      res_sae_tb <- do.call(rbind, pred_sae) %>% 
        as_tibble() 
      res_sae_mat <- rbind(res_sae_mat, res_sae_tb) # appending results 
      
      # saving in diff lists for different sample sizes
      if(i == 1){
        samp_data_ss400_list[[k]][[j]] <- samp_data
      }
      if(i == 2){
        samp_data_ss4000_list[[k]][[j]] <- samp_data
      }
    }
  }
}

saveRDS(res_mat, file=paste0("res_mat_ite", iter,".RDS"))
saveRDS(res_sae_mat, file=paste0("res_sae_mat_ite", iter,".RDS"))

saveRDS(samp_data_ss400_list, file=paste0("samp_data_ss400_list_ite", iter,".RDS"))
saveRDS(samp_data_ss4000_list, file=paste0("samp_data_ss4000_list_ite", iter,".RDS"))

