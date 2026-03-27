## generating covariates with low prevalence 
## only investigating when p(y) = 0.01

library(ggplot2)
library(arm)
library(tidyverse)
library(brms)
library(testthat)
library(here)
options(mc.cores=4)

# iteration number (ITE) from cluster
on_cluster <- F
if(on_cluster){
  wd <- '.' 
  setwd(wd)
  func_wd <- '../func/'
  res_wd <- './'
}else{
  #For testing
  iter <- 1
  func_wd <- 'code/simdata/func/'
  res_wd <- 'results/simdata/sim1/covariates/cluster/'
}

samp_size = c(400, 4000)
n_level = c(4, 10, 20, 40)
beta1_set = c(0.3, 0)
popn_size = 500000

a = -0.5 # Uniform distribution
b = 0.5 
p_set = 0.01
cov_num = 5

# solving for beta0, should be ~ logit(p_set) since E(X) = 0
(b0 = logit(p_set) - (0.3*(a+b)/2)*cov_num) 

# loading functions to store results
source(paste0(here(func_wd,'/functions.R')))

# empty lists and tibbles
popn_data_list = lapply(1:length(beta1_set), function(x)matrix(NA,  nrow=popn_size))
samp_data_ss400_list = lapply(1:length(beta1_set), function(x)lapply(1:length(n_level),function(x)(list())))
samp_data_ss4000_list = lapply(1:length(beta1_set), function(x)lapply(1:length(n_level),function(x)(list())))
res_mat = res_sae_mat = tibble()

# generate a finite population
set.seed(2468)
for(k in 1:length(beta1_set)){ # number of matrices
  x_cont_mat = matrix(NA, ncol=5, nrow=popn_size)
  
  # generate uniform x's
  for (ind in 1:5){
    x_cont_mat[,ind] = runif(popn_size,min=a,max=b)
  }
  
  p_y = invlogit(b0 + as.matrix(x_cont_mat[,1:5]) %*% rep(beta1_set[k],cov_num))
  y = rbinom(popn_size,1, p_y)
  
  popn_data = cbind(x_cont_mat, p_y, y) %>%
    as_tibble() 
  colnames(popn_data) = c('x1_cont', 'x2_cont', 'x3_cont', 'x4_cont', 'x5_cont', 'p_y', "y") 
  
  # saving in list 
  # discretizing continuous x for different number of levels and renaming
  popn_data_list[[k]] = popn_data %>% 
    mutate(across(x1_cont:x5_cont, 
                  ~cut_number(.x, n = 4, labels=F), 
                  .names= "{col}_4levels")) %>%
    mutate(across(x1_cont:x5_cont, 
                  ~cut_number(.x, n = 10, labels=F), 
                  .names= "{col}_10levels")) %>%
    mutate(across(x1_cont:x5_cont, 
                  ~cut_number(.x, n = 20, labels=F), 
                  .names= "{col}_20levels")) %>%
    mutate(across(x1_cont:x5_cont, 
                  ~cut_number(.x, n = 40, labels=F), 
                  .names= "{col}_40levels")) 
}

# fixed seed number
source(here(paste0(func_wd,'seed.R')))

missing_rows <- structure(list(n_level = c(4, 10, 20, 40, 4, 10, 20, 40, 4, 10, 
                                           20, 40, 4, 10, 20, 40, 4, 10, 20, 40),
                               beta1_set = c(0.3, 0.3, 0.3, 0.3, 0, 0, 0, 0, 0, 0,
                                             0, 0, 0.3, 0.3, 0.3, 0.3, 0, 0, 0, 0), 
                               samp_size = c(400, 400, 400, 400, 400, 400, 400, 400, 400, 400,
                                             400, 400, 400, 400, 400, 400, 400, 400, 400, 400), 
                               n = c(6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L,
                                     6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L), 
                               iter = c(1L, 1L, 1L, 1L, 9L, 9L, 9L, 9L, 56L, 56L, 56L,
                                        56L, 63L, 63L, 63L, 63L, 66L, 66L, 66L, 66L)), 
                          row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"))
                                                                                                                                                                                                                                                                  
# taking samples and fitting different models to them according to different specifications
for(index in 1:nrow(missing_rows)){
  
  iter <- missing_rows[index, 'iter'] %>% pull()
  print(iter)
  
  # setting seed using array ID
  set.seed(seed[iter]-1)
  
  
  k <- which(beta1_set == missing_rows[index,]$beta1_set)
  i <- which(samp_size == missing_rows[index,]$samp_size)
  j <- which(n_level == missing_rows[index,]$n_level)
  
  # population data 
  popn_data <- popn_data_list[[k]]
  
  # sample different sample sizes from population 
  samp_loc <- sample(popn_size, samp_size[i], replace=F)
  samp_data_ori <- popn_data[samp_loc,] 
  
  # saving different configurations in the sample data
  samp_data <- samp_data_ori %>% 
    mutate(samp_size = samp_size[i],
           n_level = n_level[j],
           beta1_set = beta1_set[k], 
           iter = iter) %>% 
    rename_with(~gsub(paste0("_cont_", n_level[j],"levels"), "", .x))   # renaming the particular level to fit the model
  
  
  ## Fitting the models ####
  # frequentist model
  fit0.lmer = glm(y ~ 1, family=binomial, data=samp_data)
  
  fit1.lmer = glmer(y ~ (1|x1), family=binomial, data=samp_data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  
  fit2.lmer = glmer(y ~ (1|x1) + (1|x2), family=binomial, data=samp_data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  
  fit3.lmer = glmer(y ~ (1|x1) + (1|x2) + (1|x3), family=binomial, data=samp_data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  
  fit4.lmer = glmer(y ~ (1|x1) + (1|x2) + (1|x3) + (1|x4), family=binomial, data=samp_data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  
  fit5.lmer = glmer(y ~ (1|x1) + (1|x2) + (1|x3) + (1|x4) + (1|x5), family=binomial, data=samp_data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  
  # Bayesian mixed effects model
  fit0.brm = brm(y ~ 1, data=samp_data, family = bernoulli, backend="cmdstanr", adapt_delta = 0.99)
  fit1.brm = brm(y ~ (1|x1), data=samp_data, family = bernoulli, backend="cmdstanr", adapt_delta = 0.99)
  fit2.brm = brm(y ~ (1|x1) + (1|x2), data=samp_data, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.99)
  fit3.brm = brm(y ~ (1|x1) + (1|x2) + (1|x3), data=samp_data, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.99)
  fit4.brm = brm(y ~ (1|x1) + (1|x2) + (1|x3) + (1|x4), data=samp_data, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.99)
  fit5.brm = brm(y ~ (1|x1) + (1|x2) + (1|x3) + (1|x4) + (1|x5), data=samp_data, family = bernoulli, backend = "cmdstanr", adapt_delta = 0.99)
  
  # getting the predictions, then posterior distribution of the mean, then take summaries of that 
  pred = lapply(list(fit0.brm, fit1.brm, 
                     fit2.brm, fit3.brm, 
                     fit4.brm, fit5.brm), multiple_func, 
                sample_data = samp_data, iteration = iter, nl = n_level[j], sample_size = samp_size[i], b1 = beta1_set[k]) 
  
  res_tb = do.call(rbind, pred) %>% 
    as_tibble() %>% 
    mutate(model = c(0:5),
           glm0_pred_mean =  predict(fit0.lmer, type="response") %>% mean(),
           glm1_pred_mean =  predict(fit1.lmer, type="response") %>% mean(),
           glm2_pred_mean =  predict(fit2.lmer, type="response") %>% mean(),
           glm3_pred_mean =  predict(fit3.lmer, type="response") %>% mean(),
           glm4_pred_mean =  predict(fit4.lmer, type="response") %>% mean(),
           glm5_pred_mean =  predict(fit5.lmer, type="response") %>% mean())
  res_mat = rbind(res_mat, res_tb) # appending results 
  
  # saving small area estimates 
  pred_sae = list(multiple_func_sae(fit0.brm, nl = n_level[j], iteration=iter, sample_size = samp_size[i], b1 = beta1_set[k]),
                  multiple_func_sae(fit1.brm, nl = n_level[j], iteration=iter, sample_size = samp_size[i], b1 = beta1_set[k]),
                  multiple_func_sae(fit2.brm, nl = n_level[j], iteration=iter, sample_size = samp_size[i], b1 = beta1_set[k]),
                  multiple_func_sae(fit3.brm, nl = n_level[j], iteration=iter, sample_size = samp_size[i], b1 = beta1_set[k]),
                  multiple_func_sae(fit4.brm, nl = n_level[j], iteration=iter, sample_size = samp_size[i], b1 = beta1_set[k]),
                  multiple_func_sae(fit5.brm, nl = n_level[j], iteration=iter, sample_size = samp_size[i], b1 = beta1_set[k])) # can't use lapply because of extracting model name
  res_sae_tb = do.call(rbind, pred_sae) %>% 
    as_tibble() 
  res_sae_mat = rbind(res_sae_mat, res_sae_tb) # appending results 
}

saveRDS(res_mat, file=paste0(res_wd, "res_mat_extra_all.RDS"))
saveRDS(res_sae_mat, file=paste0(res_wd, "res_sae_mat_extra_all.RDS"))

