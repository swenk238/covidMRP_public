# multiplying counts x10 in original data to see if same problem exists
# May 2023 
library(here)
library(cmdstanr)
source(here("code/data_list.R"), echo=F)
options(mc.cores = 8)

# randomly sample individuals to have count = 1, incl those who already are 1
set.seed(1234)
ori = which(data_list$y==1)
ind = 1:4799

# fit all models for modified counts
names = c("5percenty", "16percenty", "95percenty", "99percenty")
frac = c(0.05, 0.16, 0.95, 0.99)
ind = 1:4780

for (i in 1:4){
  set.seed(123)
  source(here("code/data_list.R"), echo=F)
  ori = which(data_list$y==1) # indicators of original positives 
  
  newcount = round(frac[i]*4780)
  
  samp1 = sample(ind[-ori], newcount - 77, replace=F) # sample newcount-77 inv to have count=1
  
  data_list$y[samp1] = 1
  sum(data_list$y)
  
  # model 5 - adapting Marnie's code to do prior predictive checks --------------------
  file_m5 <- file.path(here("code/stan/model5.stan"))
  mod_m5 <- cmdstan_model(file_m5)
  
  # sampling
  res_posterior_m5 <- mod_m5$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  res_posterior_m5$summary("p_avg1")
  
  # saving
  res_posterior_m5$save_object(file=here(paste0("results/experiments/exp6/m5_",names[i], ".RDS")))
  
  
  # model 4 - --------------------
  file_m4 <- file.path(here("code/stan/model4.stan"))
  mod_m4 <- cmdstan_model(file_m4)
  
  # sampling
  res_posterior_m4 <- mod_m4$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  res_posterior_m4$summary("p_avg1")
  
  # saving
  res_posterior_m4$save_object(file=here(paste0("results/experiments/exp6/m4_",names[i], ".RDS")))
  
  
  # model 3 ---------------------
  file_m3 <- file.path(here("code/stan/model3.stan"))
  mod_m3 <- cmdstan_model(file_m3)
  
  # sampling
  res_posterior_m3 <- mod_m3$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  res_posterior_m3$summary("p_avg1")
  
  # saving
  res_posterior_m3$save_object(file=here(paste0("results/experiments/exp6/m3_",names[i], ".RDS")))
  
  
  # model 2 ---------------------
  file_m2 <- file.path(here("code/stan/model2.stan"))
  mod_m2 <- cmdstan_model(file_m2)
  
  # sampling
  res_posterior_m2 <- mod_m2$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  res_posterior_m2$summary("p_avg1")
  
  # saving
  res_posterior_m2$save_object(file=here(paste0("results/experiments/exp6/m2_",names[i], ".RDS")))
  
  
  # model 1 ---------------------
  file_m1 <- file.path(here("code/stan/model1.stan"))
  mod_m1 <- cmdstan_model(file_m1)
  
  # sampling
  res_posterior_m1 <- mod_m1$sample(
    data = data_list, 
    seed = 2345, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 500, # print update every 500 iters
  )
  res_posterior_m1$summary("p_avg1")
  
  # saving
  res_posterior_m1$save_object(file=here(paste0("results/experiments/exp6/m1_",names[i], ".RDS")))
  
}

