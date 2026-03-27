## setting stronger priors on flipped y's 
## 17/05/2023

library(here)
source(here("code/data_list.R"), echo=F)
options(mc.cores = 8)

# flipping y's so that 1's are 0's and 0's are 1's - 
data_list$y_flip = data_list$y+1
data_list$y_flip = ifelse(data_list$y_flip == "2", "0", "1")
table(data_list$y_flip)

data_list$y = as.integer(data_list$y_flip)
table(data_list$y)
data_list["y_flip"] = NULL

# Multiplying beta values x10 for specificity only:
(data_list$TN = data_list$TN*100 )
(data_list$FP = data_list$FP*100 )

# model 5 - adapting Marnie's code to do prior predictive checks --------------------
file_m5 <- file.path(here("code/stan/model5.stan"))
mod_m5 <- cmdstan_model(file_m5)

## sampling ####
# posterior predictive checks 

res_posterior_m5 <- mod_m5$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

res_posterior_m5$summary('p_avg1')

# model 4 - without agegp -------------------------------------------------
file_m4 <- file.path(here("code/stan/model4.stan"))
mod_m4 <- cmdstan_model(file_m4)

## sampling ####
res_posterior_m4 <- mod_m4$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 3 without agegp, sex ----------------------------------------------
file_m3 <- file.path(here("code/stan/model3.stan"))
mod_m3 <- cmdstan_model(file_m3)

## sampling ####
res_posterior_m3 <- mod_m3$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 2 without agegp, sex, seifa ----------------------------------------------
file_m2 <- file.path(here("code/stan/model2.stan"))
mod_m2 <- cmdstan_model(file_m2)

## sampling ####
res_posterior_m2 <- mod_m2$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 1 without agegp, sex, seifa, postcode ----------------------------------------------
file_m1 <- file.path(here("code/stan/model1.stan"))
mod_m1 <- cmdstan_model(file_m1)

## sampling ####
res_posterior_m1 <- mod_m1$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)


save(res_posterior_m1, res_posterior_m2,
     res_posterior_m3, res_posterior_m4, 
     res_posterior_m5, file=here("results/experiments/exp3/strongprior.RData"))

