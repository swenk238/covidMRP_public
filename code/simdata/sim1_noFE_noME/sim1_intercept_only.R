## testing intercept-only model 
library(brms)
library(tidyverse)
options(mc.cores=4)

# iteration number (ITE) from cluster
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
iter = as.numeric(slurm_arrayid)

p_config = c(0.001 ,0.01, 0.1, 0.2)
sample_size = c(20, 40, 400, 4000, 40000)
popn_size = 500000

# storing the results of mean_y, glm and brms
res_mat = tibble() # should have total rows of p_config X sample size (20)
y_popn = lapply(1:length(p_config), 
                function(x)matrix(NA, 
                                  nrow=popn_size))

# generate a finite population for different p_config specifications
set.seed(2468)
for(i in 1:length(p_config)){
  y_popn[[i]] = rbinom(popn_size, 1, p_config[i])
}

# fixed seed number
pn = 100 # number of different iterations
seed = round(runif(pn, min=10, max=100000),0)

# setting seed using array ID
set.seed(seed[iter])

for(i in 1:length(p_config)){ # for all the different p(y) values
  for(j in 1:length(sample_size)){ # for all the different sample sizes
    # random sampling according to different sample sizes 
    samp_loc = sample(popn_size,sample_size[j], replace=F)
    y = y_popn[[i]][samp_loc]

    # glm model
    fit_glm = glm(y~1, family ="binomial")

    # brm estimates 
    fit_brm = brm(y ~ 1, data=list(y=y), family = bernoulli, backend = "cmdstanr", adapt_delta = 0.98)
    
    # prior_only model
    fit_brm_prior = brm(y ~ 1, data=list(y=y), family = bernoulli, backend = "cmdstanr", adapt_delta = 0.98, sample_prior = "only")
    res_tb_prior = posterior_linpred(fit_brm_prior, transform=T) %>% 
      rowMeans() %>% # taking mean for each posterior sample 
      quantile(., c(0.05,0.5,0.95)) %>% # taking quantiles of the posterior distribution of the mean 
      t() %>%
      as_tibble() %>% 
      rename(brm_prior_pred_mean_q5 = `5%`,
             brm_prior_pred_mean_q50 = `50%`,
             brm_prior_pred_mean_q95 = `95%`) 
    
    # saving results
    res_tb_row = posterior_linpred(fit_brm, transform=T) %>% 
      rowMeans() %>% # taking mean for each posterior sample 
      quantile(., c(0.05,0.5,0.95)) %>% # taking quantiles of the posterior distribution of the mean 
      t() %>%    
      as_tibble() %>% 
      rename(brm_pred_mean_q5 = `5%`,
             brm_pred_mean_q50 = `50%`,
             brm_pred_mean_q95 = `95%`) %>%
      mutate(glm_pred_mean = predict(fit_glm, type="response") %>% mean(),
             mean_y = mean(y),
             coverage_brm = ifelse(brm_pred_mean_q5 <= mean_y & brm_pred_mean_q95 >= mean_y, 1, 0),
             p_config = p_config[i],
             samp_size = sample_size[j],
             rhat_larger_1_1 = sum(rhat(fit_brm)>1.1), # diagnostics
             n_divergent = sum(subset(nuts_params(fit_brm), Parameter == "divergent__")$Value),  
             n_maxtreedepth = sum(subset(nuts_params(fit_brm), Parameter == "treedepth__")$Value ==10),
             ess_bulk =  posterior::ess_bulk(fit_brm),
             post_intercept = summary(fit_brm)$fixed[1] %>% as.numeric() %>% inv_logit_scaled(),
             prior_intercept = summary(fit_brm_prior)$fixed[1] %>% as.numeric() %>% inv_logit_scaled()) %>% 
            cbind(res_tb_prior)
    
    res_mat = rbind(res_mat, res_tb_row) # appending results 
  }
}

saveRDS(res_mat, file=paste0("res_mat_ite", iter,".RDS"))
saveRDS(y_popn, file=paste0("y_popn_ite", iter, ".RDS"))