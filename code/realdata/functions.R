# posterior mean from stan code
res_stan_func <- function(modelfit, spec_values, TN_value, FP_value){
  t1 <- modelfit$draws(format="df")  %>%
    dplyr::select(contains(c("p["))) %>%  
    rowMeans() %>% # taking mean for each posterior sample 
    quantile(., c(0.05,0.5,0.95)) %>% 
    t() %>% 
    as_tibble() %>% 
    rename(pred_mean_q5 = `5%`,
           pred_mean_q50 = `50%`,
           pred_mean_q95 = `95%`) 
  t2 <- modelfit$draws(format="df")  %>%
    dplyr::select(contains(c("p_sample["))) %>%  
    rowMeans() %>% # taking mean for each posterior sample 
    quantile(., c(0.05,0.5,0.95)) %>% 
    t() %>% 
    as_tibble() %>% 
    rename(psamp_mean_q5 = `5%`,
           psamp_mean_q50 = `50%`,
           psamp_mean_q95 = `95%`) 
  cbind(t1,t2) %>% 
    mutate(mean_y_samp = mean(sample_data$y),
           coverage = ifelse(pred_mean_q5 <= mean_y_samp & pred_mean_q95 >= mean_y_samp, 1, 0),
           rhat_larger_1_1 = sum(modelfit$summary()[,'rhat']>1.1), # diagnostics
           n_divergent = sum(modelfit$diagnostic_summary()['num_divergent'][[1]]),
           n_maxtreedepth = sum(modelfit$diagnostic_summary()['num_max_treedepth'][[1]]),
           ess_bulk =  modelfit$summary()[,'ess_bulk'] %>% min(),
           est_int = modelfit$summary('b[1]')[2] %>% as.numeric() %>% inv_logit_scaled(),
           spec = spec_values,
           TN = TN_value, 
           FP = FP_value)
}


res_stan_func_seifa <- function(modelfit, spec_values, TN_value, FP_value){
  modelfit$draws(format="df")  %>%
    dplyr::select(contains(c("p["))) %>%  
    rowMeans() %>% # taking mean for each posterior sample 
    quantile(., c(0.05,0.5,0.95)) %>% 
    t() %>% 
    as_tibble() %>% 
    rename(pred_mean_q5 = `5%`,
           pred_mean_q50 = `50%`,
           pred_mean_q95 = `95%`) %>% 
    mutate(mean_y_samp = mean(sample_data$y),
           coverage = ifelse(pred_mean_q5 <= mean_y_samp & pred_mean_q95 >= mean_y_samp, 1, 0),
           rhat_larger_1_1 = sum(modelfit$summary()[,'rhat']>1.1), # diagnostics
           n_divergent = sum(modelfit$diagnostic_summary()['num_divergent'][[1]]),
           n_maxtreedepth = sum(modelfit$diagnostic_summary()['num_max_treedepth'][[1]]),
           ess_bulk =  modelfit$summary()[,'ess_bulk'] %>% min(),
           est_int = modelfit$summary('b[1]')[2] %>% as.numeric() %>% inv_logit_scaled(),
           spec = spec_values,
           TN = TN_value, 
           FP = FP_value)
}


# extracting results directly from model fit
res_stan_summ_func <- function(modelfit, spec_values, TN_value, FP_value){
  t1 <- modelfit$summary('p_avg1')[,c('mean', 'q5', 'q95')] %>% 
    as_tibble() %>% 
    rename(pred_mean_q5 = `q5`,
           pred_mean_q50 = `mean`,
           pred_mean_q95 = `q95`) 
  t2 <- modelfit$summary('p_sample_mean')[,c('mean', 'q5', 'q95')] %>% 
    as_tibble() %>% 
    rename(psamp_mean_q5 = `q5`,
           psamp_mean_q50 = `mean`,
           psamp_mean_q95 = `q95`) 
  cbind(t1,t2) %>% 
    mutate(mean_y_samp = mean(data_list$y),
           coverage = ifelse(pred_mean_q5 <= mean_y_samp & pred_mean_q95 >= mean_y_samp, 1, 0),
           rhat_larger_1_1 = sum(modelfit$summary()[,'rhat']>1.1), # diagnostics
           n_divergent = sum(modelfit$diagnostic_summary()['num_divergent'][[1]]),
           n_maxtreedepth = sum(modelfit$diagnostic_summary()['num_max_treedepth'][[1]]),
           ess_bulk =  modelfit$summary()[,'ess_bulk'] %>% min(),
           est_int = modelfit$summary('b')[2] %>% as.numeric() %>% inv_logit_scaled(),
           spec = spec_values,
           TN = TN_value, 
           FP = FP_value)
}
