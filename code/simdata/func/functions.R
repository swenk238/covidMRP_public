## functions to use in the simulations 

# function to get the quantiles and mean
multiple_func <- function(modelfit, sample_data, iteration, nl, sample_size, b1) {
  posterior_linpred(modelfit, transform=T) %>% 
    rowMeans() %>% # taking mean for each posterior sample 
    quantile(., c(0.05,0.5,0.95)) %>% # taking quantiles of the posterior distribution of the mean 
    t() %>% 
    as_tibble() %>% 
    rename(brm_pred_mean_q5 = `5%`,
           brm_pred_mean_q50 = `50%`,
           brm_pred_mean_q95 = `95%`) %>% 
    mutate(iter = iteration,
           mean_y_samp = mean(sample_data$y),
           coverage_brm = ifelse(brm_pred_mean_q5 <= mean_y_samp & brm_pred_mean_q95 >= mean_y_samp, 1, 0),
           n_level = nl,
           samp_size = sample_size,
           beta1_set = b1,
           rhat_larger_1_1 = sum(rhat(modelfit)>1.1), # diagnostics
           n_divergent = sum(subset(nuts_params(modelfit), Parameter == "divergent__")$Value),  
           n_maxtreedepth = sum(subset(nuts_params(modelfit), Parameter == "treedepth__")$Value ==10),
           ess_bulk =  posterior::ess_bulk(modelfit),
           est_int = median(as_draws_df(modelfit)$b_Intercept) %>% as.numeric() %>% inv_logit_scaled())
}

# function to get sae estimates
multiple_func_sae <- function(modelfit, iteration, nl, sample_size, b1){
  modelname = deparse(substitute(modelfit)) # getting model name
  as_tibble(modelfit) %>% 
    dplyr::select(contains(c("b_","r_"))) %>% 
    apply(2, quantile,c(0.05,0.5,0.95)) %>%
    as_tibble() %>%
    transpose_df() %>% 
    rename(iter = iteration, 
           variable = rowname, 
           brm_pred_mean_q5 = `1`,
           brm_pred_mean_q50 = `2`,
           brm_pred_mean_q95 = `3`) %>% 
    mutate(iter = iter,
           samp_size = sample_size, 
           n_level = nl,
           beta1_set = b1,
           model = gsub(".*([0-9]).*$", "\\1", modelname))
}

# to transpose tibble while preserving name
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}

# posterior mean from stan code
res_stan_func <- function(modelfit, sample_data, iteration, nl, sample_size, spec_values, TN_value, FP_value){
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
           iter = iteration,
           n_level = nl,
           samp_size = sample_size,
           rhat_larger_1_1 = sum(modelfit$summary()[,'rhat']>1.1), # diagnostics
           n_divergent = sum(modelfit$diagnostic_summary()['num_divergent'][[1]]),
           n_maxtreedepth = sum(modelfit$diagnostic_summary()['num_max_treedepth'][[1]]),
           ess_bulk =  modelfit$summary()[,'ess_bulk'] %>% min(),
           est_int = modelfit$summary('b')[2] %>% as.numeric() %>% inv_logit_scaled(),
	   spec = spec_values,
	   TN = TN_value, 
	   FP = FP_value)
}


# to extract results from stan code
res_stan_sae_func <- function(modelfit, nl, sample_size, spec_values,TN_value, FP_value){
  modelname = deparse(substitute(modelfit)) # getting model name
 modelfit$draws(format="df")  %>%
    dplyr::select(contains(c("b", "r_"))) %>% 
    apply(2, function(x)quantile(x,c(0.05,0.5,0.95))) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble() %>%
    rename(q5 = `5%`,
           q50 = `50%`,
           q95 = `95%`,
           variable = rowname) %>%
    mutate(samp_size = sample_size, 
           n_level = nl,
           spec = spec_values,
	   TN = TN_value, 
	   FP = FP_value, 
           model = gsub(".*([0-9]).*$", "\\1", modelname))
}

# posterior mean from stan code (non-centered paramterisation)
# to extract results from stan code (non-parametrisation)
res_stan_sae_noncent_func <- function(modelfit, nl, sample_size, spec_values,TN_value, FP_value){
  modelname = deparse(substitute(modelfit)) # getting model name
  modelfit$draws(format="df")  %>%
    dplyr::select(starts_with(c("b", "a"))) %>% 
    apply(2, function(x)quantile(x,c(0.05,0.5,0.95))) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_tibble() %>%
    rename(q5 = `5%`,
           q50 = `50%`,
           q95 = `95%`,
           variable = rowname) %>%
    mutate(samp_size = sample_size, 
           n_level = nl,
           spec = spec_values,
           TN = TN_value, 
           FP = FP_value, 
           model = gsub(".*([0-9]).*$", "\\1", modelname))
}

res_stan_summ_func <- function(modelfit, sample_data, iteration, nl, sample_size, spec_values, TN_value, FP_value){
  t1 <- modelfit$summary('p_avg1')[,c('mean', 'q5', 'q95')] %>% 
    as_tibble() %>% 
    rename(pred_mean_q5 = `q5`,
           pred_mean_q50 = `mean`,
           pred_mean_q95 = `q95`) 
  t2 <- modelfit$summary('p_samp_avg1')[,c('mean', 'q5', 'q95')] %>% 
    as_tibble() %>% 
    rename(psamp_mean_q5 = `q5`,
           psamp_mean_q50 = `mean`,
           psamp_mean_q95 = `q95`) 
  cbind(t1,t2) %>% 
    mutate(mean_y_samp = mean(sample_data$y),
           coverage = ifelse(pred_mean_q5 <= mean_y_samp & pred_mean_q95 >= mean_y_samp, 1, 0),
           iter = iteration,
           n_level = nl,
           samp_size = sample_size,
           rhat_larger_1_1 = sum(modelfit$summary()[,'rhat']>1.1), # diagnostics
           n_divergent = sum(modelfit$diagnostic_summary()['num_divergent'][[1]]),
           n_maxtreedepth = sum(modelfit$diagnostic_summary()['num_max_treedepth'][[1]]),
           ess_bulk =  modelfit$summary()[,'ess_bulk'] %>% min(),
           est_int = modelfit$summary('b')[2] %>% as.numeric() %>% inv_logit_scaled(),
           spec = spec_values,
           TN = TN_value, 
           FP = FP_value)
}


res_stan_summ_func2 <- function(modelfit, sample_data, iteration, nl, sample_size, spec_values, TN_value, FP_value){
  t1 <- modelfit$summary('p_avg1')[,c('mean', 'q5', 'q95')] %>% 
    as_tibble() %>% 
    rename(pred_mean_q5 = `q5`,
           pred_mean_q50 = `mean`,
           pred_mean_q95 = `q95`) 
  t2 <- modelfit$summary('p_samp_avg1')[,c('mean', 'q5', 'q95')] %>% 
    as_tibble() %>% 
    rename(psamp_mean_q5 = `q5`,
           psamp_mean_q50 = `mean`,
           psamp_mean_q95 = `q95`) 
  cbind(t1,t2) %>% 
    mutate(mean_y_samp = mean(sample_data$y),
           coverage = ifelse(pred_mean_q5 <= mean_y_samp & pred_mean_q95 >= mean_y_samp, 1, 0),
           iter = iteration,
           n_level = nl,
           samp_size = sample_size,
           rhat_larger_1_1 = sum(modelfit$summary()[,'rhat']>1.1), # diagnostics
           n_divergent = sum(modelfit$diagnostic_summary()['num_divergent'][[1]]),
           n_maxtreedepth = sum(modelfit$diagnostic_summary()['num_max_treedepth'][[1]]),
           ess_bulk =  modelfit$summary()[,'ess_bulk'] %>% min(),
           est_int = modelfit$summary('b[1]')[2] %>% as.numeric() %>% inv_logit_scaled(),
           spec = spec_values,
           TN = TN_value, 
           FP = FP_value)
}
