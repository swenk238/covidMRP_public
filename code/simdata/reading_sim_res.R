library(tidyverse)
library(purrr)
library(here)

## reading in simdata results 
func_wd <- 'code/simdata/func/'
popn_wd <- "code/simdata/sim2_spec/data/" 
stan_wd <- 'code/simdata/sim2_spec/stan/'

# working directory
res_sim1_int_wd <- "results/simdata/sim1/intercept_only/"
res_sim1_cov_wd <- 'results/simdata/sim1/covariates/'
popn_wd <- "code/simdata/sim2_spec/data/"
res_sim2_fs_wd <- "results/simdata/sim2/noEstSpec/fixedSpec/"
res_sim2_es_wd <- "results/simdata/sim2/estSpec/"

res_sim3_wd <- "results/simdata/sim3/"


# Sim 1.1 loading the two RDS files ---------------------------
sim1_mat = readRDS(file=here(paste0(res_sim1_int_wd, "res_mat_all.RDS"))) %>% 
  do.call(rbind, .)
sim1_y_popn = readRDS(file=here(paste0(res_sim1_int_wd, "y_popn_all.RDS")))[[1]] # taking first iter as all the population across iteration is identical


# Sim 1.2 covariates --------------------------------------
# mean of y for two population for different beta1_set
popn_data_sim1_cov_list = readRDS(file=here(paste0(res_sim1_cov_wd,'cluster/popn_data_sim1_cov.RDS')))

## reading all 100 iteration results
# res_sim1_cov_list <- res_sae_sim1_cov_list <- samp_data_sim1_ss400 <- samp_data_sim1_ss4000 <- list()
# for (i in 1:100){
#   samp_data_sim1_ss400[[i]] <-  readRDS(file=here(paste0(res_sim1_cov_wd, "cluster/samp_data/samp_data_ss400_list_ite",i,".RDS")))
#   samp_data_sim1_ss4000[[i]] <-  readRDS(file=here(paste0(res_sim1_cov_wd, "cluster/samp_data/samp_data_ss4000_list_ite",i,".RDS")))
# 
#   res_sim1_cov_list[[i]]  <- readRDS(file=here(paste0(res_sim1_cov_wd, "cluster/res_mat_ite",i,".RDS"))) %>%
#     mutate(iter = i)
# 
#   res_sae_sim1_cov_list[[i]] <- readRDS(file=here(paste0(res_sim1_cov_wd, "cluster/res_sae_mat_ite",i,".RDS"))) |>
#     mutate(iter = i)
# }
# 
# # identify which iterations did not have complete estimates (only samp size 400 gives issue)
# lapply(res_sim1_cov_list, function(x)nrow(x)) %>% do.call(rbind,.) %>% as_tibble() %>% mutate(iter = row_number()) %>% filter(V1 < 96) %>% pull(iter) -> sel_iter
# 
# # missing
# missing_config <- res_sim1_cov_list[sel_iter] %>%
#   lapply(., function(x)x %>% count(n_level, beta1_set, samp_size))
# 
# # second iteration has full configuration
# full_config <- res_sim1_cov_list[[2]] %>% count(n_level, beta1_set, samp_size)
# 
# # full_config is the full 16-combination row tibble
# # missing_config is the list of tibbles to check
# missing_rows <- purrr::map(missing_config, ~
#                              anti_join(full_config, .,
#                                        by = c("n_level", "beta1_set",
#                                               "samp_size", "n"))) %>%
#   do.call(rbind, .) %>%
#   mutate(iter = rep(sel_iter, each=4)) # repeat for 4 number of levels
# 
# ## running sim1_covariates_extra based on 'missing_rows' by using seed - 1
# # extra results
# res_mat_extra <- readRDS(paste0(res_sim1_cov_wd, 'cluster/res_mat_extra_ite.RDS'))
# res_sae_mat_extra <- readRDS(paste0(res_sim1_cov_wd, 'cluster/res_sae_mat_extra_ite.RDS'))
# samp_data_ss400_extra <- readRDS(paste0(res_sim1_cov_wd, 'cluster/samp_data/samp_data_ss400_list_extra_ite.RDS'))
# 
# lapply(samp_data_ss400_extra, function(x) x %>% 
#          group_by(n_level, beta1_set, samp_size, iter) %>% 
#          summarise(pi_sampdata = mean(y),
#                    p_sampdata = mean(p_y))) %>% 
#   do.call(rbind, .) -> ss400_extra
# 
# do.call(rbind, res_sim1_cov_list) %>%
#   rbind(., res_mat_extra) %>%
#   mutate(pi_popndata_truth = case_when(beta1_set == 0.3 ~ mean(popn_data_sim1_list[[1]]$y),
#                                        beta1_set == 0 ~ mean(popn_data_sim1_list[[2]]$y),
#                                        TRUE ~ NA_real_),
#          p_popndata_truth = case_when(beta1_set == 0.3 ~ mean(popn_data_sim1_list[[1]]$p_y),
#                                       beta1_set == 0 ~ mean(popn_data_sim1_list[[2]]$p_y),
#                                       TRUE ~ NA_real_),
#          pi_sampdata = mean_y_samp, # saved in the iterations as mean_y_samp
#          p_sampdata = pmap_dbl(list(samp_size, beta1_set, iter), ~{
#            samp_data <- get(paste0("samp_data_sim1_ss", ..1))
#            beta1_index <- ifelse(..2 == 0.3, 1, 2)
#            mean(samp_data[[..3]][[beta1_index]]$p_y)}),
#          post_intercept_bias = est_int - 0.01) %>% 
#   select(-mean_y_samp) %>% 
#   rows_update(ss400_extra, by=c('n_level', 'beta1_set', 'samp_size', 'iter')) %>%  #update the sample data results for missing iter
#   saveRDS(., file=paste0(res_sim1_cov_wd,'cluster/res_sim1_cov_mat_complete.RDS'))
# 
# do.call(rbind, res_sae_sim1_cov_list) %>%
#   rbind(., res_sae_mat_extra)  %>%
#   mutate(pi_popndata_truth = case_when(beta1_set == 0.3 ~ mean(popn_data_sim1_cov_list[[1]]$y),
#                                        beta1_set == 0 ~ mean(popn_data_sim1_cov_list[[2]]$y),
#                                        TRUE ~ NA_real_),
#          p_popndata_truth = case_when(beta1_set == 0.3 ~ mean(popn_data_sim1_cov_list[[1]]$p_y),
#                                       beta1_set == 0 ~ mean(popn_data_sim1_cov_list[[2]]$p_y),
#                                       TRUE ~ NA_real_),
#          pi_sampdata = pmap_dbl(list(samp_size, beta1_set, iter), ~ {
#            samp_data <- get(paste0("samp_data_sim1_ss", ..1))
#            beta1_index <- ifelse(..2 == 0.3, 1, 2)
#            mean(samp_data[[..3]][[beta1_index]]$y) # getting y_pi_y according to sample size and beta1
#          }),
#          p_sampdata = pmap_dbl(list(samp_size, beta1_set, iter), ~{
#            samp_data <- get(paste0("samp_data_sim1_ss", ..1))
#            beta1_index <- ifelse(..2 == 0.3, 1, 2)
#            mean(samp_data[[..3]][[beta1_index]]$p_y)})) %>%  # getting mean of y* according to sample size and beta1_set
#   rows_update(ss400_extra, by=c('n_level', 'beta1_set', 'samp_size', 'iter')) %>%  #update the sample data results for missing iter
#   saveRDS(., file=paste0(res_sim1_cov_wd,'cluster/res_sae_sim1_cov_mat_complete.RDS'))

res_sim1_cov_mat <- readRDS(file=here(paste0(res_sim1_cov_wd, "cluster/res_sim1_cov_mat_complete.RDS")))  # results for 2 beta1_set, 2 samp_size, 4 n_levels, 6 models
res_sae_sim1_cov_mat <- readRDS(file=here(paste0(res_sim1_cov_wd,"cluster/res_sae_sim1_cov_mat_complete.RDS"))) 


# Sim 2.1 noEstSpec -------------------------------------------------------------------
popn_data_TNFP = readRDS(file=here(paste0(popn_wd,"popn_data_TNFP_stableN.rds")))
mean_y_popn_vec = apply(popn_data_TNFP %>%dplyr::select(y_spec0.980:y_spec1.000), 2, function(x)mean(x))
samp_data_nes_ss400_list <- readRDS(file=here(paste0(res_sim2_fs_wd,"samp_data_ss400_list.RDS")))
samp_data_nes_ss4000_list <- readRDS(file=here(paste0(res_sim2_fs_wd,"samp_data_ss4000_list.RDS")))


# # reading in all iterations
# res_sim2_nes_mat <- res_sae_sim2_nes_mat <- list()
# for (i in 1:100){
#   res_sim2_nes_mat[[i]] <- readRDS(file=here(paste0(res_sim2_fs_wd, "res_mat_ite",i,".RDS"))) |>
#     mutate(iter = i,
#            pi_popndata_truth = mean(popn_data_TNFP$y_pi_y_gen),
#            p_spec_popndata_truth = mean_y_popn_vec[paste0("y_spec", format(spec, nsmall = 3))], # uses spec column to choose which y_spec
#            pi_sampdata = map_dbl(samp_size, ~ mean(get(paste0("samp_data_nes_ss", .x, "_list"))[[i]]$y_pi_y_gen)), # getting y_pi_y according to sample size
#            p_spec_sampdata = map2_dbl(samp_size, spec, ~{
#              samp_data <- get(paste0("samp_data_nes_ss", .x, "_list"))
#              spec_col <- paste0("y_spec", format(.y, nsmall = 3))
#              mean(samp_data[[i]][[spec_col]])}), # getting mean of y* according to sample size and specificity
#            post_intercept_bias = est_int - 0.01) # bias for intercept
# 
#   res_sae_sim2_nes_mat[[i]] <- readRDS(file=here(paste0(res_sim2_fs_wd, "res_sae_mat_ite",i,".RDS"))) |>
#     mutate(iter = i,
#            pi_truth = mean(popn_data_TNFP$y_pi_y_gen),
#            p_truth_spec = mean_y_popn_vec[paste0("y_spec", format(spec, nsmall = 3))], # uses spec column to choose which y_spec
#            pi_sample = map_dbl(samp_size, ~ mean(get(paste0("samp_data_nes_ss", .x, "_list"))[[i]]$y_pi_y_gen)), # getting y_pi_y according to sample size
#            p_spec_sample = map2_dbl(samp_size, spec, ~{
#              samp_data <- get(paste0("samp_data_nes_ss", .x, "_list"))
#              spec_col <- paste0("y_spec", format(.y, nsmall = 3))
#              mean(samp_data[[i]][[spec_col]])}))
# }
# 
# 
# res_sim2_nes_mat %>% do.call(rbind,.) %>% select(-mean_y_samp) %>% saveRDS(., file=here(paste0(res_sim2_fs_wd, "res_mat_all.RDS")))
# res_sae_sim2_nes_mat %>% do.call(rbind, .) %>% saveRDS(., file=here(paste0(res_sim2_fs_wd,"res_sae_mat_all.RDS")))

# fixed spec 
res_sim2_nes_fixedspec_all <- readRDS(file=here(paste0(res_sim2_fs_wd, "res_mat_all.RDS"))) 
res_sae_sim2_nes_fixedspec_all <- readRDS(file=here(paste0(res_sim2_fs_wd, "res_sae_mat_all.RDS"))) 

# Sim 2.2 EstSpec -----------------------------------------------------------
# res_sim2_es_mat_all <- readRDS(file=here(paste0(res_sim2_es_wd,"res_mat_list_all.RDS"))) %>%
#   do.call(rbind, .)
# res_sae_sim2_es_mat_all = readRDS(file=here(paste0(res_sim2_es_wd,"res_sae_mat_list_all.RDS"))) %>%
#   do.call(rbind, .)
# 
# ## resaving the popn and sample truth
# res_sim2_es_mat <- res_sim2_es_mat_all %>%
#   select(-c(mean_y_samp, mean_y_popn:est_int_truth_bias)) %>%
#   mutate(pi_popndata_truth = mean(popn_data_TNFP$y_pi_y_gen),
#          p_spec_popndata_truth = mean_y_popn_vec[paste0("y_spec", format(spec, nsmall = 3))], # uses spec column to choose which y_spec
#          pi_sampdata = map2_dbl(samp_size, iter, ~ {
#            mean(get(paste0("samp_data_nes_ss", .x, "_list"))[[.y]]$y_pi_y_gen)}), # getting y_pi_y according to sample size
#          p_spec_sampdata = pmap_dbl(list(samp_size, spec, iter), ~{
#            samp_data <- get(paste0("samp_data_nes_ss", ..1, "_list"))
#            spec_col <- paste0("y_spec", format(..2, nsmall = 3))
#            mean(samp_data[[..3]][[spec_col]])}), # getting mean of y* according to sample size and specificity
#          post_intercept_bias = est_int - 0.01) # bias for intercept
# 
# res_sae_sim2_es_mat <- res_sae_sim2_es_mat_all %>%
#   select(-c(mean_y_popn, mean_y_pi_y_gen)) %>%
#   mutate(pi_popndata_truth = mean(popn_data_TNFP$y_pi_y_gen),
#          p_spec_popndata_truth = mean_y_popn_vec[paste0("y_spec", format(spec, nsmall = 3))], # uses spec column to choose which y_spec
#          pi_sampdata = map2_dbl(samp_size, iter, ~ {
#            mean(get(paste0("samp_data_nes_ss", .x, "_list"))[[.y]]$y_pi_y_gen)}), # getting y_pi_y according to sample size
#          p_spec_sampdata = pmap_dbl(list(samp_size, spec, iter), ~{
#            samp_data <- get(paste0("samp_data_nes_ss", ..1, "_list"))
#            spec_col <- paste0("y_spec", format(..2, nsmall = 3))
#            mean(samp_data[[..3]][[spec_col]])})) # getting mean of y* according to sample size and specificity
# 
# res_sim2_es_mat %>% saveRDS(., file=here(paste0(res_sim2_es_wd,"res_mat_list.RDS")))
# res_sae_sim2_es_mat %>% saveRDS(., file=here(paste0(res_sim2_es_wd,"res_sae_mat_list.RDS")))

res_sim2_es_mat <- readRDS(file=here(paste0(res_sim2_es_wd,"res_mat_list.RDS")))
res_sae_sim2_es_mat <- readRDS(file=here(paste0(res_sim2_es_wd,"res_sae_mat_list.RDS")))

# Sim 3 -------------------------------------------------------------------
# FE results
samp_data_sim3_ss400_list <- samp_data_sim3_ss4000_list <- list()

for (i in 1:100){
    samp_data_sim3_ss400_list[[i]] <-  readRDS(file=here(paste0(res_sim3_wd, "bprior/samp_data/samp_data_ss400_ite",i,".RDS")))
    samp_data_sim3_ss4000_list[[i]] <-  readRDS(file=here(paste0(res_sim3_wd, "bprior/samp_data/samp_data_ss4000_ite",i,".RDS")))
}

res_sim3_FE_mat <- readRDS(file=here(paste0(res_sim3_wd,"bprior/res_FE_mat_list_all.RDS"))) %>%
  select(-c(mean_y_samp, mean_y_popn:est_int_truth_bias)) %>%
  mutate(pi_popndata_truth = mean(popn_data_TNFP$y_pi_y_gen),
         p_spec_popndata_truth = mean_y_popn_vec[paste0("y_spec", format(spec, nsmall = 3))], # uses spec column to choose which y_spec
         pi_sampdata = pmap_dbl(list(samp_size, iter, spec), function(samp_size, iter, spec){
           samp_data <- get(paste0("samp_data_sim3_ss", samp_size, "_list"))
           spec_index <- case_when(
             spec == '0.980' ~ 1,
             spec == '0.990' ~ 5,
             spec == '0.995' ~ 9,
             spec == '1.000' ~ 13,
             TRUE ~ NA_integer_
           )
           mean(samp_data[[iter]][[spec_index]]$y_pi_y_gen)}), # getting y_pi_y according to sample size
         p_spec_sampdata = pmap_dbl(list(samp_size, iter, spec),  function(samp_size, iter, spec){
           samp_data <- get(paste0("samp_data_sim3_ss", samp_size, "_list"))
           spec_index <- case_when(
             spec == '0.980' ~ 1,
             spec == '0.990' ~ 5,
             spec == '0.995' ~ 9,
             spec == '1.000' ~ 13,
             TRUE ~ NA_integer_
           )
           mean(samp_data[[iter]][[spec_index]]$y)}), # getting mean of y* according to sample size and specificity
         post_intercept_bias = est_int - 0.01) # bias for intercept


results_sim_list = list()
results_sim_list = list(sim1_mat = sim1_mat,
                        sim1_y_popn = sim1_y_popn,
                        res_sim1_cov_mat = res_sim1_cov_mat,
                        res_sae_sim1_cov_mat = res_sae_sim1_cov_mat, 
                        popn_data_sim1_cov_list = popn_data_sim1_cov_list,
                        popn_data_TNFP = popn_data_TNFP,
                        res_sim2_nes_fixedspec_all = res_sim2_nes_fixedspec_all,
                        res_sae_sim2_nes_fixedspec_all = res_sae_sim2_nes_fixedspec_all,
                        res_sim2_es_mat = res_sim2_es_mat,
                        res_sae_sim2_es_mat = res_sae_sim2_es_mat,
                        res_sim3_FE_mat = res_sim3_FE_mat)

save(results_sim_list, file=here("results/results_sim_list.RData"), compress=T)
