library(here)
library(tidyverse)

# reading in all the results
res_ori_wd <- 'results/realdata/original/'
res_noFEnoME_wd <- 'results/realdata/sim1_noFE_noME/'
res_nopsamp_wd <- 'results/realdata/sim1_noFE_noME/no_p_samp/'
res_sim2_wd <- 'results/realdata/sim2_spec/'
res_seifaInv_wd <- 'results/realdata/sim3_noFE/'

# original MRP estimates ----------------------------------------------------------------
res_posterior = list()
for(i in 1:8){
  res_posterior[[i]] =  readRDS(file=here(paste0(res_ori_wd, "m",i-1,"_postpred.RDS")))
}

res_est_full = lapply(res_posterior, function(x)x$summary(c("p_avg1","b[1]"))[,c('mean', 'q5', 'q95')]) %>% 
  do.call(rbind, .) |> 
  mutate(model = rep(0:7,each=2),
         modification = "Original", 
         type=rep(c("predicted estimate", "intercept"),8)) 
res_est = res_est_full |> filter(model %in% 0:5)


# # Sim 1.1 equivalent  ------------------------------------------------------------
# no ME no FE models 
# results for p_samp (measurement error) and FEs removed for models 3-5
res_noFE_noME_list = list()
for(i in 1:3){
  res_noFE_noME_list[[i]] <- readRDS(file=here(paste0(res_noFEnoME_wd,"res_m",i+2,"a_noFE_noME.RDS")))
}

res_noME_noFE <- lapply(res_noFE_noME_list, function(x)x$summary(c("p_avg1", "b[1]"))[,c('mean', 'q5', 'q95')]) %>%
  do.call(rbind, .) |> 
  mutate(model = rep(3:5,each=2),
         modification = "No", 
         type=rep(c("predicted estimate", "intercept"), 3))

# results for p_samp (measurement error) removed for all models
res_nopsamp_list = list()
for(i in 1:7){
  res_nopsamp_list[[i]] = readRDS(file=here(paste0(res_nopsamp_wd,"m",i-1,"nopsamp.RDS")))
}

res_est_nopsamp = lapply(res_nopsamp_list, function(x)x$summary(c("p_avg1", "b[1]"))[,c('mean', 'q5', 'q95')]) %>%
  do.call(rbind, .) |> 
  mutate(model = rep(0:6,each=2),
         modification = "No", 
         type=rep(c("predicted estimate", "intercept"), 7))


# Sim 2.1 equivalent ------------------------------------------------------------
## fixed spec ####
# spec = 0.980
res_0.980_list = list()
for(i in 1:6){
  res_0.980_list[[i]] = readRDS(file=here(paste0(res_sim2_wd,"model",i-1,"spec0.98sens1_noFE.RDS")))
}

res_est_0.980spec = lapply(res_0.980_list, function(x)x$summary(c("p_avg1", "b[1]"))[,c('mean', 'q5', 'q95')]) %>%
  do.call(rbind,.) |> 
  mutate(model = rep(0:5, each=2),
         modification = "0.980", 
         type=rep(c("predicted estimate", "intercept"), 6))

# fixed spec = 0.99
res_0.990_list = list()
for(i in 1:6){
  res_0.990_list[[i]] = readRDS(file=here(paste0(res_sim2_wd,"model",i-1,"spec0.99sens1_noFE.RDS")))
}

res_est_0.990spec = lapply(res_0.990_list, function(x)x$summary(c("p_avg1", "b[1]"))[,c('mean', 'q5', 'q95')]) %>%
  do.call(rbind,.) |> 
  mutate(model = rep(0:5, each=2),
         modification = "0.990", 
         type=rep(c("predicted estimate", "intercept"), 6))

# fixed spec = 0.995 
res_0.995_list = list()
for(i in 1:6){
  res_0.995_list[[i]] = readRDS(file=here(paste0(res_sim2_wd,"model",i-1,"spec0.995sens1_noFE.RDS")))
}

res_est_0.995spec = lapply(res_0.995_list, function(x)x$summary(c("p_avg1", "b[1]"))[,c('mean', 'q5', 'q95')]) %>%
  do.call(rbind,.) |> 
  mutate(model = rep(0:5, each=2),
         modification = "0.995", 
         type=rep(c("predicted estimate", "intercept"), 6))

# fixed spec=1
res_1.000_list = list()
for(i in 1:6){
  res_1.000_list[[i]] = readRDS(file=here(paste0(res_sim2_wd,"model",i-1,"spec1sens1_noFE.RDS")))
}

res_est_1.000spec = lapply(res_1.000_list, function(x)x$summary(c("p_avg1", "b[1]"))[,c('mean', 'q5', 'q95')]) %>%
  do.call(rbind,.)  |> 
  mutate(model = rep(0:5, each=2),
         modification = "1.000", 
         type=rep(c("predicted estimate", "intercept"), 6))

# Sim 2.2 equivalent ------------------------------------------------------
## estimating specificity ####
# 16 different TN, FP values
total_vec <- rep(c(400,800,1200,8000),4)
FP_vec <- c(8,16,24,160, 4,8,12,80, 2,4,6,40, rep(0,4))
TN_vec <- total_vec - FP_vec # corresponding to specificity values of 0.98, 0.99, 0.995 and 1 

## reading in the results for various TN and FP
res_est_fixedsens_noFE <- tibble()
for(i in 1:6){
  for(j in 1:length(TN_vec)){
  file <- readRDS(file=here(paste0(res_sim2_wd,"m",i-1,"_TN", TN_vec[j], "FP", FP_vec[j],"fixedsens_noFE.RDS"))) 
  res_tb <-  file$summary(c("p_avg1", "b[1]"))[,c('mean', 'q5', 'q95')] %>% 
    mutate(model = i-1, 
           type= c("predicted estimate", "intercept"),
           TN = TN_vec[j], # repeat for 6 models
           FP = FP_vec[j],
           modification = "Est spec fixed sens", 
           N_xi = TN + FP,
           spec = TN / N_xi)
  res_est_fixedsens_noFE <- rbind(res_est_fixedsens_noFE, res_tb) # appending results 
  }
}

## estimating sens and spec as well
res_est_estsens_noFE <- tibble()
for(i in 1:6){
  for(j in 1:length(TN_vec)){
    file <- readRDS(file=here(paste0(res_sim2_wd,"estsens/m",i-1,"_TN", TN_vec[j], "FP", FP_vec[j],"estsens_noFE.RDS"))) 
    res_tb <-  file$summary(c("p_avg1", "b[1]"))[,c('mean', 'q5', 'q95')] %>% 
      mutate(model = i-1, 
             type= c("predicted estimate", "intercept"),
             TN = TN_vec[j], # repeat for 6 models
             FP = FP_vec[j],
             modification = "Est spec est sens", 
             N_xi = TN + FP,
             spec = TN / N_xi)
    res_est_estsens_noFE <- rbind(res_est_estsens_noFE, res_tb) # appending results 
  }
}

# Sim 3 equivalent --------------------------------------------------------
## no fixed effects
res_m3a <-  readRDS(file=here(paste0(res_seifaInv_wd,"res_m3a_rmSeifaCont.RDS")))
res_m4d <- readRDS(file=here(paste0(res_seifaInv_wd,"res_m4d_rmFE.RDS")))
res_m5d <- readRDS(file=here(paste0(res_seifaInv_wd,"res_m5d_rmFE.RDS")))

lapply(list(res_m3a, res_m4d, res_m5d), function(x)x$summary(c("p_avg1", 'b[1]'))[,c('mean', 'q5', 'q95')]) %>% # combine duplicate of res_est (models 0-3), and oneFE and noFE for models 3-5
  do.call(rbind, .) %>% 
  mutate(FEornot = rep('No',3*2),
         model = rep(c(3:5), each=2), 
         modification = 'both FE removed',
         type= rep(c("predicted estimate", "intercept"),3)) -> bothFE_rm_mat

## one fixed effect
res_m4e <- readRDS(file=here(paste0(res_seifaInv_wd,"res_m4e_rmSex.RDS")))
res_m5e <- readRDS(file=here(paste0(res_seifaInv_wd,"res_m5e_rmSex.RDS")))


lapply(list(res_m4e, res_m5e), function(x)x$summary(c("p_avg1", 'b[1]'))[,c('mean', 'q5', 'q95')]) %>% # combine duplicate of res_est (models 0-3), and oneFE and noFE for models 3-5
  do.call(rbind, .) %>% 
  mutate(FEornot = rep('One',2*2),
         model = rep(c(4,5),each=2), 
         modification = 'sex var. removed',
         type= rep(c("predicted estimate", "intercept"),2)) -> oneFE_rm_mat

res_est %>%
  mutate(FEornot = 'Two') %>% 
  rbind(oneFE_rm_mat) %>% 
  rbind(bothFE_rm_mat) -> res_est_FErm


# R2D2 --------------------------------------------------------------------
# R2_mean0.1_prec20
file_mat_r2_0.1_mat <- here(paste0(res_ori_wd, "R2_0.1/m", 0:7, "_r2_0.1_postpred.RDS"))

res_r2_0.1_mat <- lapply(1:8, function(iter){
  readRDS(file_mat_r2_0.1_mat[iter])$summary('p_avg1')[c('mean', 'q5', 'q95')] %>% 
    mutate(model = iter-1)}) %>% 
  do.call(rbind, .) %>% 
  mutate(modification = "R2(mean, prec) = (0.1, 20)")

# 'base' R2_mean0.5_prec2
file_mat_r2_0.5_mat <- here(paste0(res_ori_wd, "R2_0.5/m", 0:7, "_r2_0.5_postpred.RDS"))
res_r2_0.5_mat <- lapply(1:8, function(iter){
  readRDS(file_mat_r2_0.5_mat[iter])$summary('p_avg1')[c('mean', 'q5', 'q95')] %>% 
    mutate(model = iter-1)}) %>% 
  do.call(rbind, .) %>% 
  mutate(modification = "R2(mean, prec) = (0.5, 2)")

# Aki's recommendation
# 'base'
file_mat_r2_0.3_mat <- here(paste0(res_ori_wd, "R2_0.3/m", 0:7, "_r2_0.3_postpred.RDS"))
res_r2_0.3_mat <- lapply(1:8, function(iter){
  readRDS(file_mat_r2_0.3_mat[iter])$summary('p_avg1')[c('mean', 'q5', 'q95')] %>% 
    mutate(model = iter-1)}) %>% 
  do.call(rbind, .) %>% 
  mutate(modification = "R2(mean, prec) = (0.333, 3)")

results_list = list(res_est = res_est, 
                    res_est_full = res_est_full,
                    res_noME_noFE = res_noME_noFE,
                    res_est_0.980spec = res_est_0.980spec,
                    res_est_0.990spec = res_est_0.990spec,
                    res_est_0.995spec = res_est_0.995spec,
                    res_est_1.000spec = res_est_1.000spec, 
                    res_est_fixedsens_noFE = res_est_fixedsens_noFE,
                    res_est_estsens_noFE = res_est_estsens_noFE,
                    res_est_nopsamp = res_est_nopsamp,
                    res_est_FErm = res_est_FErm,
                    res_r2_0.5_mat = res_r2_0.5_mat,
                    res_r2_0.3_mat = res_r2_0.3_mat,
                    res_r2_0.1_mat = res_r2_0.1_mat)

save(results_list, file=here("results/results_list.RData"), compress=T)
