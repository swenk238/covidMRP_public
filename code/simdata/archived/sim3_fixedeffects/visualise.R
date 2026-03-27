# visualising -------------------------------------------------------------
# simdata
load("~/GitHub/covidMRP/results/simdata/sim3/res_mat_allFE_priorpred.RData")

res_mat = lapply(list(res_m3_FE, res_m4_FE, res_m5_FE,
            res_m3_twoFE, res_m4_twoFE, res_m5_twoFE), function(x)x$draws(variables=c("p_avg1","y_rep")))

load("~/GitHub/covidMRP/results/simdata/sim3/R2D2M2/res_priorpred_ite1.RData")
res_R2D2_mat = lapply(list(res_m3_FE_R2D2M2, res_m4_FE_R2D2M2, res_m5_FE_R2D2M2,
                           res_m3_twoFE_R2D2M2, res_m4_twoFE_R2D2M2, res_m5_twoFE_R2D2M2), function(x)x$draws(variables=c("p_avg1","y_rep")))





p_prior <- res_priorpred_m5_r2prior$draws(variables = "p_avg1")
ori_prior <- res_prior_m5$draws(variables = "p_avg1")[c(seq(1,nrow(ori_prior), by=1000)),,]


ori_prior  %>%
  as.data.frame() %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(., aes(x = value)) +
  geom_density(alpha = 0.5, fill = "blue") +
  labs(title = "Prior Predictive Distribution (ori model)", x = "p_avg1", y = "Density")

p_prior  %>%
  as.data.frame() %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(., aes(x = value)) +
  geom_density(alpha = 0.5, fill = "blue") +
  labs(title = "Prior Predictive Distribution (R2 prior)", x = "p_avg1", y = "Density")
