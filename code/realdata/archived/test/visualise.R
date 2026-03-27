## visualising results 
library(here)
source(here("code/data_list.R"), echo=F)

res_prior_m1 = readRDS(file=here("results/model files/m1_priorpred.RDS"))
res_post_m1 = readRDS(file=here("results/model files/m1_postpred.RDS"))

res_prior_m2 = readRDS(file=here("results/model files/m2_priorpred.RDS"))
res_post_m2 = readRDS(file=here("results/model files/m2_postpred.RDS"))

res_prior_m3 = readRDS(file=here("results/model files/m3_priorpred.RDS"))
res_post_m3 = readRDS(file=here("results/model files/m3_postpred.RDS"))

res_prior_m4 = readRDS(file=here("results/model files/m4_priorpred.RDS"))
res_post_m4 = readRDS(file=here("results/model files/m4_postpred.RDS"))

res_prior_m5 = readRDS(file=here("results/model files/m5_priorpred.RDS"))
res_post_m5 = readRDS(file=here("results/model files/m5_postpred.RDS"))

# loading calculating yrep
load("~/GitHub/covidMRP/results/model files/yrep_prior_post_mat.RData")


# model 1 -----------------------------------------------------------------
# yrep and y plots
priormean_m1 = as.numeric(res_prior_m1$summary("p_avg1")['mean'])*100
priorquant_m1 = as.numeric(res_prior_m1$summary("p_avg1")[c("q5","q95")])*100 
y = data_list$y
truemean = mean(y)*100 # in percentages

SEQ = seq(1,4000, length.out=1000)
color_scheme_set("brightblue")
p1_m1 = ppc_dens_overlay(yrep_m1_prior_mat[SEQ, ]) +
  labs(title="Prior predictive densities - model1") +
  legend_move(c(0.95,0.5)) 

# MRP estimation 
c1_m1 = mcmc_intervals(res_prior_m1$draws(), pars=c("p_avg1"), transformations = function(x)x*100) + 
  scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+0.5, y=1.25, label="truth")) +
    geom_text(aes(x=priormean_m1+0.05, y=0.7,  label=paste0("MRP estimate = \n", round(priormean_m1,2), "% (", round(priorquant_m1[1],3), ", ", round(priorquant_m1[2],3), ")")), col="blue") 
  
( plot_m1_prior = bayesplot_grid(p1_m1, c1_m1, grid_args = list(nrow=2, heights=c(7,2))) )

ggsave(file=here("results/plot_m1_prior.png"), plot_m1_prior, device = "png")

# posterior pred
postmean_m1 = as.numeric(res_post_m1$summary("p_avg1")['mean'])*100 
postquant_m1 = as.numeric(res_post_m1$summary("p_avg1")[c("q5","q95")])*100 

color_scheme_set("purple")
p2_m1 = ppc_dens_overlay(y, yrep_m1_post_mat[SEQ, ]) +
  labs(title="Posterior predictive densities - model1") +
  legend_move(c(0.95,0.5)) 

# MRP estimation
( c2_m1 = mcmc_intervals(res_post_m1$draws(), pars=c("p_avg1"), point_size=15, outer_size=13, inner_size=13, transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+0.02, y=1.25, label="truth")) +
    geom_text(aes(x=postmean_m1+12, y=1, label=paste0("MRP estimate = \n", round(postmean_m1,2), "% (", round(postquant_m1[1],3), ", ", round(postquant_m1[2],3), ")")), col="purple") )

( plot_m1_post = bayesplot_grid(p2_m1, c2_m1, grid_args = list(nrow=2, heights=c(7,2))) )

ggsave(file=here("results/plot_m1_post.png"), plot_m1_post, device = "png")

# model 2 -----------------------------------------------------------------
# prior pred 
priormean_m2 = as.numeric(res_prior_m2$summary("p_avg1")['mean'])*100
priorquant_m2 = as.numeric(res_prior_m2$summary("p_avg1")[c("q5","q95")])*100 
y = data_list$y
truemean = mean(y)100

# yrep and y plots
color_scheme_set("brightblue")
SEQ = seq(1,4000, length.out=1000)
p1_m2 = ppc_dens_overlay(y, yrep_m2_prior_mat[SEQ, ]) + # thinned
  labs(title="Prior predictive densities - model2") +
  legend_move(c(0.95,0.5)) 

# MRP estimation 
( c1_m2 = mcmc_intervals(res_prior_m2$draws(), pars=c("p_avg1"), transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+0.02, y=1.25, label="truth")) +
    geom_text(aes(x=priormean_m2+0.15, y=0.7,  label=paste0("MRP estimate = \n", round(priormean_m2,2), "% (", round(priorquant_m2[1],3), ", ", round(priorquant_m2[2],3), ")")), col="blue") )

( plot_m2_prior = bayesplot_grid(p1_m2, c1_m2, grid_args = list(nrow=2, heights=c(7,2))) )

ggsave(file=here("results/plot_m2_prior.png"), plot_m2_prior, device = "png")

# posterior pred - constructing yrep matrix
postmean_m2 = as.numeric(res_post_m2$summary("p_avg1")['mean'])*100 
postquant_m2 = as.numeric(res_post_m2$summary("p_avg1")[c("q5","q95")])*100 

color_scheme_set("purple")
( p2_m2 = ppc_dens_overlay(y, yrep_m2_post_mat[SEQ, ]) +
  labs(title="Posterior predictive densities - model2") +
  legend_move(c(0.95,0.5)) )

# MRP estimation
( c2_m2 = mcmc_intervals(res_post_m2$draws(), pars=c("p_avg1"), point_size=15, outer_size=13, inner_size=13, transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+0.02, y=1.25, label="truth")) +
    geom_text(aes(x=postmean_m2+12, y=1, label=paste0("MRP estimate = \n", round(postmean_m2,2), "% (", round(postquant_m2[1],3), ", ", round(postquant_m2[2],3), ")")), col="purple") )

( plot_m2_post = bayesplot_grid(p2_m2, c2_m2, grid_args = list(nrow=2, heights=c(7,2))) )

ggsave(file=here("results/plot_m2_post.png"), plot_m2_post, device = "png")

# model 3 -----------------------------------------------------------------
# prior pred 
priormean_m3 = as.numeric(res_prior_m3$summary("p_avg1")['mean'])*100
priorquant_m3 = as.numeric(res_prior_m3$summary("p_avg1")[c("q5","q95")])*100 

y = data_list$y
truemean = mean(y)*100

# yrep and y plots
color_scheme_set("brightblue")
SEQ = seq(1,4000, length.out=1000)
p1_m3 = ppc_dens_overlay(y, yrep_m3_prior_mat[SEQ, ]) + # thinned
  labs(title="Prior predictive densities - model3") +
  legend_move(c(0.95,0.5)) 

# priorpred - MRP estimation 
( c1_m3 = mcmc_intervals(res_prior_m3$draws(), pars=c("p_avg1"), transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+0.02, y=1.2, label="truth")) +
    geom_text(aes(x=priormean_m3+12, y=0.7,  label=paste0("MRP estimate = \n", round(priormean_m3,2), "% (", round(priorquant_m3[1],3), ", ", round(priorquant_m3[2],3), ")")), col="blue") )

plot_m3_prior = bayesplot_grid(p1_m3, c1_m3, grid_args = list(nrow=2, heights=c(7,2)))

ggsave(file=here("results/plot_m3_prior.png"), plot_m3_prior, device = "png")

# posterior pred 
postmean_m3 = as.numeric(res_post_m3$summary("p_avg1")['mean'])*100 
postquant_m3 = as.numeric(res_post_m3$summary("p_avg1")[c("q5","q95")])*100 
color_scheme_set("purple")
( p2_m3 = ppc_dens_overlay(y, yrep_m3_post_mat[SEQ, ]) +
    labs(title="Posterior predictive densities - model3") +
    legend_move(c(0.95,0.5)) )

# MRP estimation
( c2_m3 = mcmc_intervals(res_post_m3$draws(), pars=c("p_avg1"), point_size=15, outer_size=13, inner_size=13, transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+1, y=1.25, label="truth")) +
    geom_text(aes(x=postmean_m3+15, y=1, label=paste0("MRP estimate = \n", round(postmean_m3,2), "% (", round(postquant_m3[1],3), ", ", round(postquant_m3[2],3), ")")), col="purple") )

( plot_m3_post = bayesplot_grid(p2_m3, c2_m3, grid_args = list(nrow=2, heights=c(7,2))) )

ggsave(file=here("results/plot_m3_post.png"), plot_m3_post, device = "png")

# model 4 -----------------------------------------------------------------
# prior pred
priormean_m4 = as.numeric(res_prior_m4$summary("p_avg1")['mean'])*100
priorquant_m4 = as.numeric(res_prior_m4$summary("p_avg1")[c("q5","q95")])*100 
y = data_list$y
truemean = mean(y)*100

# yrep and y plots
color_scheme_set("brightblue")
SEQ = seq(1,4000, length.out=1000)
p1_m4 = ppc_dens_overlay(y, yrep_m4_prior_mat[SEQ, ]) + # thinned
  labs(title="Prior predictive densities - model4") +
  legend_move(c(0.95,0.5)) 

# MRP estimation 
( c1_m4 = mcmc_intervals(res_prior_m4$draws(), pars=c("p_avg1"), transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+1.1, y=1.25, label="truth")) +
    geom_text(aes(x=priormean_m4+0.05, y=0.65,  label=paste0("MRP estimate = \n", round(priormean_m4,2), "% (", round(priorquant_m4[1],3), ", ", round(priorquant_m4[2],3), ")")), col="blue") )

( plot_m4_prior = bayesplot_grid(p1_m4, c1_m4, grid_args = list(nrow=2, heights=c(7,2))) )

ggsave(file=here("results/plot_m4_prior.png"), plot_m4_prior, device = "png")

# posterior pred 
postmean_m4 = as.numeric(res_post_m4$summary("p_avg1")['mean'])*100 
postquant_m4 = as.numeric(res_post_m4$summary("p_avg1")[c("q5","q95")])*100 
color_scheme_set("purple")
( p2_m4 = ppc_dens_overlay(y, yrep_m4_post_mat[SEQ, ]) +
    labs(title="Posterior predictive densities - model4") +
    legend_move(c(0.95,0.5)) )

# MRP estimation
( c2_m4 = mcmc_intervals(res_post_m4$draws(), pars=c("p_avg1"), point_size=15, outer_size=13, inner_size=13, transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+1, y=1.25, label="truth")) +
    geom_text(aes(x=postmean_m4+15, y=1, label=paste0("MRP estimate = \n", round(postmean_m4,2), "% (", round(postquant_m4[1],3), ", ", round(postquant_m4[2],3), ")")), col="purple") )

( plot_m4_post = bayesplot_grid(p2_m4, c2_m4, grid_args = list(nrow=2, heights=c(7,2))) )

ggsave(file=here("results/plot_m4_post.png"), plot_m4_post, device = "png")

# model 5 -----------------------------------------------------------------
# prior pred
priormean_m5 = as.numeric(res_prior_m5$summary("p_avg1")['mean'])*100
priorquant_m5 = as.numeric(res_prior_m5$summary("p_avg1")[c("q5","q95")])*100 
y = data_list$y
truemean = mean(y)*100

# yrep and y plots
color_scheme_set("brightblue")
SEQ = seq(1,4000, length.out=1000)
p1_m5 = ppc_dens_overlay(y, yrep_m5_prior_mat[SEQ, ]) + # thinned
  labs(title="Prior predictive densities - model5") +
  legend_move(c(0.95,0.5)) 

# MRP estimation 
( c1_m5 = mcmc_intervals(res_prior_m5$draws(), pars=c("p_avg1"), transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+1, y=1.25, label="truth")) +
    geom_text(aes(x=priormean_m5+0.05, y=0.75,  label=paste0("MRP estimate = \n", round(priormean_m5,2), "% (", round(priorquant_m5[1],3), ", ", round(priorquant_m5[2],3), ")")), col="blue") )

plot_m5_prior = bayesplot_grid(p1_m5, c1_m5, grid_args = list(nrow=2, heights=c(7,2)))

ggsave(file=here("results/plot_m5_prior.png"), plot_m5_prior, device = "png")

# posterior pred 
postmean_m5 = as.numeric(res_post_m5$summary("p_avg1")['mean'])*100 
postquant_m5 = as.numeric(res_post_m5$summary("p_avg1")[c("q5","q95")])*100 

color_scheme_set("purple")
( p2_m5 = ppc_dens_overlay(y, yrep_m5_post_mat[SEQ, ]) +
    labs(title="Posterior predictive densities - model5") +
    legend_move(c(0.95,0.5)) )

# MRP estimation
( c2_m5 = mcmc_intervals(res_post_m5$draws(), pars=c("p_avg1"), point_size=15, outer_size=13, inner_size=13, transformations = function(x)x*100) + 
    scale_y_discrete(labels="") +
    scale_x_continuous(limits=c(0,100), expand=c(0,0)) +
    geom_segment(aes(x=truemean, xend=truemean, y=0.1, yend=1.2), linewidth=1.5, col="black") +
    geom_text(aes(x=truemean+1.1, y=1.25, label="truth")) +
    geom_text(aes(x=postmean_m5+15, y=1, label=paste0("MRP estimate = \n", round(postmean_m5,2), "% (", round(postquant_m5[1],3), ", ", round(postquant_m5[2],3), ")")), col="purple") )

( plot_m5_post = bayesplot_grid(p2_m5, c2_m5, grid_args = list(nrow=2, heights=c(7,2))) )

ggsave(file=here("results/plot_m5_post.png"), plot_m5_post, device = "png")


# # extracting y_pred and saving them 
# # m1 prior pred - constructing yrep matrix
# yrep_m1_prior_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m1_prior_mat[,i] = res_prior_m1$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# # m1 posterior pred - constructing yrep matrix
# yrep_m1_post_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m1_post_mat[,i] = res_post_m1$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# # m2 prior pred - constructing yrep matrix
# yrep_m2_prior_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m2_prior_mat[,i] = res_prior_m2$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# # m2 posterior pred - constructing yrep matrix
# yrep_m2_post_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m2_post_mat[,i] = res_post_m2$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# # m3 prior pred - constructing yrep matrix
# yrep_m3_prior_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m3_prior_mat[,i] = res_prior_m3$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# m3 posterior pred - constructing yrep matrix
yrep_m3_post_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
for(i in 1:length(data_list$y)){
  yrep_m3_post_mat[,i] = res_post_m3$draws(paste0("y_rep[",i,"]"), format="matrix")
}
# 
# # m4 prior pred - constructing yrep matrix
# yrep_m4_prior_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m4_prior_mat[,i] = res_prior_m4$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# # m4 posterior pred - constructing yrep matrix
# yrep_m4_post_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m4_post_mat[,i] = res_post_m4$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# # m5 prior pred - constructing yrep matrix
# yrep_m5_prior_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m5_prior_mat[,i] = res_prior_m5$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# # m5 posterior pred - constructing yrep matrix
# yrep_m5_post_mat = matrix(NA, nrow=4000, ncol=length(data_list$y))
# for(i in 1:length(data_list$y)){
#   yrep_m5_post_mat[,i] = res_post_m5$draws(paste0("y_rep[",i,"]"), format="matrix")
# }
# 
# save(yrep_m1_prior_mat, yrep_m1_post_mat, yrep_m2_prior_mat,  yrep_m2_post_mat, yrep_m3_prior_mat, yrep_m3_post_mat, yrep_m4_prior_mat, yrep_m4_post_mat, yrep_m5_prior_mat, yrep_m5_post_mat, file="yrep_prior_post_mat.RData")
