## replicating models 1-5 from Marnie and doing prior and posterior predictive checks
## Jan 2023
library(here)
library(cmdstanr)
source(here("code/realdata/data_list.R"), echo=F)
options(mc.cores = 8)

stan_wd <- "code/realdata/stan/modified_order/"
res_wd <- "results/realdata/original/modified_order/"

# model 4 - no postcode -------------------------------------------------
file_m4_nopostcode <- file.path(here(paste0(stan_wd,"model4_nopostcode.stan")))
mod_m4_nopostcode <- cmdstan_model(file_m4_nopostcode)

## sampling ####
res_m4_nopostcode <- mod_m4_nopostcode$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_m4_nopostcode$summary("p_avg1")

# saving
res_m4_nopostcode$save_object(file=here(paste0(res_wd,"m4_nopostcode.RDS")), compress=T)

# model 3 with sex, seifa, agegp ----------------------------------------------
file_m3_sex_seifa_agegp <- file.path(here(paste0(stan_wd,"model3_sex_seifa_agegp.stan")))
mod_m3_sex_seifa_agegp <- cmdstan_model(file_m3_sex_seifa_agegp)

## sampling ####
res_m3_sex_seifa_agegp <- mod_m3_sex_seifa_agegp$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_m3_sex_seifa_agegp$summary("b[2]")

# saving
res_m3_sex_seifa_agegp$save_object(file=here(paste0(res_wd, "m3_sex_seifa_agegp.RDS")), compress=T)


# model 3 with sex, seifa, strata ----------------------------------------------
file_m3_sex_seifa_strata <- file.path(here(paste0(stan_wd,"model3_sex_seifa_strata.stan")))
mod_m3_sex_seifa_strata <- cmdstan_model(file_m3_sex_seifa_strata)

## sampling ####
res_m3_sex_seifa_strata <- mod_m3_sex_seifa_strata$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_m3_sex_seifa_strata$summary("b[2]")

# saving
res_m3_sex_seifa_strata$save_object(file=here(paste0(res_wd,"m3_sex_seifa_strata.RDS")), compress=T)

# model 2 with sex, seifa (fixed effects only) ----------------------------------------------
file_m2_seifa <- file.path(here(paste0(stan_wd,"model2_seifa.stan")))
mod_m2_seifa <- cmdstan_model(file_m2_seifa)

## sampling ####
res_m2_seifa <- mod_m2_seifa$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_m2_seifa$summary("p_avg1")

# saving
res_m2_seifa$save_object(file=here(paste0(res_wd, "m2_seifa.RDS")), compress=T)


# model 2 with sex, strata ----------------------------------------------
file_m2_strata <- file.path(here(paste0(stan_wd,"model2_strata.stan")))
mod_m2_strata <- cmdstan_model(file_m2_strata)

## sampling ####
res_m2_strata <- mod_m2_strata$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_m2_strata$summary("p_avg1")

# saving
res_m2_strata$save_object(file=here(paste0(res_wd, "m2_strata.RDS")), compress=T)

# model 2 with sex, postcode ----------------------------------------------
file_m2_postcode <- file.path(here(paste0(stan_wd,"model2_postcode.stan")))
mod_m2_postcode <- cmdstan_model(file_m2_postcode)

## sampling ####
res_m2_postcode <- mod_m2_postcode$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_m2_postcode$summary("p_avg1")

# saving
res_m2_postcode$save_object(file=here(paste0(res_wd, "m2_postcode.RDS")), compress=T)


# model 1 sex only ----------------------------------------------
file_m1_sex <- file.path(here(paste0(stan_wd, "model1_sex.stan")))
mod_m1_sex <- cmdstan_model(file_m1_sex)

# posterior predictive
res_m1_sex <- mod_m1_sex$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)
res_m1_sex$summary("p_avg1")

# saving
res_m1_sex$save_object(file=here(paste0(res_wd, "m1_sex.RDS")), compress=T)



# reading in results ------------------------------------------------------



res_list = list()
res_list =  lapply(list.files(here(res_wd), pattern = ".RDS$", full.names = TRUE), readRDS)
  

res_est_full = lapply(res_list, function(x)x$summary(c("p_avg1","b[1]"))[,c('mean', 'q5', 'q95')]) %>% 
  do.call(rbind, .) |> 
  mutate(model = rep(c('1_sex','2_postcode', '2_seifa', '2_strata', '3_sex_seifa_agegp','3_sex_seifa_strata', '4_postcode'),each=2),
         modification = "Modified order", 
         type=rep(c("predicted estimate", "intercept"),7)) 


res_est_full |> 
  as_tibble() |> 
  filter(type == "intercept") |> 
  ggplot(aes(x=model, y=mean, group=1)) + 
  geom_point(size=2) +
  geom_line() +
  geom_text(aes(x='2_strata', y = 0.0165), label="Sample average", col = "gray47") +
  geom_errorbar(position=position_dodge(width=0.5),aes(x=model,ymin=q5, ymax=q95), width=.2, alpha=0.4) +
  scale_color_manual( values=c("#000000",  "#F28E2B", "#E15759",  "#FDE333", "#59A14F", "#59E01F", "#756BB1", "#54278F")  ) +
  geom_hline(yintercept=0.016, lty=2, col="gray47") +
  labs(y = 'MRP estimate',
       x = "Model",
       title = "Original problem: real data (modified order)") +
  ylim(0,0.019) +
  theme_bw() + 
  theme(plot.title = element_text(size=15, colour="#59A14F"))


