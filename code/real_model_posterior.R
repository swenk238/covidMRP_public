#Stan model prior check
library(cmdstanr)
library(tidyverse)
library(posterior)
library(ggridges)
library(colorspace)
library(pROC)
library(here)

####Real world data

## data from Machelek (2020) ####
## reading in population stratification table
melb_popn_tab = read_csv(here('files from Marnie/melb_resident_pop_full.csv'), show_col_types=F) |> 
  rename(seifa = IRSD_Decile) |> 
  mutate(postcode = factor(postcode), 
         agegp = factor(agegp),
         seifa = fct_collapse(factor(seifa),
                              "1" = "2",
                              "2" = c("3", "4"),
                              "3" = c("5", "6"),
                              "4" = c("7", "8"),
                              "5" = c("9", "10"))) %>% 
  filter(!is.na(seifa))

# melb_popn_tab %>% group_by(strata) %>% summarise(sum(N_pop)) # 1-1244776, 2-916034, 3-517722

# reading data from s2 file (sample)
tab = read_csv(here('files from Marnie/s2 table.csv'), show_col_types=F) |> 
  rename(seifa = irsd_decile,
         strata = sampling_group) |> 
  mutate(postcode = factor(postcode), 
         agegp = factor(agegp),
         strata = fct_relevel(factor(strata), c("Low incidence", "Medium incidence", "High incidence")),
         strata = fct_recode(factor(strata), 
                             "1" = "Low incidence", 
                             "2" = "Medium incidence", 
                             "3" = "High incidence"),
         seifa = fct_collapse(factor(seifa),
                              "1" = "2",
                              "2" = c("3", "4"),
                              "3" = c("5", "6"),
                              "4" = c("7", "8"),
                              "5" = c("9", "10")))

dat = tab |> 
  mutate(seifa_cont = as.numeric(seifa),
         seifa_sex = as.factor(as.numeric(seifa)*sex),
         seifa_agegp = as.factor(as.numeric(seifa) * as.numeric(agegp)),
         sex_agegp = as.factor(as.numeric(sex) * as.numeric(agegp)),
         seifa_strata = as.factor(as.numeric(seifa) * as.numeric(strata)),
         sex_strata = as.factor(as.numeric(sex) * as.numeric(strata)),
         agegp_strata = as.factor(as.numeric(agegp) * as.numeric(strata))) |> 
  uncount(total)

data_list_popn = list( n = nrow(dat),                                                    # number of tests in sample
                       n_postcode = dat$postcode |> unique() |> length(),                # number of sampled postcodes
                       n_strata = dat$strata |> unique() |> length(),                    # number of sampling strata 
                       n_age = dat$agegp |> unique() |> length(),                        # number of age groups
                       n_seifa = as.factor(dat$seifa_cont) |> levels() |> length(),      # number of SEIFA categories
                       y = dat$positive,                                                 # Wantai test result, 1=positive, 0=negative
                       postcode = dat$postcode,                                          # postcode
                       strata = as.numeric(as.factor(dat$strata)),                       # converting levels to numbers  # sampling strata
                       sex = dat$sex,                                                    # sex, 1=male, 0=female
                       agegp = dat$agegp,                                                # age group, 1=20-29, 2=30-39,3=40-49, 4=50-59, 5=60-69 
                       seifa_cont = as.numeric(as.factor(dat$seifa_cont)),               # continuous seifa for linear term
                       seifa = dat$seifa,                                                # SEFIA decile, 1-10 
                       TP = 97,                                                          # number of true positives (for test sens) # sensitivity = 97/102 = 0.9509
                       FN = 5,                                                           # number of false negatives (for test sens) 
                       FP = 3,                                                     # number of false positives (for test spec) # specificity = 797/800 = 0.9963
                       TN = 797,                                                       # number of true negatives (for test spec)
                       coef_prior_scale = 0.5,    # also tried alternative 0.1, 1                                       # prior for scale parameter of random effects of age, seifa, =0.5 as per Gelman paper
                       sd_seifa_postcode = sd(dat$seifa_cont),                           # adjustment to scale parameter for linear term for postcode, as per Gelman)
                       J = nrow(melb_popn_tab),                                                    # using sample as the population to poststratify to 
                       N_pop1 = melb_popn_tab$N_pop,
                       postcode_pop = as.numeric(as.factor(melb_popn_tab$postcode)),
                       strata_pop = as.numeric(as.factor(melb_popn_tab$strata)),
                       sex_pop = melb_popn_tab$sex,
                       agegp_pop = melb_popn_tab$agegp,
                       seifa_pop = melb_popn_tab$seifa, 
                       seifa_sex_pop = as.factor(as.numeric(melb_popn_tab$seifa) * (melb_popn_tab$sex)), ## for model 6 ##
                       seifa_agegp_pop = as.factor(as.numeric(melb_popn_tab$seifa) * as.numeric(melb_popn_tab$agegp)),
                       sex_agegp_pop = as.factor(melb_popn_tab$sex * as.numeric(melb_popn_tab$agegp)),
                       seifa_sex = as.factor(as.numeric(dat$seifa)*dat$sex), # for interaction effects
                       seifa_agegp = as.factor(as.numeric(dat$seifa) * as.numeric(dat$agegp)),
                       sex_agegp = as.factor(as.numeric(dat$sex) * as.numeric(dat$agegp)),
                       n_seifa_sex = as.factor(dat$seifa_sex) |> levels() |> length(),      # number of SEIFA x sex categories
                       n_seifa_agegp = as.factor(dat$seifa_agegp) |> levels() |> length(),      # number of SEIFA x agegp categories
                       n_sex_agegp = as.factor(dat$sex_agegp) |> levels() |> length(),      # number of sex x agegp categories
                       seifa_strata_pop = as.factor(as.numeric(melb_popn_tab$seifa) * as.numeric(melb_popn_tab$strata)), ## for model 7 ##
                       sex_strata_pop = as.factor(as.numeric(melb_popn_tab$sex) * as.numeric(melb_popn_tab$strata)),
                       agegp_strata_pop = as.factor(as.numeric(melb_popn_tab$agegp) * as.numeric(melb_popn_tab$strata)),
                       seifa_strata = as.factor(as.numeric(dat$seifa)*as.numeric(dat$strata)), # for interaction effects
                       sex_strata = as.factor(as.numeric(dat$sex) * as.numeric(dat$strata)),
                       agegp_strata = as.factor(as.numeric(dat$agegp) * as.numeric(dat$strata)),
                       n_seifa_strata = as.factor(dat$seifa_strata) |> levels() |> length(),      # number of SEIFA x strata categories
                       n_sex_strata = as.factor(dat$sex_strata) |> levels() |> length(),      # number of sex x strata categories
                       n_agegp_strata = as.factor(dat$ agegp_strata) |> levels() |> length(),      # number of agegp x strata categories
                       prior_only=0, 
                       sens = 0.95,
                       spec = 0.99)

mod0 <- cmdstan_model(here("code/realdata/stan/model0.stan"))

fit0 <- mod0$sample(
  data = data_list_popn,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

mod0_ppc <- posterior::as_draws_df(fit0$draws(c("p_sample","p", 
                                                "sens","spec")))%>%
  dplyr::select(-c('.chain',".draw",".iteration"))%>%
  pivot_longer(everything(),names_to = "parameter", values_to = "posterior_draw") %>%
  mutate(model = "model0")

mod1 <- cmdstan_model(here("code/realdata/stan/model1.stan"))

fit1 <- mod1$sample(
  data = data_list_popn,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

mod1_ppc <- posterior::as_draws_df(fit1$draws(c("p_sample_mean","p_mean", 
                                                "sens","spec")))%>%
  dplyr::select(-c('.chain',".draw",".iteration"))%>%
  pivot_longer(everything(),names_to = "parameter", values_to = "posterior_draw")%>%
  mutate(model = "model1") %>%
  mutate(parameter =  gsub("_mean","",parameter))

mod2 <- cmdstan_model(here("code/realdata/stan/model2.stan"))

fit2 <- mod2$sample(
  data = data_list_popn,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

mod2_ppc <- posterior::as_draws_df(fit2$draws(c("p_sample_mean","p_mean", 
                                                "sens","spec")))%>%
  dplyr::select(-c('.chain',".draw",".iteration"))%>%
  pivot_longer(everything(),names_to = "parameter", values_to = "posterior_draw")%>%
  mutate(model = "model2")%>%
  mutate(parameter =  gsub("_mean","",parameter))

mod3 <- cmdstan_model(here("code/realdata/stan/model3.stan"))

fit3 <- mod3$sample(
  data = data_list_popn,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

mod3_ppc <- posterior::as_draws_df(fit3$draws(c("p_sample_mean","p_mean", 
                                                "sens","spec")))%>%
  dplyr::select(-c('.chain',".draw",".iteration"))%>%
  pivot_longer(everything(),names_to = "parameter", values_to = "posterior_draw")%>%
  mutate(model = "model3")%>%
  mutate(parameter =  gsub("_mean","",parameter))

mod4 <- cmdstan_model(here("code/realdata/stan/model4.stan"))

fit4 <- mod4$sample(
  data = data_list_popn,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

mod4_ppc <- posterior::as_draws_df(fit4$draws(c("p_sample_mean","p_mean", 
                                                "sens","spec")))%>%
  dplyr::select(-c('.chain',".draw",".iteration"))%>%
  pivot_longer(everything(),names_to = "parameter", values_to = "posterior_draw")%>%
  mutate(model = "model4")%>%
  mutate(parameter =  gsub("_mean","",parameter))

mod5 <- cmdstan_model(here("code/realdata/stan/model5.stan"))

fit5 <- mod5$sample(
  data = data_list_popn,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

mod5_ppc <- posterior::as_draws_df(fit5$draws(c("p_sample_mean","p_mean", 
                                                "sens","spec")))%>%
  dplyr::select(-c('.chain',".draw",".iteration"))%>%
  pivot_longer(everything(),names_to = "parameter", values_to = "posterior_draw")%>%
  mutate(model = "model5")%>%
  mutate(parameter =  gsub("_mean","",parameter))

final_df <- rbind(mod0_ppc,mod1_ppc,mod2_ppc,mod3_ppc,mod4_ppc,mod5_ppc) %>%
  group_by(model,parameter) %>%
  summarise(median_posterior = median (posterior_draw),
            low10_posterior = quantile(posterior_draw,.1),
            up90_posterior = quantile(posterior_draw,.9))%>%
  ungroup()%>%
  mutate(model = parse_number(model))

facet_label <- c("p" = "Sample~hat(pi)", "p_sample" = "Sample~Pr(y*'*' == 1)",
                 "sens" = "Sensitivity~hat(delta)", "spec" = "Specificity~hat(gamma)") 

final_df %>%
  ggplot(., aes(y = median_posterior, ymin = low10_posterior, ymax = up90_posterior, x = model))+
  geom_point()+
  geom_hline(data.frame(yintercept = c(0.99625,0.9509804,0.01611,0.01611),parameter = c("spec","sens","p_sample","p")), mapping = aes(yintercept = yintercept), colour = "darkgrey")+
  geom_line()+
  geom_errorbar(width = .1)+
  theme_light()+
  facet_wrap(parameter~., scales = "free_y",  labeller = labeller(parameter =as_labeller(facet_label, default = label_parsed)))+
  xlab("Model") + ylab("Posterior Estimate")

ggsave("code/figures/parameter_values_real_example.png", width =15, height = 10, units = "cm")
