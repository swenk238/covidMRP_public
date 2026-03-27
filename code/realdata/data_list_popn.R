## data from Machelek (2020) ####
## reading in population stratification table
library(tidyverse)

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
                  n_postcode = melb_popn_tab$postcode |> levels() |> length(),                # number of sampled postcodes
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

