## datalist
library(here)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(arm)

## simulated data ####
## May 2023
## attempt to simulate data 
## simulate data so we can see if we recover the values
set.seed(1234)

gen_data_func_prior = function(mu_a, sigma_a){
  n_j = 40
  J = 100
  N = n_j*J
  a_sens = 98 # TP + 1
  b_sens = 6 # FN + 1 
  a_spec = 798 # TN + 1
  b_spec = 4 # FP + 1 
  sens = rbeta(1, a_sens, b_sens)
  spec = rbeta(1, a_spec, b_spec) 
  a_j = rnorm(J, mu_a, sigma_a)
  p_i = invlogit(rep(a_j, n_j))
  group = rep(1:J, n_j) # create n_j=40 samples in each Jth group
  p_obs <- sens*p_i + (1 - p_i)*(1-spec)
  y <- rbinom(N, 1, p_obs)
  data_1 <- list(J=J, group=group, N=N, y=y, p_i=p_i, sens=sens, spec=spec, a_sens=a_sens, b_sens=b_sens, a_spec=a_spec, b_spec=b_spec, mu_a=mu_a, sigma_a=sigma_a)
  return(data_1)
}

gen_data_func_cov = function(p_sim){
  n_j = 40 # number in each postcode
  J = 100 # total number in poststrat table (same as postcode)
  group = rep(1:J, rep(n_j,J)) # create n_j=40 samples in each Jth group
  N = length(group)
  n_X1 = 3
  n_X2 = 102
  n_X3 = 5
  n_X4 = 2
  n_X5 = 5
  X1 = rep(1:n_X1, length.out=n_j*J)
  X2 = rep(1:n_X2, length.out=n_j*J)
  X3 = rep(1:n_X3, length.out=n_j*J)
  X4 = rep(1:n_X4, length.out=n_j*J)
  X5 = rep(1:n_X5, length.out=n_j*J)
  p <- invlogit(rnorm(n_j*J, logit(p_sim), 0.5))
  p <- p*p_sim/mean(p) # to force p to have mean = p_sim
  sens <- 0.95
  spec <- 0.99
  p_obs <- sens*p[group] + (1 - p[group])*(1-spec)
  y <- rbinom(N, 1, p_obs)
  data_1 <- list(J=J, group=group, N=N, y=y, p=p, n_X1=n_X1, n_X2=n_X2, n_X3=n_X3, n_X4=n_X4, n_X5=n_X5, X1=X1, X2=X2, X3=X3, X4=X4, X5=X5)
  return(data_1)
}

gen_data_func = function(p_sim){
  n_j = 40 # number in each postcode
  J = 100 # total number in poststrat table (same as postcode)
  group = rep(1:J, rep(n_j,J)) # create n_j=40 samples in each Jth group
  N = length(group)
  p <- invlogit(rnorm(J, logit(p_sim), 0.5))
  p <- p*p_sim/mean(p) # to force p to have mean = p_sim
  sens <- 0.95
  spec <- 0.99
  p_obs <- sens*p[group] + (1 - p[group])*(1-spec)
  y <- rbinom(N, 1, p_obs)
  data_1 <- list(J=J, group=group, N=N, y=y, p=p, p_obs=p_obs)
  return(data_1)
}

## data from Machelek (2020) ####
# reading data from s2 file
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

data_list = list( n = nrow(dat),                                                    # number of tests in sample
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
                  J = nrow(tab),                                                    # using sample as the population to poststratify to 
                  N_pop1 = tab$total,
                  postcode_pop = as.numeric(as.factor(tab$postcode)),
                  strata_pop = as.numeric(as.factor(tab$strata)),
                  sex_pop = tab$sex,
                  agegp_pop = tab$agegp,
                  seifa_pop = tab$seifa, 
                  seifa_sex_pop = as.factor(as.numeric(tab$seifa) * (tab$sex)), ## for model 6 ##
                  seifa_agegp_pop = as.factor(as.numeric(tab$seifa) * as.numeric(tab$agegp)),
                  sex_agegp_pop = as.factor(tab$sex * as.numeric(tab$agegp)),
                  seifa_sex = as.factor(as.numeric(dat$seifa)*dat$sex), # for interaction effects
                  seifa_agegp = as.factor(as.numeric(dat$seifa) * as.numeric(dat$agegp)),
                  sex_agegp = as.factor(as.numeric(dat$sex) * as.numeric(dat$agegp)),
                  n_seifa_sex = as.factor(dat$seifa_sex) |> levels() |> length(),      # number of SEIFA x sex categories
                  n_seifa_agegp = as.factor(dat$seifa_agegp) |> levels() |> length(),      # number of SEIFA x agegp categories
                  n_sex_agegp = as.factor(dat$sex_agegp) |> levels() |> length(),      # number of sex x agegp categories
                  seifa_strata_pop = as.factor(as.numeric(tab$seifa) * as.numeric(tab$strata)), ## for model 7 ##
                  sex_strata_pop = as.factor(as.numeric(tab$sex) * as.numeric(tab$strata)),
                  agegp_strata_pop = as.factor(as.numeric(tab$agegp) * as.numeric(tab$strata)),
                  seifa_strata = as.factor(as.numeric(dat$seifa)*as.numeric(dat$strata)), # for interaction effects
                  sex_strata = as.factor(as.numeric(dat$sex) * as.numeric(dat$strata)),
                  agegp_strata = as.factor(as.numeric(dat$agegp) * as.numeric(dat$strata)),
                  n_seifa_strata = as.factor(dat$seifa_strata) |> levels() |> length(),      # number of SEIFA x strata categories
                  n_sex_strata = as.factor(dat$sex_strata) |> levels() |> length(),      # number of sex x strata categories
                  n_agegp_strata = as.factor(dat$ agegp_strata) |> levels() |> length(),      # number of agegp x strata categories
                  prior_only=0, 
                  sens = 0.95,
                  spec = 0.99)


set.seed(1234)
dat_subset <- dat[sample(nrow(dat),2000),]

data_list_subset = list( n = nrow(dat_subset),                                                    # number of tests in sample
                  n_postcode = dat_subset$postcode |> levels() |> length(),                # number of sampled postcodes
                  n_strata = dat_subset$strata |> unique() |> length(),                    # number of sampling strata 
                  n_age = dat_subset$agegp |> unique() |> length(),                        # number of age groups
                  n_seifa = as.factor(dat_subset$seifa_cont) |> levels() |> length(),      # number of SEIFA categories
                  y = dat_subset$positive,                                                 # Wantai test result, 1=positive, 0=negative
                  postcode = dat_subset$postcode,                                          # postcode
                  strata = as.numeric(as.factor(dat_subset$strata)),                       # converting levels to numbers  # sampling strata
                  sex = dat_subset$sex,                                                    # sex, 1=male, 0=female
                  agegp = dat_subset$agegp,                                                # age group, 1=20-29, 2=30-39,3=40-49, 4=50-59, 5=60-69 
                  seifa_cont = as.numeric(as.factor(dat_subset$seifa_cont)),               # continuous seifa for linear term
                  seifa = dat_subset$seifa,                                                # SEFIA decile, 1-10 
                  TP = 97,                                                          # number of true positives (for test sens)
                  FN = 5,                                                           # number of false negatives (for test sens)
                  FP = 0, # original = 3                                                          # number of false positives (for test spec)
                  TN = 10,   #  original = 797                                                      # number of true negatives (for test spec)
                  coef_prior_scale = 0.5,                                           # prior for scale parameter of random effects of age, seifa, =0.5 as per Gelman paper
                  sd_seifa_postcode = sd(dat_subset$seifa_cont),                           # adjustment to scale parameter for linear term for postcode, as per Gelman)
                  J = nrow(tab),                                                    # using sample as the population to poststratify to 
                  N_pop1 = tab$total,
                  postcode_pop = as.numeric(as.factor(tab$postcode)),
                  strata_pop = as.numeric(as.factor(tab$strata)),
                  sex_pop = tab$sex,
                  agegp_pop = tab$agegp,
                  seifa_pop = tab$seifa, 
                  seifa_sex_pop = as.factor(as.numeric(tab$seifa) * (tab$sex)), ## for model 6 ##
                  seifa_agegp_pop = as.factor(as.numeric(tab$seifa) * as.numeric(tab$agegp)),
                  sex_agegp_pop = as.factor(tab$sex * as.numeric(tab$agegp)),
                  seifa_sex = as.factor(as.numeric(dat_subset$seifa)*dat_subset$sex), # for interaction effects
                  seifa_agegp = as.factor(as.numeric(dat_subset$seifa) * as.numeric(dat_subset$agegp)),
                  sex_agegp = as.factor(as.numeric(dat_subset$sex) * as.numeric(dat_subset$agegp)),
                  n_seifa_sex = as.factor(dat_subset$seifa_sex) |> levels() |> length(),      # number of SEIFA x sex categories
                  n_seifa_agegp = as.factor(dat_subset$seifa_agegp) |> levels() |> length(),      # number of SEIFA x agegp categories
                  n_sex_agegp = as.factor(dat_subset$sex_agegp) |> levels() |> length(),      # number of sex x agegp categories
                  seifa_strata_pop = as.factor(as.numeric(tab$seifa) * as.numeric(tab$strata)), ## for model 7 ##
                  sex_strata_pop = as.factor(as.numeric(tab$sex) * as.numeric(tab$strata)),
                  agegp_strata_pop = as.factor(as.numeric(tab$agegp) * as.numeric(tab$strata)),
                  seifa_strata = as.factor(as.numeric(dat_subset$seifa)*as.numeric(dat_subset$strata)), # for interaction effects
                  sex_strata = as.factor(as.numeric(dat_subset$sex) * as.numeric(dat_subset$strata)),
                  agegp_strata = as.factor(as.numeric(dat_subset$agegp) * as.numeric(dat_subset$strata)),
                  n_seifa_strata = as.factor(dat_subset$seifa_strata) |> levels() |> length(),      # number of SEIFA x strata categories
                  n_sex_strata = as.factor(dat_subset$sex_strata) |> levels() |> length(),      # number of sex x strata categories
                  n_agegp_strata = as.factor(dat_subset$ agegp_strata) |> levels() |> length(),      # number of agegp x strata categories
                  prior_only=0)
