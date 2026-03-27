
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #1: STRATA

data {
  int<lower=0> n;                                      // number of tests in sample
  int<lower=0> n_strata;                               // number of sampling strata 
  array[n] int<lower=0,upper=1> y;                           // Wantai test result, 1=positive, 0=negative
  array[n] int<lower=1,upper=n_strata> strata;               // sampling strata
  int<lower=0> TP;                                     // number of true positives (for test sens)
  int<lower=0> FN;                                     // number of false negatives (for test sens)
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  int<lower=0> J;                                      // number of poststratification cells (2*5*10=100)
  vector<lower=0>[J] N_pop1;                           // treating sample as The population
  array[J] int<lower=1,upper=n_strata> strata_pop;           // poststratification cells, strata
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects, =0.5 as per Gelman paper
  int prior_only;  // should the likelihood be ignored?
  real<lower=0,upper=1> sens;                          // fixed test sensitivity
  real<lower=0,upper=1> spec;                          // fixed test specificity, 
}
parameters{
  vector[1] b;                                         // intercept
  real<lower=0> sigma_strata;
  vector<multiplier=sigma_strata>[n_strata] a_strata;
}
model{
  if (!prior_only) {
  vector[n] p = inv_logit(b[1] + 
                          a_strata[strata]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  y ~ bernoulli(p_sample);
  }
  a_strata ~ normal(0, sigma_strata);
  b[1] ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
  sigma_strata ~ normal(0, coef_prior_scale);
}
generated quantities{
  vector[n] p = inv_logit(b[1] +
                          a_strata[strata] );
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  array[n] real y_rep = bernoulli_rng(p_sample);
  
  real p_avg1;                                                                         // estimated overall population prevalence - Lifeblood
  vector[J] p_pop;                                                                     // estimated population prevalences in J cells
  for (j in 1:J){
      p_pop[j] = inv_logit(b[1] +
                           a_strata[strata_pop[j]]);
  }
  p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}

