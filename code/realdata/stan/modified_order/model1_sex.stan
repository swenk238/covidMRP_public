
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #4: STRATA, POSTCODE, SEIFA, SEX

data {
  int<lower=0> n;                                      // number of tests in sample
  array[n] int<lower=0,upper=1> y;                           // Wantai test result, 1=positive, 0=negative
  vector<lower=0,upper=1>[n] sex;                      // sex, 1=male, 0=female
  int<lower=0> TP;                                     // number of true positives (for test sens)
  int<lower=0> FN;                                     // number of false negatives (for test sens)
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  int<lower=0> J;                                      // number of poststratification cells
  vector<lower=0>[J] N_pop1;                           // treating sample as The population
  array[J] int<lower=0,upper=1> sex_pop;                     // poststratification cells, sex
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects of seifa, =0.5 as per Gelman paper
  int prior_only;  // should the likelihood be ignored?
}
parameters{
  real<lower=0,upper=1> sens;                          // estimated test sensitivity
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  vector[2] b;                                         // intercept, coef for sex
}
model{
  if (!prior_only) {
  vector[n] p = inv_logit(b[1] + 
                          b[2] * sex);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  y ~ bernoulli(p_sample);
  }
  sens ~ beta(TP+1, FN+1);                             // or use half-normal priors like N+(1,0.05) as per Andrew suggestion?
  spec ~ beta(TN+1, FP+1);                             // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?
  b[1] + b[2] * mean(sex) ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
  b[2] ~ normal(0, coef_prior_scale);
}
generated quantities{
  vector[n] p = inv_logit(b[1] + 
                          b[2] * sex);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  array[n] real y_rep = bernoulli_rng(p_sample);
  real p_mean = mean(p);
  real p_sample_mean = mean(p_sample);
  
  real p_avg1;                                                                         // estimated overall population prevalence - Lifeblood
  vector[J] p_pop;                                                                     // estimated population prevalences in J cells
  for (j in 1:J){
      p_pop[j] = inv_logit(b[1] +
                           b[2] * sex_pop[j]);
  }
  p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}

