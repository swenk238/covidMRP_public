
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #0: intercept-only

data {
  int<lower=0> n;                                      // number of tests in sample
  array[n] int<lower=0,upper=1> y;                     // Wantai test result, 1=positive, 0=negative
  int<lower=0> TP;                                     // number of true positives (for test sens)
  int<lower=0> FN;                                     // number of false negatives (for test sens)
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  int<lower=0> J;                                      // number of poststratification cells (2*5*10=100)
  vector<lower=0>[J] N_pop1;                           // treating sample as The population
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects, =0.5 as per Gelman paper
}
parameters{
  real<lower=0,upper=1> sens;                          // estimated test sensitivity
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  vector[1] b;                                         // intercept
}
model{
  real p = inv_logit(b[1]);
  real p_sample = p + (1 - spec); // p * sens + (1 - p) * (1 - spec);
  y ~ bernoulli(p_sample);

  sens ~ beta(TP+1, FN+1);                             // or use half-normal priors like N+(1,0.05) as per Andrew suggestion?
  spec ~ beta(TN+1, FP+1);                             // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?
  b[1] ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
}
generated quantities{
  real p = inv_logit( b[1] );
  real p_sample = p * sens + (1 - p) * (1 - spec);
  // int<lower=0,upper=1> y_rep = bernoulli_rng(p_sample);
  // real p_mean = p;
  // real p_sample_mean = p_sample;
  
  real p_avg1;                                                                         // estimated overall population prevalence - Lifeblood
  vector[J] p_pop;                                                                     // estimated population prevalences in J cells
  for (j in 1:J){
    p_pop[j] = inv_logit(b[1]);
  }
  p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
  
  array[n] int y_pred;
  real joint_log_lik = 0;
  for (i in 1:n) {
    y_pred[i] = bernoulli_rng(p_sample);
    joint_log_lik += bernoulli_lpmf(y[i]|p_sample);
    joint_log_lik += bernoulli_lpmf(y_pred[i]|p_sample);
  }
}

