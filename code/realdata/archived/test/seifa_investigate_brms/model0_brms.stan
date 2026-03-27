
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #0: NO COVARIATES

data {
  int<lower=1> N; // number of tests in sample
  array[N] int Y; // response variable
  int<lower=0> TP;           // number of true positives (for test sens)
  int<lower=0> FN;           // number of false negatives (for test sens)
  int<lower=0> FP;           // number of false positives (for test spec)
  int<lower=0> TN;           // number of true negatives (for test spec)
  int<lower=0> J;            // number of poststratification cells 
  vector<lower=0>[J] N_pop1; // treating sample as The population
  int prior_only; // should the likelihood be ignored?
}
transformed data {
}
parameters{
  real<lower=0,upper=1> sens;                          // estimated test sensitivity
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  real b; // temporary intercept for centered predictors
}
transformed parameters {
  real lprior = 0; // prior contributions to the log posterior
  lprior += student_t_lpdf(b | 3, 0, 2.5);
}
model {
// likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += b;
    vector[N] p = inv_logit(mu);
    vector[N] p_sample = p * sens + (1 - p) * (1 - spec);
    vector[N] logit_p_samp = logit(p_sample);
    Y ~ bernoulli_logit(logit_p_samp);
  }
  
  // priors including constants
  sens ~ beta(TP+1, FN+1);     // or use half-normal priors like N+(1,0.05) as per Andrew suggestion?
  spec ~ beta(TN+1, FP+1);     // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?

  target += lprior;
}
generated quantities {
   vector[N] mu = rep_vector(0.0, N);   // estimated population prevalences in J cells
    for (n in 1:N){
      mu[n] += b;
    } 
    vector[N] p = inv_logit(mu);
    vector[N] p_sample = p * sens + (1 - p) * (1 - spec);

   vector[J] mu_pop = rep_vector(0.0, J);   // estimated population prevalences in J cells
    for (j in 1:J){
      mu_pop[j] += b;
    } 
    vector[J] p_pop = inv_logit(mu_pop);
    real p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}

