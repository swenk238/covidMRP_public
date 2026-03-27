// generated with brms 2.19.0
functions {
  
}
data {
  int<lower=1> n; // total number of observations
  array[n] int y; // response variable
  int<lower=1> K; // number of population-level effects

  int<lower=1> n_strata;                        
  array[n] int<lower=1> strata;      
  int<lower=0> FP; // number of false positives (for test spec)
  int<lower=0> TN; // number of true negatives (for test spec)
  int<lower=0> TP; // number of true positives (for test sens)
  int<lower=0> FN; // number of false negatives (for test sens)

  int prior_only; // should the likelihood be ignored?
  
  int<lower=0> J;                           
  vector<lower=0>[J] N_pop1; // treating sample as The population
  matrix[J, K] X_pop1; // poststrat design matrix
  array[J] int<lower=1,upper=n_strata> strata_pop;           // poststratification cells, strata
}
transformed data {
  
}
parameters {
  vector[1] b; // population-level effects  
  real<lower=0,upper=1> spec;  // estimated test specificity, 
  real<lower=0,upper=1> sens;  // estimated test sensitivity

  vector<lower=0>[1] sigma_strata;                             // scale parameters for random effects
  array[1] vector[n_strata] z_strata;
}
transformed parameters {
  vector[n_strata] a_strata;
  a_strata = sigma_strata[1] * z_strata[1];
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[n] p = inv_logit(b[1] + 
                            a_strata[strata]);
    vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
    vector[n] logit_p_samp = logit(p_sample);
    y ~ bernoulli_logit(logit_p_samp);
  }
  // priors including constants
  spec ~ beta(TN+1, FP+1);                           
  sens ~ beta(TP+1, FN+1);    
  z_strata[1] ~ normal(0, 1);                        // random effects
  sigma_strata[1] ~ student_t(3, 0, 2.5) T[0,];
}
generated quantities {
  vector[n] p = inv_logit(b[1] + 
                            a_strata[strata]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  real p_mean = mean(p);
  real p_samp_mean = mean(p_sample);
  
  // poststrat estimate
  vector[J] mu_pop = inv_logit(b[1] + 
  a_strata[strata_pop]);
  vector[J] p_pop = inv_logit(mu_pop + X_pop1*b);
  
  real p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}
