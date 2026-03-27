// Model 3 - dropping age, sex
functions {
  
}
data {
  int<lower=1> n; // total number of observations
  array[n] int y; // response variable
  int<lower=1> K; // number of population-level effects
  matrix[n, K] X; // population-level design matrix
  
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  int<lower=0> TP;                                     // number of true positives (for test sens)
  int<lower=0> FN;                                     // number of false negatives (for test sens)
  
  int prior_only; // should the likelihood be ignored?
  
  // data for group-level effects of ID 2
  int<lower=1> n_postcode; // number of grouping levels
  array[n] int<lower=1> postcode; // grouping indicator per observation
  
  // data for group-level effects of ID 3
  int<lower=1> n_strata; // number of grouping levels
  array[n] int<lower=1> strata; // grouping indicator per observation
  
  // data for group-level effects of ID 4
  int<lower=1> n_seifa; // number of grouping levels
  array[n] int<lower=1> seifa; // grouping indicator per observation
  
  int<lower=0> J;                           
  vector<lower=0>[J] N_pop1; // treating sample as The population
  matrix[J, K] X_pop_seifa; // poststrat design matrix
  
  array[J] int<lower=1,upper=n_postcode> postcode_pop;       // poststratification cells, postcode
  array[J] int<lower=1,upper=n_strata> strata_pop;           // poststratification cells, strata
  array[J] int<lower=1,upper=n_age> agegp_pop;               // poststratification cells, agegp
  array[J] int<lower=1,upper=n_seifa> seifa_pop;             // poststratification cells, SEIFA decile
}
transformed data {
  
}
parameters {
  vector[K] b; // population-level effects  
  real<lower=0,upper=1> spec;  // estimated test specificity
  real<lower=0,upper=1> sens;  // estimated test sensitivity
  
  array[1] vector[n_postcode] z_postcode; // standardized group-level effects
  array[1] vector[n_strata] z_strata; // standardized group-level effects
  array[1] vector[n_seifa] z_seifa; // standardized group-level effects
  
  vector<lower=0>[1] sigma_postcode; // group-level standard deviations
  vector<lower=0>[1] sigma_strata; // group-level standard deviations
  vector<lower=0>[1] sigma_seifa; // group-level standard deviations
}
transformed parameters {
  vector[n_postcode] a_postcode; // actual group-level effects
  vector[n_strata] a_strata; // actual group-level effects
  vector[n_seifa] a_seifa; // actual group-level effects
  a_postcode = sigma_postcode[1] * z_postcode[1];
  a_strata = sigma_strata[1] * z_strata[1];
  a_seifa = sigma_seifa[1] * z_seifa[1];
}
model {
  // likelihood including constants
  if (!prior_only) {
    
    vector[n] mu = rep_vector(0.0, n);
    mu = inv_logit(b[1]  +
    a_postcode[postcode] +
    a_strata[strata] +
    a_seifa[seifa]);
    vector[n] p = inv_logit(mu + X*b);
    vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
    // vector[n] logit_p_samp = logit(p_sample);
    y ~ bernoulli(p_sample);
  }
  // priors including constants
  spec ~ beta(TN+1, FP+1); 
  sens ~ beta(TP+1, FN+1);    
  b[1] ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
  z_postcode[1] ~ normal(0, 1);                        // random effects
  z_strata[1] ~ normal(0, 1);                        // random effects
  z_seifa[1] ~ normal(0, 1);                        // random effects
  
  sigma_postcode ~ student_t(3, 0, 2.5);
  sigma_strata ~ student_t(3, 0, 2.5);
  sigma_seifa ~ student_t(3, 0, 2.5);
}
generated quantities {
  vector[n] mu = inv_logit(b[1]  +
  a_postcode[postcode] +
  a_strata[strata] +
  a_seifa[seifa]);
  vector[n] p = inv_logit(mu + X*b);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  
  real p_mean = mean(p);
  real p_samp_mean = mean(p_sample);
  
  // poststrat estimate
  vector[J] mu_pop = inv_logit(b[1] + 
  a_postcode[postcode_pop] + 
  a_strata[strata_pop] + 
  a_seifa[seifa_pop]);
  vector[J] p_pop = inv_logit(mu_pop + X_pop_sex_seifa*b);
  
  real p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}
