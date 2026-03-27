// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #5: STRATA, POSTCODE, SEIFA, SEX, AGEGP

functions {
  /* compute scale parameters of the R2D2 prior
  * Args:
  *   phi: local weight parameters
  *   tau2: global scale parameter
  * Returns:
  *   scale parameter vector of the R2D2 prior
  */
  vector scales_R2D2(vector phi, real tau2) {
    return sqrt(phi * tau2);
  }
  
}

data {
  int<lower=1> Kscales;  // number of local scale parameters
  // data for the R2D2 prior
  real<lower=0> R2D2_mean_R2;  // mean of the R2 prior
  real<lower=0> R2D2_prec_R2;  // precision of the R2 prior
  // concentration vector of the D2 prior
  vector<lower=0>[Kscales] R2D2_cons_D2;
  int<lower=1> K;  // number of population-level effects
  
  int<lower=0> n;                                      // number of tests in sample
  int<lower=0> n_postcode;                             // number of sampled postcodes
  int<lower=0> n_strata;                               // number of sampling strata 
  int<lower=0> n_age;                                  // number of age groups
  int<lower=0> n_seifa;                                // number of SEIFA categories
  array[n] int<lower=0,upper=1> y;                           // Wantai test result, 1=positive, 0=negative
  array[n] int<lower=1,upper=n_postcode> postcode;           // postcode
  array[n] int<lower=1,upper=n_strata> strata;               // sampling strata
  vector<lower=0,upper=1>[n] sex;                      // sex, 1=male, 0=female
  array[n] int<lower=1,upper=n_age> agegp;                   // age group, 1=20-29, 2=30-39,3=40-49, 4=50-59, 5=60-69 
  vector<lower=1,upper=n_seifa>[n] seifa_cont;         // continuous seifa for linear term
  array[n] int<lower=1,upper=n_seifa> seifa;                 // SEIFA decile, 1-10 
  int<lower=0> TP;                                     // number of true positives (for test sens)
  int<lower=0> FN;                                     // number of false negatives (for test sens)
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  int<lower=0> J;                                      // number of poststratification cells
  vector<lower=0>[J] N_pop1;                           // treating sample as The population
  array[J] int<lower=1,upper=n_postcode> postcode_pop;       // poststratification cells, postcode
  array[J] int<lower=1,upper=n_strata> strata_pop;           // poststratification cells, strata
  array[J] int<lower=0,upper=1> sex_pop;                     // poststratification cells, sex
  array[J] int<lower=1,upper=n_age> agegp_pop;               // poststratification cells, agegp
  array[J] int<lower=1,upper=n_seifa> seifa_pop;             // poststratification cells, SEIFA decile
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects
  real<lower=0> sd_seifa_postcode;                     // adjustment to scale parameter for linear term for postcode
  int prior_only;  // should the likelihood be ignored?
}

parameters{
  vector[K] zb;  // unscaled coefficients
  // parameters of the R2D2 prior
  real<lower=0,upper=1> R2D2_R2;
  simplex[Kscales] R2D2_phi;
  
  real<lower=0,upper=1> sens;                          // estimated test sensitivity
  real<lower=0,upper=1> spec;                          // estimated test specificity
  
  real<lower=0> sigma_age;                             // scale parameters for random effects
  real<lower=0> sigma_postcode;
  real<lower=0> sigma_strata;
  real<lower=0> sigma_seifa;   

  vector<multiplier=sigma_age>[n_age] a_age;           // varying intercepts for random effects
  vector<multiplier=sigma_postcode>[n_postcode] a_postcode;  
  vector<multiplier=sigma_strata>[n_strata] a_strata;
  vector<multiplier=sigma_seifa>[n_seifa] a_seifa; 
}

transformed parameters {
  vector[K] b;  // scaled coefficients intercept, coef for sex, seifa_cont
  vector<lower=0>[K] sdb;  // SDs of the coefficients
  real R2D2_tau2;  // global R2D2 scale parameter
  vector<lower=0>[Kscales] scales;  // local R2D2 scale parameters
  
  // compute R2D2 scale parameters
  R2D2_tau2 = R2D2_R2 / (1 - R2D2_R2);
  scales = scales_R2D2(R2D2_phi, R2D2_tau2);
  sdb = scales[(1):(K)];
  b = zb .* sdb;  // scale coefficients
}
model {
  if (!prior_only) {
    vector[n] p = inv_logit(b[1] +
    b[2] * sex +
    b[3] * seifa_cont +
    a_age[agegp] +
    a_postcode[postcode] +
    a_strata[strata] +
    a_seifa[seifa]);
    vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
    y ~ bernoulli(p_sample);
  }
  
  // Prior for the global R2D2 scale parameter
  R2D2_R2 ~ beta(R2D2_mean_R2 * R2D2_prec_R2, (1 - R2D2_mean_R2) * R2D2_prec_R2);
  
  // Priors
  zb ~ normal(0, 1);                  // Standard normal prior for unscaled coefficients
  R2D2_phi ~ dirichlet(R2D2_cons_D2); // Dirichlet prior for R2D2_phi
  
  sens ~ beta(TP+1, FN+1);
  spec ~ beta(TN+1, FP+1);
  
  a_age ~ normal(0, sigma_age);
  a_postcode ~ normal(0, sigma_postcode);
  a_strata ~ normal(0, sigma_strata);
  a_seifa ~ normal(0, sigma_seifa);
  
  b[1] + b[2] * mean(sex) + b[3] * mean(seifa_cont) ~ logistic(0,1);
  b[2] ~ normal(0, coef_prior_scale);
  b[3] ~ normal(0, coef_prior_scale/sd_seifa_postcode);
}

generated quantities {
  vector[n] p = inv_logit(b[1] +
  b[2] * sex +
  b[3] * seifa_cont +
  a_age[agegp] +
  a_postcode[postcode] +
  a_strata[strata] +
  a_seifa[seifa]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  array[n] real y_rep = bernoulli_rng(p_sample);
  real p_mean = mean(p);
  real p_sample_mean = mean(p_sample);
  
  real p_avg1;
  vector[J] p_pop;
  for (j in 1:J) {
    p_pop[j] = inv_logit(b[1] +
    b[2] * sex_pop[j] +
    b[3] * seifa_pop[j] +
    a_age[agegp_pop[j]] +
    a_postcode[postcode_pop[j]] +
    a_strata[strata_pop[j]] +
    a_seifa[seifa_pop[j]]);
  }
  p_avg1 = sum(N_pop1 .* p_pop) / sum(N_pop1);
}
