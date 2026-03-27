
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #2: STRATA, POSTCODE

// adding R2D2M2 prrior
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
  array[n] int<lower=0,upper=1> y;                           // Wantai test result, 1=positive, 0=negative
  array[n] int<lower=1,upper=n_postcode> postcode;           // postcode
  array[n] int<lower=1,upper=n_strata> strata;               // sampling strata
  int<lower=0> TP;                                     // number of true positives (for test sens)
  int<lower=0> FN;                                     // number of false negatives (for test sens)
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  int<lower=0> J;                                      // number of poststratification cells (2*5*10=100)
  vector<lower=0>[J] N_pop1;                           // treating sample as The population
  array[J] int<lower=1,upper=n_postcode> postcode_pop;       // poststratification cells, postcode
  array[J] int<lower=1,upper=n_strata> strata_pop;           // poststratification cells, strata
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects, =0.5 as per Gelman paper
  int prior_only;  // should the likelihood be ignored?
}
parameters{
  vector[1] zb;  // unscaled coefficients
  // parameters of the R2D2 prior
  real<lower=0,upper=1> R2D2_R2;
  simplex[Kscales] R2D2_phi;
  
  real<lower=0,upper=1> sens;                          // estimated test sensitivity
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  real<lower=0> sigma_postcode;
  real<lower=0> sigma_strata;
  vector<multiplier=sigma_postcode>[n_postcode] a_postcode;    // varying intercepts for random effects
  vector<multiplier=sigma_strata>[n_strata] a_strata;
}
transformed parameters {
  vector[1] b;  // scaled coefficients
  vector<lower=0>[1] sdb;  // SDs of the coefficients
  real R2D2_tau2;  // global R2D2 scale parameter
  vector<lower=0>[Kscales] scales;  // local R2D2 scale parameters
  
  // compute R2D2 scale parameters
  R2D2_tau2 = R2D2_R2 / (1 - R2D2_R2);
  scales = scales_R2D2(R2D2_phi, R2D2_tau2);
  sdb = scales[(1):(1)];
  b = zb .* sdb;  // scale coefficients
}
model{
  if (!prior_only) {
    vector[n] p = inv_logit(b[1] +
    a_postcode[postcode] +
    a_strata[strata]);
    vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
    y ~ bernoulli(p_sample);
  }
  
  // Prior for the global R2D2 scale parameter
  R2D2_R2 ~ beta(R2D2_mean_R2 * R2D2_prec_R2, (1 - R2D2_mean_R2) * R2D2_prec_R2);
  
  // Priors
  zb ~ normal(0, 1);                  // Standard normal prior for unscaled coefficients
  R2D2_phi ~ dirichlet(R2D2_cons_D2); // Dirichlet prior for R2D2_phi
  
  sens ~ beta(TP+1, FN+1);                             // or use half-normal priors like N+(1,0.05) as per Andrew suggestion?
  spec ~ beta(TN+1, FP+1);                             // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?
  a_postcode ~ normal(0, sigma_postcode);             // random effects
  a_strata ~ normal(0, sigma_strata);
  b[1] ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
  sigma_postcode ~ normal(0, coef_prior_scale);
  sigma_strata ~ normal(0, coef_prior_scale);
}
generated quantities{
  vector[n] p = inv_logit(b[1] +
  a_postcode[postcode] +
  a_strata[strata] );
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  array[n] real y_rep = bernoulli_rng(p_sample);
  real p_mean = mean(p);
  real p_sample_mean = mean(p_sample);
  
  real p_avg1;                                                                         // estimated overall population prevalence - Lifeblood
  vector[J] p_pop;                                                                     // estimated population prevalences in J cells
  for (j in 1:J){
    p_pop[j] = inv_logit(b[1] +
    a_postcode[postcode_pop[j]] +
    a_strata[strata_pop[j]]);
  }
  p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}

