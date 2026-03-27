
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #2: STRATA, POSTCODE

data {
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
  vector[1] b;                                         // intercept
  real<lower=0> sigma_postcode;
  real<lower=0> sigma_strata;
  vector<multiplier=sigma_postcode>[n_postcode] a_postcode;    // varying intercepts for random effects
  vector<multiplier=sigma_strata>[n_strata] a_strata;
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  real<lower=0,upper=1> sens;                          // estimating test sensitivity
}
model{
  if (!prior_only) {
  vector[n] p = inv_logit(b[1] + 
                          a_postcode[postcode] +
                          a_strata[strata]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  y ~ bernoulli(p_sample);
  }
  spec ~ beta(TN+1, FP+1);                             // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?
  sens ~ beta(TP+1, FN+1);                             // or use half-normal priors like N+(1,0.05) as per Andrew suggestion?
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

  real p_avg1;                                                                         // estimated overall population prevalence - Lifeblood
  vector[J] p_pop;                                                                     // estimated population prevalences in J cells
  for (j in 1:J){
      p_pop[j] = inv_logit(b[1] +
                           a_postcode[postcode_pop[j]] +
                           a_strata[strata_pop[j]]);
  }
  p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}

