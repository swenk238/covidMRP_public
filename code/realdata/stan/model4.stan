
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #4: STRATA, POSTCODE, SEIFA, SEX

data {
  int<lower=0> n;                                      // number of tests in sample
  int<lower=0> n_postcode;                             // number of sampled postcodes
  int<lower=0> n_strata;                               // number of sampling strata 
  int<lower=0> n_seifa;                                // number of SEIFA categories
  array[n] int<lower=0,upper=1> y;                           // Wantai test result, 1=positive, 0=negative
  array[n] int<lower=1,upper=n_postcode> postcode;           // postcode
  array[n] int<lower=1,upper=n_strata> strata;               // sampling strata
  vector<lower=0,upper=1>[n] sex;                      // sex, 1=male, 0=female
  vector<lower=1,upper=n_seifa>[n] seifa_cont;         // continuous seifa for linear term
  array[n] int<lower=1,upper=n_seifa> seifa;                 // SEFIA decile, 1-10 
  int<lower=0> TP;                                     // number of true positives (for test sens)
  int<lower=0> FN;                                     // number of false negatives (for test sens)
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  int<lower=0> J;                                      // number of poststratification cells
  vector<lower=0>[J] N_pop1;                           // treating sample as The population
  // array[J] int<lower=1,upper=n_postcode> postcode_pop;       // poststratification cells, postcode
  // array[J] int<lower=1,upper=n_strata> strata_pop;           // poststratification cells, strata
  // array[J] int<lower=0,upper=1> sex_pop;                     // poststratification cells, sex
  // array[J] int<lower=1,upper=n_seifa> seifa_pop;             // poststratification cells, SEIFA decile
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects of seifa, =0.5 as per Gelman paper
  real<lower=0> sd_seifa_postcode;                     // adjustment to scale parameter for linear term for postcode, as per Gelman
  int prior_only;  // should the likelihood be ignored?
}
parameters{
  real<lower=0,upper=1> sens;                          // estimated test sensitivity
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  vector[3] b;                                         // intercept, coef for sex, seifa_cont
  real<lower=0> sigma_postcode;
  real<lower=0> sigma_strata;
  real<lower=0> sigma_seifa;   
  vector<multiplier=sigma_postcode>[n_postcode] a_postcode;    // varying intercepts for random effects
  vector<multiplier=sigma_strata>[n_strata] a_strata;
  vector<multiplier=sigma_seifa>[n_seifa] a_seifa; 
}
model{
  if (!prior_only) {
  vector[n] p = inv_logit(b[1] + 
                          b[2] * sex + 
                          b[3] * seifa_cont +
                          a_postcode[postcode] +
                          a_strata[strata] +
                          a_seifa[seifa]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  y ~ bernoulli(p_sample);
  }
  sens ~ beta(TP+1, FN+1);                             // or use half-normal priors like N+(1,0.05) as per Andrew suggestion?
  spec ~ beta(TN+1, FP+1);                             // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?
  a_postcode ~ normal(0, sigma_postcode);             // random effects
  a_strata ~ normal(0, sigma_strata);
  a_seifa ~ normal(0, sigma_seifa);
  b[1] + b[2] * mean(sex) + b[3] * mean(seifa_cont) ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
  b[2] ~ normal(0, coef_prior_scale);
  b[3] ~ normal(0, coef_prior_scale/sd_seifa_postcode);
  sigma_postcode ~ normal(0, coef_prior_scale);
  sigma_strata ~ normal(0, coef_prior_scale);
  sigma_seifa ~ normal(0, coef_prior_scale);  
}
generated quantities{
  vector[n] p = inv_logit(b[1] + 
                          b[2] * sex + 
                          b[3] * seifa_cont +
                          a_postcode[postcode] +
                          a_strata[strata] +
                          a_seifa[seifa]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  array[n] real y_rep = bernoulli_rng(p_sample);
  real p_mean = mean(p);
  real p_sample_mean = mean(p_sample);

  // real p_avg1;                                                                         // estimated overall population prevalence - Lifeblood
  // vector[J] p_pop;                                                                     // estimated population prevalences in J cells
  // for (j in 1:J){
  //     p_pop[j] = inv_logit(b[1] +
  //                          b[2] * sex_pop[j] +
  //                          b[3] * seifa_pop[j] +
  //                          a_postcode[postcode_pop[j]] +
  //                          a_strata[strata_pop[j]] +
  //                          a_seifa[seifa_pop[j]]);
  // }
  // p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}

