
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #6: Model 5 + two way interaction between SEIFA, SEX, AGEGP

data {
  int<lower=0> n;                                      // number of tests in sample
  int<lower=0> n_postcode;                             // number of sampled postcodes
  int<lower=0> n_strata;                               // number of sampling strata 
  int<lower=0> n_age;                                  // number of age groups
  int<lower=0> n_seifa;                                // number of SEIFA categories
  int<lower=0,upper=1> y[n];                           // Wantai test result, 1=positive, 0=negative
  int<lower=1,upper=n_postcode> postcode[n];           // postcode
  int<lower=1,upper=n_strata> strata[n];               // sampling strata
  vector<lower=0,upper=1>[n] sex;                      // sex, 1=male, 0=female
  int<lower=1,upper=n_age> agegp[n];                   // age group, 1=20-29, 2=30-39,3=40-49, 4=50-59, 5=60-69 
  vector<lower=1,upper=n_seifa>[n] seifa_cont;         // continuous seifa for linear term
  int<lower=1,upper=n_seifa> seifa[n];                 // SEFIA decile, 1-10 
  int<lower=0> TP;                                     // number of true positives (for test sens)
  int<lower=0> FN;                                     // number of false negatives (for test sens)
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  int<lower=0> J;                                      // number of poststratification cells (2*5*10=100)
  vector<lower=0>[J] N_pop1;                           // treating sample as The population
  int<lower=1,upper=n_postcode> postcode_pop[J];       // poststratification cells, postcode
  int<lower=1,upper=n_strata> strata_pop[J];           // poststratification cells, strata
  int<lower=0,upper=1> sex_pop[J];                     // poststratification cells, sex
  int<lower=1,upper=n_age> agegp_pop[J];               // poststratification cells, agegp
  int<lower=1,upper=n_seifa> seifa_pop[J];             // poststratification cells, SEIFA decile
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects of age, seifa, =0.5 as per Gelman paper
  real<lower=0> sd_seifa_postcode;                     // adjustment to scale parameter for linear term for postcode, as per Gelman
  int prior_only;  // should the likelihood be ignored?
  int<lower=0> n_seifa_sex;                            // number of SEIFA categories x sex categories
  int<lower=0> n_seifa_agegp;                          // number of SEIFA categories x agegp categories
  int<lower=0> n_sex_agegp;                            // number of sex categories x agegp 
  int seifa_sex[n];  # for interaction effects
  int seifa_agegp[n];
  int sex_agegp[n];
  int seifa_sex_pop[J];  # interaction effects for poststrat table
  int seifa_agegp_pop[J];  # interaction effects for poststrat table
  int sex_agegp_pop[J];  # interaction effects for poststrat table
  real<lower=0,upper=1> spec;                          // fixed test specificity, 
}
parameters{
  vector[3] b;                                         // intercept, coef for sex, seifa_cont
  real<lower=0> sigma_age;                             // scale parameters for random effects
  real<lower=0> sigma_postcode;
  real<lower=0> sigma_strata;
  real<lower=0> sigma_seifa;   
  real<lower=0> sigma_seifa_sex;   
  real<lower=0> sigma_seifa_agegp;   
  real<lower=0> sigma_sex_agegp;   
  vector<multiplier=sigma_age>[n_age] a_age;           // varying intercepts for random effects
  vector<multiplier=sigma_postcode>[n_postcode] a_postcode;  
  vector<multiplier=sigma_strata>[n_strata] a_strata;
  vector<multiplier=sigma_seifa>[n_seifa] a_seifa; 
  vector<multiplier=sigma_seifa_sex>[n_seifa_sex] a_seifa_sex; 
  vector<multiplier=sigma_seifa_agegp>[n_seifa_agegp] a_seifa_agegp; 
  vector<multiplier=sigma_sex_agegp>[n_sex_agegp] a_sex_agegp; 
  real<lower=0,upper=1> sens;                          // estimated test sensitivity
}
model{
  if (!prior_only) {
  vector[n] p = inv_logit(b[1] + 
                          b[2] * sex + 
                          b[3] * seifa_cont +
                          a_age[agegp] + 
                          a_postcode[postcode] +
                          a_strata[strata] +
                          a_seifa[seifa] +
                          a_seifa_sex[seifa_sex] +
                          a_seifa_agegp[seifa_agegp] +
                          a_sex_agegp[sex_agegp]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  y ~ bernoulli(p_sample);
}
  sens ~ beta(TP+1, FN+1);                             // or use half-normal priors like N+(1,0.05) as per Andrew suggestion?
  a_age ~ normal(0, sigma_age);                        // random effects
  a_postcode ~ normal(0, sigma_postcode);
  a_strata ~ normal(0, sigma_strata);
  a_seifa ~ normal(0, sigma_seifa);
  a_seifa_sex ~ normal(0, sigma_seifa_sex);
  a_seifa_agegp ~ normal(0, sigma_seifa_agegp);
  a_sex_agegp ~ normal(0, sigma_sex_agegp);
  b[1] + b[2] * mean(sex) + b[3] * mean(seifa_cont) ~ logistic(0,1);    // prior on centred intercept, corresponds to U(0,1)
  b[2] ~ normal(0, coef_prior_scale);
  b[3] ~ normal(0, coef_prior_scale/sd_seifa_postcode);
  sigma_age ~ normal(0, coef_prior_scale);                                             // prior for scale parameters for random effects
  sigma_postcode ~ normal(0, coef_prior_scale);
  sigma_strata ~ normal(0, coef_prior_scale);
  sigma_seifa ~ normal(0, coef_prior_scale);  
  sigma_seifa_sex ~ normal(0, coef_prior_scale);  
  sigma_seifa_agegp ~ normal(0, coef_prior_scale);  
  sigma_sex_agegp ~ normal(0, coef_prior_scale);  
}
generated quantities{
  // vector[n] p = inv_logit(b[1] + 
  //                         b[2] * sex + 
  //                         b[3] * seifa_cont +
  //                         a_age[agegp] + 
  //                         a_postcode[postcode] +
  //                         a_strata[strata] +
  //                         a_seifa[seifa] +
  //                         a_seifa_sex[seifa_sex] +
  //                         a_seifa_agegp[seifa_agegp] +
  //                         a_sex_agegp[sex_agegp]);
  // vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  // real y_rep[n] = bernoulli_rng(p_sample);
  real p_avg1;                                                                         // estimated overall population prevalence - Lifeblood
  vector[J] p_pop;                                                                     // estimated population prevalences in J cells
  for (j in 1:J){
      p_pop[j] = inv_logit(b[1] +
                           b[2] * sex_pop[j] +
                           b[3] * seifa_pop[j] +
                           a_age[agegp_pop[j]] +
                           a_postcode[postcode_pop[j]] +
                           a_strata[strata_pop[j]] +
                           a_seifa[seifa_pop[j]] +
                           a_seifa_sex[seifa_sex_pop[j]] +
                           a_seifa_agegp[seifa_agegp_pop[j]] +
                           a_sex_agegp[sex_agegp_pop[j]] );
  }
  p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}

