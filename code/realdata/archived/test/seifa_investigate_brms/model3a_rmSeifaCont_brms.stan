
// Stan model for COVID-19 Melbourne serosurvey
// Bayesian estimation of COVID-19 seroprevalence using multilevel regression and poststratification (MRP)
// Based on Gelman & Carpenter 2020 paper:
// "Bayesian analysis of tests with unknown sensitivity and specificity"

// Test model #3: STRATA, POSTCODE, SEIFA

data {
  int<lower=1> N; // number of tests in sample
  array[N] int Y; // response variable
  int<lower=1> K; // number of population-level effects
  // data for group-level effects of ID 1
  int<lower=1> N_1; // number of grouping levels
  int<lower=1> M_1; // number of coefficients per level
  array[N] int<lower=1> J_1; // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2; // number of grouping levels
  int<lower=1> M_2; // number of coefficients per level
  array[N] int<lower=1> J_2; // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  int<lower=0> TP;           // number of true positives (for test sens)
  int<lower=0> FN;           // number of false negatives (for test sens)
  int<lower=0> FP;           // number of false positives (for test spec)
  int<lower=0> TN;           // number of true negatives (for test spec)
  int prior_only;  // should the likelihood be ignored?
  int<lower=0> J;            // number of poststratification cells (2*5*10=100)
  int<lower=0> n_postcode;   // number of sampled postcodes
  int<lower=0> n_strata;     // number of sampling strata 
  vector<lower=0>[J] N_pop1; // treating sample as The population
  array[J] int<lower=1,upper=n_postcode> postcode_pop;       // poststratification cells, postcode
  array[J] int<lower=1,upper=n_strata> strata_pop;           // poststratification cells, strata
}
transformed data {
}
parameters{
  real<lower=0,upper=1> sens;                          // estimated test sensitivity
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  vector[K] b; // regression coefficients
  vector<lower=0>[M_1] sd_1; // group-level standard deviations
  array[M_1] vector[N_1] z_1; // standardized group-level effects
  vector<lower=0>[M_2] sd_2; // group-level standard deviations
  array[M_2] vector[N_2] z_2; // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1; // actual group-level effects
  vector[N_2] r_2_1; // actual group-level effects
  real lprior = 0; // prior contributions to the log posterior
  r_1_1 = sd_1[1] * z_1[1];
  r_2_1 = sd_2[1] * z_2[1];
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
   // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    for (n in 1 : N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
    }
    target += bernoulli_logit_lpmf(Y | mu);
  }

  // priors including constants
  sens ~ beta(TP+1, FN+1);     // or use half-normal priors like N+(1,0.05) as per Andrew suggestion?
  spec ~ beta(TN+1, FP+1);     // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?

  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {
    vector[N] mu = rep_vector(0.0, N);
    for (n in 1 : N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
    }
    vector[N] p = inv_logit(mu);
    vector[N] p_sample = p * sens + (1 - p) * (1 - spec);

    real p_avg1;                                                                         // estimated overall population prevalence - Lifeblood
    vector[J] mu_pop = rep_vector(0.0, J);
    vector[J] p_pop;                                                                     // estimated population prevalences in J cells
    for (j in 1:J){
      mu_pop[j] += r_1_1[postcode_pop[j]] * Z_1_1[j] + r_2_1[strata_pop[j]] * Z_2_1[j];
    }
    p_pop = inv_logit(mu_pop);

    p_avg1 = sum(N_pop1 .* p_pop)/sum(N_pop1);
}

