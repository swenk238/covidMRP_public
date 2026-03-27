data {
  int<lower=1> n; // total number of observations
  array[n] int y; // response variable
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects of age, seifa, =0.5 as per Gelman paper
  int prior_only; // should the likelihood be ignored?
  real<lower=0,upper=1> sens;  // estimated test sensitivity
}
transformed data {
  
}
parameters {
  real b; // population-level effects  
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
}
model {
  // likelihood including constants
  if (!prior_only) {
    real p = inv_logit(b);
    real p_sample = p * sens + (1 - p) * (1 - spec);
    y ~ bernoulli(p_sample);
  }
  // priors including constants
  spec ~ beta(TN+1, FP+1);                             // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?
  b ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
}
generated quantities {
  real p = inv_logit(b);
  real p_sample = p * sens + (1 - p) * (1 - spec);
  real p_avg1 = p;
  real p_samp_avg1 = p_sample;
}
