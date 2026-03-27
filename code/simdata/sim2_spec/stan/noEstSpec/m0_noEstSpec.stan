data {
  int<lower=1> n; // total number of observations
  array[n] int y; // response variable
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects of age, seifa, =0.5 as per Gelman paper
  int prior_only; // should the likelihood be ignored?
  real<lower=0,upper=1> sens;  // fixed test sensitivity
  real<lower=0,upper=1> spec;  // fixed test specificity, 
}
transformed data {
  
}
parameters {
  real b; // population-level effects  
}
model {
  // likelihood including constants
  if (!prior_only) {
    real p = inv_logit(b);
    real p_sample = p * sens + (1 - p) * (1 - spec);
    y ~ bernoulli(p_sample);
  }
  // priors including constants
  b ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
}
generated quantities {
  real p = inv_logit(b);
  real p_sample = p * sens + (1 - p) * (1 - spec);
  array[n] int y_rep;  // Array of integers for sampled data

  for (i in 1:n) {
    y_rep[i] = bernoulli_rng(p_sample);  // Sampling from Bernoulli distribution
  }
  
  real p_avg1 = p;
  real p_samp_avg1 = p_sample;
}
