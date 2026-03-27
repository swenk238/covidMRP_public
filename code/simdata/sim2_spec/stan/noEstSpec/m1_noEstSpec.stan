data {
  int<lower=1> n; // total number of observations
  array[n] int y; // response variable
  int<lower=0> n_x1;                        
  array[n] int<lower=1> x1;      
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects of age, seifa, =0.5 as per Gelman paper
  int prior_only; // should the likelihood be ignored?
  real<lower=0,upper=1> sens;  // fixed test sensitivity
  real<lower=0,upper=1> spec; // fixed test specificity, 
}
transformed data {
  
}
parameters {
  vector[1] b; // population-level effects  
  real<lower=0> sigma_x1;                             // scale parameters for random effects
  vector<multiplier=sigma_x1>[n_x1] a_x1;           // varying intercepts for random effects
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[n] p = inv_logit(b[1] + 
                            a_x1[x1]);
    vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
    y ~ bernoulli(p_sample);
  }
  // priors including constants
  a_x1 ~ normal(0, sigma_x1);                        // random effects
  b[1] ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
  sigma_x1 ~ normal(0, coef_prior_scale);  // prior for scale parameters for random effects
}
generated quantities {
  vector[n] p = inv_logit(b[1] + 
                            a_x1[x1]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  real p_avg1 = mean(p);
  real p_samp_avg1 = mean(p_sample);
}
