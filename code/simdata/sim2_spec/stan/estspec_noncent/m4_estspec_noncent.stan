data {
  int<lower=1> n; // total number of observations
  array[n] int y; // response variable
  int<lower=0> n_x1;                        
  int<lower=0> n_x2;                           
  int<lower=0> n_x3;                           
  int<lower=0> n_x4;                              
  array[n] int<lower=1> x1;      
  array[n] int<lower=1> x2;           
  array[n] int<lower=1> x3;           
  array[n] int<lower=1> x4;           
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects of age, seifa, =0.5 as per Gelman paper
  int prior_only; // should the likelihood be ignored?
  real<lower=0,upper=1> sens;  // estimated test sensitivity
}
transformed data {
  
}
parameters {
  vector[1] b; // population-level effects  
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  real<lower=0> sigma_x1;                             // scale parameters for random effects
  real<lower=0> sigma_x2;                             // scale parameters for random effects
  real<lower=0> sigma_x3;                             // scale parameters for random effects
  real<lower=0> sigma_x4;                             // scale parameters for random effects
  vector<multiplier=sigma_x1>[n_x1] a_x1;           // varying intercepts for random effects
  vector<multiplier=sigma_x2>[n_x2] a_x2;           // varying intercepts for random effects
  vector<multiplier=sigma_x3>[n_x3] a_x3;           // varying intercepts for random effects
  vector<multiplier=sigma_x4>[n_x4] a_x4;           // varying intercepts for random effects
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[n] p = inv_logit(b[1] + 
                            a_x1[x1] + 
                            a_x2[x2] + 
                            a_x3[x3] + 
                            a_x4[x4]);
    vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
    y ~ bernoulli(p_sample);
  }
  // priors including constants
  spec ~ beta(TN+1, FP+1);                             // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?
  a_x1 ~ normal(0, sigma_x1);                        // random effects
  a_x2 ~ normal(0, sigma_x2);                        // random effects
  a_x3 ~ normal(0, sigma_x3);                        // random effects
  a_x4 ~ normal(0, sigma_x4);                        // random effects
  b[1] ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
  sigma_x1 ~ normal(0, coef_prior_scale);  // prior for scale parameters for random effects
  sigma_x2 ~ normal(0, coef_prior_scale);           
  sigma_x3 ~ normal(0, coef_prior_scale);           
  sigma_x4 ~ normal(0, coef_prior_scale);           
}
generated quantities {
  vector[n] p = inv_logit(b[1] + 
                            a_x1[x1] + 
                            a_x2[x2] + 
                            a_x3[x3] + 
                            a_x4[x4]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  real p_avg1 = mean(p);
  real p_samp_avg1 = mean(p_sample);
}
