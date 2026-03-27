// adapted from exp4/estSpec/non-centered/m3_estSpec_noncent.stan
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
  
  int<lower=1> n; // total number of observations
  array[n] int y; // response variable
  int<lower=0> n_x1;
  int<lower=0> n_x2;
  int<lower=0> n_x3;      
  array[n] int<lower=1> x1;      
  array[n] int<lower=1> x2;           
  array[n] int<lower=1> x3;
  int<lower=0> FP;                                     // number of false positives (for test spec)
  int<lower=0> TN;                                     // number of true negatives (for test spec)
  real<lower=0> coef_prior_scale;                      // prior for scale parameter of random effects of age, seifa, =0.5 as per Gelman paper
  real<lower=0> sd_x3_con;                     // adjustment to scale parameter for linear term for 'postcode'
  int prior_only; // should the likelihood be ignored?
  real<lower=0,upper=1> sens;  // estimated test sensitivity
  vector<lower=1,upper=n_x3>[n] x3_con;
  vector<lower=0,upper=1>[n] x4_bin;
}
parameters{
  vector[K] zb;  // unscaled coefficients
  // parameters of the R2D2 prior
  real<lower=0,upper=1> R2D2_R2;
  simplex[Kscales] R2D2_phi;
  
  real<lower=0,upper=1> spec;                          // estimated test specificity, 
  real<lower=0> sigma_x1;                             // scale parameters for random effects
  real<lower=0> sigma_x2;                             // scale parameters for random effects
  real<lower=0> sigma_x3;                             // scale parameters for random effects
  vector<multiplier=sigma_x1>[n_x1] a_x1;           // varying intercepts for random effects
  vector<multiplier=sigma_x2>[n_x2] a_x2;           // varying intercepts for random effects
  vector<multiplier=sigma_x3>[n_x3] a_x3;           // varying intercepts for random effects
}
transformed parameters {
  vector[K] b;  // scaled coefficients
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
  // likelihood including constants
  if (!prior_only) {
    vector[n] p = inv_logit(b[1] + b[2]*x3_con +
    b[3]*x4_bin +
    a_x1[x1] + 
    a_x2[x2] +
    a_x3[x3]);
    vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
    y ~ bernoulli(p_sample);
  }
  
  // Prior for the global R2D2 scale parameter
  R2D2_R2 ~ beta(R2D2_mean_R2 * R2D2_prec_R2, (1 - R2D2_mean_R2) * R2D2_prec_R2);
  
  // Priors
  zb ~ normal(0, 1);                  // Standard normal prior for unscaled coefficients
  R2D2_phi ~ dirichlet(R2D2_cons_D2); // Dirichlet prior for R2D2_phi
  
  
  // priors including constants
  spec ~ beta(TN+1, FP+1);                             // or use half-normal priors like N+(1,0.02) as per Andrew suggestion?
  a_x1 ~ normal(0, sigma_x1);                        // random effects
  a_x2 ~ normal(0, sigma_x2);                        // random effects
  a_x3 ~ normal(0, sigma_x3);                        // random effects
  b[1] + b[2]*mean(x3_con) + b[3]*mean(x4_bin) ~ logistic(0,1);                // prior on centred intercept, corresponds to U(0,1)
  b[2] ~ normal(0, coef_prior_scale/sd_x3_con);
  b[3] ~ normal(0, coef_prior_scale);
  sigma_x1 ~ normal(0, coef_prior_scale);  // prior for scale parameters for random effects
  sigma_x2 ~ normal(0, coef_prior_scale);           
  sigma_x3 ~ normal(0, coef_prior_scale);           
}
generated quantities {
  vector[n] p = inv_logit(b[1] + b[2]*x3_con +
  b[3]*x4_bin +
  a_x1[x1] + 
  a_x2[x2] + 
  a_x3[x3]);
  vector[n] p_sample = p * sens + (1 - p) * (1 - spec);
  array[n] real y_rep = bernoulli_rng(p_sample);

  real p_avg1 = mean(p);
  real p_samp_avg1 = mean(p_sample);
}
