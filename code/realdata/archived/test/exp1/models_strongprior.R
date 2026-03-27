## multiplying beta values for sens and spec

# model 6 - with all two-way interactions between SEIFA, sex, agegp  --------------------
file_m6 <- file.path(here("code/stan/model6.stan"))
mod_m6 <- cmdstan_model(file_m6)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m6 <- mod_m6$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m6$save_object(file = here("results/experiments/exp5/m6_postpred_b10_sens.RDS"), compress=T)


# model 5 - adapting Marnie's code to do prior predictive checks --------------------
file_m5 <- file.path(here("code/stan/model5.stan"))
mod_m5 <- cmdstan_model(file_m5)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m5 <- mod_m5$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m5$save_object(file = here("results/experiments/exp5/m5_postpred_b10_sens.RDS"), compress=T)

# model 4 - without agegp -------------------------------------------------
file_m4 <- file.path(here("code/stan/model4.stan"))
mod_m4 <- cmdstan_model(file_m4)

## sampling ####
# post pred
data_list$prior_only = 0

res_posterior_m4 <- mod_m4$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m4$save_object(file=here("results/experiments/exp5/m4_postpred_b10_sens.RDS"), compress=T)


# model 3 without agegp, sex ----------------------------------------------
file_m3 <- file.path(here("code/stan/model3.stan"))
mod_m3 <- cmdstan_model(file_m3)

## sampling ####
## posterior pred
data_list$prior_only = 0

res_posterior_m3 <- mod_m3$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m3$save_object(file=here("results/experiments/exp5/m3_postpred_b10_sens.RDS"), compress=T)


# model 2 without agegp, sex, seifa ----------------------------------------------
file_m2 <- file.path(here("code/stan/model2.stan"))
mod_m2 <- cmdstan_model(file_m2)

## sampling ####
## post pred
data_list$prior_only = 0

res_posterior_m2 <- mod_m2$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m2$save_object(file=here("results/experiments/exp5/m2_postpred_b10_sens.RDS"), compress=T)


# model 1 without agegp, sex, seifa, postcode ----------------------------------------------
file_m1 <- file.path(here("code/stan/model1.stan"))
mod_m1 <- cmdstan_model(file_m1)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_m1 <- mod_m1$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m1$save_object(file=here("results/experiments/exp5/m1_postpred_b10_sens.RDS"), compress=T)
```


```{r, include=F}
# Multiplying beta values x10 for specificity only:
(data_list$TN = data_list$TN*10 )
(data_list$FP = data_list$FP*10 )
(data_list$TP = data_list$TP/10 ) # converting back to original
(data_list$FN = data_list$FN/10 ) # converting back to original
```

```{r spec x10 models, include=F, cache=T}
# model 6 - with all two-way interactions between SEIFA, sex, agegp  --------------------
file_m6 <- file.path(here("code/stan/model6.stan"))
mod_m6 <- cmdstan_model(file_m6)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m6 <- mod_m6$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m6$save_object(file = here("results/experiments/exp5/m6_postpred_b10_spec.RDS"), compress=T)


# model 5 - adapting Marnie's code to do prior predictive checks --------------------
file_m5 <- file.path(here("code/stan/model5.stan"))
mod_m5 <- cmdstan_model(file_m5)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m5 <- mod_m5$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m5$save_object(file = here("results/experiments/exp5/m5_postpred_b10_spec.RDS"), compress=T)

# model 4 - without agegp -------------------------------------------------
file_m4 <- file.path(here("code/stan/model4.stan"))
mod_m4 <- cmdstan_model(file_m4)

## sampling ####
# post pred
data_list$prior_only = 0

res_posterior_m4 <- mod_m4$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m4$save_object(file=here("results/experiments/exp5/m4_postpred_b10_spec.RDS"), compress=T)


# model 3 without agegp, sex ----------------------------------------------
file_m3 <- file.path(here("code/stan/model3.stan"))
mod_m3 <- cmdstan_model(file_m3)

## sampling ####
## posterior pred
data_list$prior_only = 0

res_posterior_m3 <- mod_m3$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m3$save_object(file=here("results/experiments/exp5/m3_postpred_b10_spec.RDS"), compress=T)


# model 2 without agegp, sex, seifa ----------------------------------------------
file_m2 <- file.path(here("code/stan/model2.stan"))
mod_m2 <- cmdstan_model(file_m2)

## sampling ####
## post pred
data_list$prior_only = 0

res_posterior_m2 <- mod_m2$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m2$save_object(file=here("results/experiments/exp5/m2_postpred_b10_spec.RDS"), compress=T)


# model 1 without agegp, sex, seifa, postcode ----------------------------------------------
file_m1 <- file.path(here("code/stan/model1.stan"))
mod_m1 <- cmdstan_model(file_m1)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_m1 <- mod_m1$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m1$save_object(file=here("results/experiments/exp5/m1_postpred_b10_spec.RDS"), compress=T)


```

```{r, include=F}
# Multiplying beta values x10 for both specificity and sensitivity:
(data_list$TP = data_list$TP*10 ) # multiply sens by 10 
(data_list$FN = data_list$FN*10 ) # multiply sens by 10 

```

```{r spec and sens x10 models, include=F, cache=T}
# model 6 - with all two-way interactions between SEIFA, sex, agegp  --------------------
file_m6 <- file.path(here("code/stan/model6.stan"))
mod_m6 <- cmdstan_model(file_m6)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m6 <- mod_m6$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m6$save_object(file = here("results/experiments/exp5/m6_postpred_b10_both.RDS"), compress=T)


# model 5 - adapting Marnie's code to do prior predictive checks --------------------
file_m5 <- file.path(here("code/stan/model5.stan"))
mod_m5 <- cmdstan_model(file_m5)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m5 <- mod_m5$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m5$save_object(file = here("results/experiments/exp5/m5_postpred_b10_both.RDS"), compress=T)

# model 4 - without agegp -------------------------------------------------
file_m4 <- file.path(here("code/stan/model4.stan"))
mod_m4 <- cmdstan_model(file_m4)

## sampling ####
# post pred
data_list$prior_only = 0

res_posterior_m4 <- mod_m4$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m4$save_object(file=here("results/experiments/exp5/m4_postpred_b10_both.RDS"), compress=T)


# model 3 without agegp, sex ----------------------------------------------
file_m3 <- file.path(here("code/stan/model3.stan"))
mod_m3 <- cmdstan_model(file_m3)

## sampling ####
## posterior pred
data_list$prior_only = 0

res_posterior_m3 <- mod_m3$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m3$save_object(file=here("results/experiments/exp5/m3_postpred_b10_both.RDS"), compress=T)


# model 2 without agegp, sex, seifa ----------------------------------------------
file_m2 <- file.path(here("code/stan/model2.stan"))
mod_m2 <- cmdstan_model(file_m2)

## sampling ####
## post pred
data_list$prior_only = 0

res_posterior_m2 <- mod_m2$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m2$save_object(file=here("results/experiments/exp5/m2_postpred_b10_both.RDS"), compress=T)


# model 1 without agegp, sex, seifa, postcode ----------------------------------------------
file_m1 <- file.path(here("code/stan/model1.stan"))
mod_m1 <- cmdstan_model(file_m1)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_m1 <- mod_m1$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m1$save_object(file=here("results/experiments/exp5/m1_postpred_b10_both.RDS"), compress=T)

```

```{r, include=F}
# Multiplying beta values x100: 
# beta values * 100 -------------------------------------------------------
(data_list$TP = data_list$TP*10 )
(data_list$FN = data_list$FN*10 )
(data_list$TN = data_list$TN*10 )
(data_list$FP = data_list$FP*10 )
```

```{r x100 models, include=F, cache=T}
# model 6 - with all two-way interactions between SEIFA, sex, agegp  --------------------
file_m6 <- file.path(here("code/stan/model6.stan"))
mod_m6 <- cmdstan_model(file_m6)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m6 <- mod_m6$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 5 - adapting Marnie's code to do prior predictive checks --------------------
file_m5 <- file.path(here("code/stan/model5.stan"))
mod_m5 <- cmdstan_model(file_m5)

## sampling ####
# posterior predictive checks 
data_list$prior_only = 0

res_posterior_m5 <- mod_m5$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 4 - without agegp -------------------------------------------------
file_m4 <- file.path(here("code/stan/model4.stan"))
mod_m4 <- cmdstan_model(file_m4)

## sampling ####
# post pred
data_list$prior_only = 0

res_posterior_m4 <- mod_m4$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 3 without agegp, sex ----------------------------------------------
file_m3 <- file.path(here("code/stan/model3.stan"))
mod_m3 <- cmdstan_model(file_m3)

## sampling ####
## posterior pred
data_list$prior_only = 0

res_posterior_m3 <- mod_m3$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 2 without agegp, sex, seifa ----------------------------------------------
file_m2 <- file.path(here("code/stan/model2.stan"))
mod_m2 <- cmdstan_model(file_m2)

## sampling ####
## post pred
data_list$prior_only = 0

res_posterior_m2 <- mod_m2$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# model 1 without agegp, sex, seifa, postcode ----------------------------------------------
file_m1 <- file.path(here("code/stan/model1.stan"))
mod_m1 <- cmdstan_model(file_m1)

## sampling ####
# posterior predictive
data_list$prior_only = 0

res_posterior_m1 <- mod_m1$sample(
  data = data_list, 
  seed = 2345, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500, # print update every 500 iters
)

# saving
res_posterior_m6$save_object(file = here("results/experiments/exp5/m6_postpred_betavalues*100.RDS"), compress=T)
res_posterior_m5$save_object(file = here("results/experiments/exp5/m5_postpred_betavalues*100.RDS"), compress=T)
res_posterior_m4$save_object(file=here("results/experiments/exp5/m4_postpred_betavalues*100.RDS"), compress=T)
res_posterior_m3$save_object(file=here("results/experiments/exp5/m3_postpred_betavalues*100.RDS"), compress=T)
res_posterior_m2$save_object(file=here("results/experiments/exp5/m2_postpred_betavalues*100.RDS"), compress=T)
res_posterior_m1$save_object(file=here("results/experiments/exp5/m1_postpred_betavalues*100.RDS"), compress=T)