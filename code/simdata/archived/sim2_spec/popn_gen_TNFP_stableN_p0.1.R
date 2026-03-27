library(rlang) # eval_tidy, parse_expr()
library(arm)
library(tidyverse)
library(here)

# generate a finite population
p_set = 0.1
samp_size = c(400, 4000)
n_level = c(4, 10, 20, 40)
beta1_set = 0.3

popn_size = 500000

a = -0.5 # Uniform distribution
b = 0.5 
cov_num = 5

# solving for beta0, should be ~ logit(p_set) since E(X) = 0
# b0 are fixed given the uniform upp and low and number of covariates 
(b0 = logit(p_set) - (0.3*(a+b)/2)*cov_num) # beta1 = 0.3， # b0 = -4.59

# specifying sens and spec 
sens = 1 

# generate spec from TN and FP values
total_vec <- rep(c(400,800,1200,8000),4)
FP_vec <- c(8,16,24,160, 4,8,12,80, 2,4,6,40, rep(0,4))
TN_vec <- total_vec - FP_vec # corresponding to specificity values of 0.98, 0.99, 0.995 and 1 
spec_dup <- TN_vec / (TN_vec + FP_vec)
spec <- unique(spec_dup)

# generating 
set.seed(2468)
x_cont_mat = matrix(NA, ncol=5, nrow=popn_size)

# generate uniform x's
for (ind in 1:5){
  x_cont_mat[,ind] = runif(popn_size,min=a,max=b)
}
  
# pi_y_gen
pi_y_gen = invlogit(b0 + as.matrix(x_cont_mat[,1:5]) %*% rep(beta1_set,cov_num))

# p_y and y's
# naming by specificity
p_y_name <- paste0('p_y_spec', format(spec, digits=4))
y_name <- paste0('y_spec', format(spec, digits=4))

# assigning values to each of the name
for (k in 1:length(spec)){
  assign(p_y_name[k], pi_y_gen * sens + (1-pi_y_gen) * (1-spec[k]))
  assign(y_name[k], rbinom(popn_size,1, eval_tidy(parse_expr(p_y_name[k]))) )
}

y_pi_y_gen <- rbinom(popn_size,1, pi_y_gen)

# saving as popn 
popn_data = cbind(x_cont_mat, pi_y_gen, y_pi_y_gen,
                  sapply(p_y_name, function(x)eval_tidy(parse_expr(x))),  #extracting variables 
                  sapply(y_name, function(x)eval_tidy(parse_expr(x)))) %>%
  as_tibble()

colnames(popn_data) = c('x1_cont', 'x2_cont', 'x3_cont', 'x4_cont', 
                        'x5_cont', 'pi_y_gen', 'y_pi_y_gen', p_y_name, y_name)

popn_data = popn_data %>% 
  mutate(across(x1_cont:x5_cont, 
                ~cut_number(.x, n = 4, labels=F), 
                .names= "{col}_4levels")) %>%
  mutate(across(x1_cont:x5_cont, 
                ~cut_number(.x, n = 10, labels=F), 
                .names= "{col}_10levels")) %>%
  mutate(across(x1_cont:x5_cont, 
                ~cut_number(.x, n = 20, labels=F), 
                .names= "{col}_20levels")) %>%
  mutate(across(x1_cont:x5_cont, 
                ~cut_number(.x, n = 40, labels=F), 
                .names= "{col}_40levels")) 


# # checking the pi_y_gen and p_y values 
# popn_data |> summarise(mean(pi_y_gen), 
#                        mean(eval_tidy(parse_expr(p_y_name[1]))),
#                        mean(eval_tidy(parse_expr(y_name[4]))))

saveRDS(popn_data, here("code/simdata/sim2_spec/data/popn_data_TNFP_stableN_p0.1.rds"))
