## using BART instead of MR in the MRP 
## as per Qixuan's suggestion 
## 16/02/2023
library(dbarts)
library(here)
source(here("code/data_list.R"), echo=F)

y = data_list$y 
N_pop1 = data_list$N_pop1

## model 1 ####
x_m1 = cbind(data_list$strata)  |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

x_pop_m1 = cbind(data_list$strata_pop) |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

bartfit_m1 = bart(x_m1, y, keeptrees=T) 
pop_m1 = predict(bartfit_m1, newdata = x_pop_m1)
p_pop_m1_quant =  apply(pop_m1, 2, quantile, c(0.05, 0.5, 0.95)) # getting posterior quantiles 

## model 2 ####
x_m2 = cbind(data_list$strata, data_list$postcode)  |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

x_pop_m2 = cbind(data_list$strata_pop, data_list$postcode_pop) |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

bartfit_m2 = bart(x_m2, y, keeptrees=T) 
pop_m2 = predict(bartfit_m2, newdata = x_pop_m2)
p_pop_m2_quant =apply(pop_m2, 2, quantile, c(0.05, 0.5, 0.95)) # getting posterior quantiles 

## model 3 ####
x_m3 = cbind(data_list$strata, data_list$postcode, data_list$seifa)  |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

x_pop_m3 = cbind(data_list$strata_pop, data_list$postcode_pop, data_list$seifa_pop) |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

bartfit_m3 = bart(x_m3, y, keeptrees=T) 
pop_m3 = predict(bartfit_m3, newdata = x_pop_m3)
p_pop_m3_quant =  apply(pop_m3, 2, quantile, c(0.05, 0.5, 0.95)) # getting posterior quantiles 

## model 4 ####
# data
x_m4 = cbind(data_list$strata, data_list$postcode, data_list$seifa, data_list$sex) |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

x_pop_m4 = cbind(data_list$strata_pop, data_list$postcode_pop, data_list$seifa_pop, data_list$sex_pop) |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

# bart
bartfit_m4 = bart(x_m4, y, keeptrees=T) 
pop_m4 = predict(bartfit_m4, newdata = x_pop_m4)
p_pop_m4_quant = apply(pop_m4, 2, quantile, c(0.05, 0.5, 0.95)) # getting posterior quantiles 

## model 5 ####
x_m5 = cbind(data_list$strata, data_list$postcode, data_list$seifa, data_list$sex, data_list$agegp) |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

x_pop_m5 = cbind(data_list$strata_pop, data_list$postcode_pop, data_list$seifa_pop, data_list$sex_pop, data_list$agegp_pop) |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

# bart
bartfit_m5 = bart(x_m5, y, keeptrees=T) 
pop_m5 = predict(bartfit_m5, newdata = x_pop_m5) # predict for poststrat table 
p_pop_m5_quant = apply(pop_m5, 2, quantile, c(0.05, 0.5, 0.95)) # getting posterior quantiles 

## model 6 ####
x_m6 = as.matrix(cbind(data_list$strata, data_list$postcode, data_list$seifa, data_list$sex, data_list$agegp,
                       data_list$seifa_sex, data_list$seifa_agegp, data_list$sex_agegp))|> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

x_pop_m6 = cbind(data_list$strata_pop, data_list$postcode_pop, data_list$seifa_pop, data_list$sex_pop, data_list$agegp_pop,
                 data_list$seifa_sex_pop, data_list$seifa_agegp_pop, data_list$sex_agegp_pop) |> 
  as_tibble() |>
  mutate_if(is.numeric,as.factor) # factorise 

bartfit_m6 = bart(x_m6, y, keeptrees=T)
pop_m6 = predict(bartfit_m6, newdata = x_pop_m6) # predict for poststrat table 
p_pop_m6_quant = apply(pop_m6, 2, quantile, c(0.05, 0.5, 0.95)) # getting posterior quantiles 

# poststrat estimates 
( p_avg_m1 = c(sum(N_pop1 * p_pop_m1_quant[1,])/sum(N_pop1) *100, sum(N_pop1 * p_pop_m1_quant[3,])/sum(N_pop1) *100) )
( p_avg_m2 = c(sum(N_pop1 * p_pop_m2_quant[1,])/sum(N_pop1) *100, sum(N_pop1 * p_pop_m2_quant[3,])/sum(N_pop1) *100) )
( p_avg_m3 = c(sum(N_pop1 * p_pop_m3_quant[1,])/sum(N_pop1) *100, sum(N_pop1 * p_pop_m3_quant[3,])/sum(N_pop1) *100) )
( p_avg_m4 = c(sum(N_pop1 * p_pop_m4_quant[1,])/sum(N_pop1) *100, sum(N_pop1 * p_pop_m4_quant[3,])/sum(N_pop1) *100) )
( p_avg_m5 = c(sum(N_pop1 * p_pop_m5_quant[1,])/sum(N_pop1) *100, sum(N_pop1 * p_pop_m5_quant[3,])/sum(N_pop1) *100) )
( p_avg_m6 = c(sum(N_pop1 * p_pop_m6_quant[1,])/sum(N_pop1) *100, sum(N_pop1 * p_pop_m6_quant[3,])/sum(N_pop1) *100) )
