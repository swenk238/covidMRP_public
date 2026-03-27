# to test if we flip the y's, it will be sensitivity that we care about instead of specificity 
library(here)
library(cmdstanr)
# data 
source(here("code/data_list.R"), echo=F)
options(mc.cores = 8)

table(data_list$y)

# flipping y's so that 1's are 0's and 0's are 1's - 
data_list$y_flip = data_list$y+1
data_list$y_flip = ifelse(data_list$y_flip == "2", "0", "1")
table(data_list$y_flip)

data_list$y = as.integer(data_list$y_flip)
table(data_list$y)
data_list["y_flip"] = NULL

source(here("code/experiments/exp3/models_fixedsens.R", echo=TRUE))

source(here("code/experiments/exp3/models_fixedspec.R", echo=TRUE))

source(here("code/experiments/exp3/models_fixedsensspec.R", echo=TRUE))