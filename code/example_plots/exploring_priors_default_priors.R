library(tidyverse)
library(brms)
library(posterior)
library(ggridges)
library(colorspace)
library(pROC)
df <- data.frame(y = rbinom(1000,1,.1),
                 y_10 = rbinom(1000,1,.1),
                 y_25 = rbinom(1000,1,.25),
                 y_50 = rbinom(1000,1,.5),
                 y_75 = rbinom(1000,1,.75),
                 y_90 = rbinom(1000,1,.9),
                 x1 = factor(sample.int(5,1000,replace = TRUE)),
                 x2 = factor(sample.int(5,1000,replace = TRUE)),
                 x3 = factor(sample.int(5,1000,replace = TRUE)),
                 x4 = factor(sample.int(5,1000,replace = TRUE)),
                 x5 = factor(sample.int(5,1000,replace = TRUE)),
                 x5_many = factor(sample.int(50,1000,replace = TRUE)))

model_formulas <- c("Intercept" = formula(y~1),
                    "One_Predictor" = formula(y~(1|x1)),
                    "Two_Predictors" = formula(y~(1|x1) + (1|x2)),
                    "Three_Predictors" = formula(y~(1|x1) + (1|x2) +(1|x3)),
                    "Four_Predictors" = formula(y~(1|x1) + (1|x2) + (1|x3) + (1|x4)),
                    "Five_Predictors" = formula(y~(1|x1) + (1|x2) + (1|x3) + (1|x4) +(1|x5)))

fit_model <- map(model_formulas, function(x){ brm(x, family = bernoulli(link = "logit"), sample_prior = "only", data = df, backend = "cmdstanr")})

make_prediction <- map(fit_model, posterior_predict)
make_prediction_prob <- map(fit_model, posterior_linpred, transform = TRUE)

#### Prior probability of average p(y = 1) ####
predicted_probability <- map(make_prediction, function(x) apply(x, 1, mean))

unit_predicted_probability <- map(make_prediction_prob, function(x) x[,1])

unit_predicted_probability_df <- as.data.frame(unit_predicted_probability,col.names = names(model_formulas)) %>%
  pivot_longer(everything(), names_to = "model", values_to = "draw") %>%
  mutate(model = ordered(model, levels = c("Intercept",
                                           "One_Predictor",
                                           "Two_Predictors",
                                           "Three_Predictors",
                                           "Four_Predictors",
                                           "Five_Predictors")))%>%
  mutate(ppc = "Unit PPC")

predicted_probability_df <- as.data.frame(predicted_probability,col.names = names(model_formulas)) %>%
  pivot_longer(everything(), names_to = "model", values_to = "draw") %>%
  mutate(model = ordered(model, levels = c("Intercept",
                                           "One_Predictor",
                                           "Two_Predictors",
                                           "Three_Predictors",
                                           "Four_Predictors",
                                           "Five_Predictors")))%>%
  mutate(ppc = "Population Expectation PPC")

ggplot(predicted_probability_df %>%
         rbind(unit_predicted_probability_df), aes(x = draw))+
  geom_histogram()+
  facet_grid(ppc~model)+
  xlab("Expected mean p(y=1)") +
  ggtitle("Default Priors Intercept ~ t(3, 0, 2.5)")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("code/example_plots/default_prior_predicted_prob.png", width =15, height = 10, units = "cm")

#### Score prior ####
calculate_specificity <- function(z, yobs = yobs_variable) {
  return(mean(z[df[[yobs]] == 0]))
}

calculate_sensitivity <- function(z, yobs = yobs_variable) {
  return(mean(z[df[[yobs]] == 1]))
}

calculate_accuracy <- function(z, yobs = yobs_variable) {
  return(mean(z == df[[yobs]]))
}

calculate_auc <- function(z, yobs = yobs_variable) {
  return(auc(roc(df[[yobs]], z)))
}



calculate_prior_score <- function(make_prediction, model_formulas, yobs_variable, function_type, metric_type){
  prior_check <- map(make_prediction, function(x) apply(x, 1, function_type,yobs = yobs_variable))
  prior_df <- as.data.frame(prior_check,col.names = names(model_formulas)) %>%
    pivot_longer(everything(), names_to = "model", values_to = "draw") %>%
    mutate(model = ordered(model, levels = c("Intercept",
                                             "One_Predictor",
                                             "Two_Predictors",
                                             "Three_Predictors",
                                             "Four_Predictors",
                                             "Five_Predictors")))%>%
    mutate(metric = metric_type)
}
calculate_suite_of_scores <- function(make_prediction, model_formulas, name_y_variable){
  test_specificity <- calculate_prior_score(make_prediction,model_formulas, 
                                            yobs_variable = name_y_variable, 
                                            function_type = calculate_specificity,
                                            metric_type = "Specificity")
  test_sensitivity <- calculate_prior_score(make_prediction,model_formulas, 
                                            yobs_variable = name_y_variable, 
                                            function_type = calculate_sensitivity,
                                            metric_type = "Sensitivity")
  test_accuracy <- calculate_prior_score(make_prediction,model_formulas, 
                                         yobs_variable = name_y_variable, 
                                         function_type = calculate_accuracy,
                                         metric_type = "Accuracy")
  test_auc <- calculate_prior_score(make_prediction,model_formulas, 
                                    yobs_variable = name_y_variable, 
                                    function_type = calculate_auc,
                                    metric_type = "AUC")
  full_suite <- rbind(test_specificity, test_sensitivity, test_accuracy, test_auc) %>%
    mutate(y_variable = name_y_variable)
  return(full_suite)
}         

get_all_scores_y_10 <- calculate_suite_of_scores(make_prediction, model_formulas, name_y_variable = "y_10")
get_all_scores_y_25 <- calculate_suite_of_scores(make_prediction, model_formulas, name_y_variable = "y_25")
get_all_scores_y50 <- calculate_suite_of_scores(make_prediction, model_formulas, name_y_variable = "y_50")
get_all_scores_y_75 <- calculate_suite_of_scores(make_prediction, model_formulas, name_y_variable = "y_75")
get_all_scores_y90 <- calculate_suite_of_scores(make_prediction, model_formulas, name_y_variable = "y_90")

all_scores <- rbind(get_all_scores_y_10,get_all_scores_y_25, get_all_scores_y50,get_all_scores_y_75,  get_all_scores_y90)

ggplot(all_scores %>% 
         filter(metric== "Accuracy"),
       aes(x = draw))+
  geom_histogram()+
  facet_grid(y_variable~model)+
  xlab("Expected accuracy p(y = hat(y))") +
  xlim(c(0,1))+
  ggtitle("Default Priors Intercept ~ t(3, 0, 2.5)")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("code/example_plots/default_prior_accuracy.png", width =15, height = 10, units = "cm")

ggplot(all_scores %>% 
         filter(metric== "Sensitivity"),
       aes(x = draw))+
  geom_histogram()+
  facet_grid(y_variable~model)+
  xlab("Expected Sensitivity p(hat(y) = 1|y=1)")+
  ggtitle("Default Priors Intercept ~ t(3, 0, 2.5)")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("code/example_plots/default_prior_predicted_sensitivity.png", width =15, height = 10, units = "cm")

ggplot(all_scores %>% 
         filter(metric== "Specificity"),
       aes(x = draw))+
  geom_histogram()+
  facet_grid(y_variable~model)+
  xlab("Expected Specificity p(hat(y) = 0|y=0))")+
  ggtitle("Default Priors Intercept ~ t(3, 0, 2.5)")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("code/example_plots/default_prior_predicted_specificity.png", width =15, height = 10, units = "cm")

ggplot(all_scores %>% 
         filter(metric== "AUC"),
       aes(x = draw))+
  geom_histogram()+
  geom_vline(xintercept = 0.5) +
  facet_grid(y_variable~model)+
  xlab("Expected AUC")+
  ggtitle("Default Priors Intercept ~ t(3, 0, 2.5)")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("code/example_plots/default_prior_predicted_auc.png", width =15, height = 10, units = "cm")
