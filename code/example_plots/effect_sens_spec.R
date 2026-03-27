library(tidyverse)
library(colorspace)

df = expand.grid(sensitivity_delta = seq(0,1,.1),
                 specificity_gamma = seq(0,1,.1),
                 true_prevalence = c(.001,.01,.1,.4,.5,.6,.9,.99,.999))

df <- df %>%
  mutate(observed_test = (1-specificity_gamma)*(1-true_prevalence) +sensitivity_delta*true_prevalence,
         bias = observed_test - true_prevalence)

ggplot(df, aes(x = specificity_gamma, y = sensitivity_delta, fill = bias))+
  geom_tile()+
  facet_wrap(.~true_prevalence) +
  scale_fill_continuous_diverging(palette = "Purple-Brown")+
  xlab((expression("Specificity ("~gamma~")")))+
  ylab((expression("Sensitivity ("~delta~")")))+
  labs(fill = expression("Test - Disease prevalence: Pr(y* = 1) - "~pi~""))+
  theme(legend.position = "bottom")

ggsave(filename = "code/example_plots/figures/effect_sens_spec.png", width =13, height = 12, units = "cm")
