res_wd <- "results/experiments/exp1/"
# reading in data 
res_m0_nopsamp = readRDS(file=here(paste0(res_wd,"m0nopsamp.RDS")))
res_m1_nopsamp = readRDS(file=here(paste0(res_wd,"m1nopsamp.RDS")))
res_m2_nopsamp = readRDS(file=here(paste0(res_wd,"m2nopsamp.RDS")))
res_m3_nopsamp = readRDS(file=here(paste0(res_wd,"m3nopsamp.RDS")))
res_m4_nopsamp = readRDS(file=here(paste0(res_wd,"m4nopsamp.RDS")))
res_m5_nopsamp = readRDS(file=here(paste0(res_wd,"m5nopsamp.RDS")))
res_m6_nopsamp = readRDS(file=here(paste0(res_wd,"m6nopsamp.RDS")))


res_est_nopsamp = rbind(res_m0_nopsamp$summary("p_avg1")[,c('mean', 'q5', 'q95')],
                        res_m1_nopsamp$summary("p_avg1")[,c('mean', 'q5', 'q95')],
                        res_m2_nopsamp$summary("p_avg1")[,c('mean', 'q5', 'q95')],
                        res_m3_nopsamp$summary("p_avg1")[,c('mean', 'q5', 'q95')],
                        res_m4_nopsamp$summary("p_avg1")[,c('mean', 'q5', 'q95')],
                        res_m5_nopsamp$summary("p_avg1")[,c('mean', 'q5', 'q95')],
                        res_m6_nopsamp$summary("p_avg1")[,c('mean', 'q5', 'q95')])
res_est_nopsamp <- res_est_nopsamp |> 
  mutate(model = 0:6) |> 
  rename(mrp_mean = mean, mrp_q5 = q5, mrp_q95 = q95)

# sample_predictions and intercept
func_summ <- function(x){
  x$summary(c("p", "b[1]"))[,c('variable', 'mean', 'q5', 'q95')]}

pop_pred <- lapply(list(res_m1_nopsamp, res_m2_nopsamp,
                        res_m3_nopsamp, res_m4_nopsamp,
                        res_m5_nopsamp, res_m6_nopsamp), func_summ)

pop_list = lapply(pop_pred, function(x) {x |> 
    mutate(var = gsub("(p)(\\[[0-9]*\\])", "\\1", variable)) |>
    select(-variable) |> 
    group_by(var) |> 
    summarise_all(mean)})  

tab_nop_samp = do.call(rbind,pop_list) |> 
  as_tibble() |> 
  mutate(model = rep(1:6, each=2)) |> 
  pivot_wider(names_from = var, values_from = c(mean, q5, q95)) |> 
  rename(int_mean = `mean_b[1]`, int_q5 = `q5_b[1]`, int_q95 = `q95_b[1]`,
         indPred_mean = mean_p, indPred_q5 = q5_p, indPred_q95 = q95_p) |> 
  right_join(res_est_nopsamp, ., by='model') 


t1 = tab_nop_samp |> 
  pivot_longer(c(indPred_q5, mrp_q5, int_q5), names_to = c("type", "quantile"),
               names_sep = "_",
               values_to = "q5") |>
  select(model,type,q5) 

t2 = tab_nop_samp |> 
  pivot_longer(c(indPred_q95, mrp_q95, int_q95), names_to = c("type", "quantile"),
               names_sep = "_",
               values_to = "q95") |> 
  select(model,type,q95) 

t3 = left_join(t1,t2, by=c('model', "type"))

tab_nop_samp |> 
  pivot_longer(c(indPred_mean, mrp_mean, int_mean), names_to = c("type","measure"),
               names_sep = "_",
               values_to = "mean") %>%
  select(-measure) |> 
  left_join(t3, by=c('model', "type")) |> 
  select(model, type:q95) |> 
  filter(type != "int") |> 
  ggplot(aes(x=model, y = mean, col=type)) + 
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_line(position=position_dodge(width=0.5)) +
  geom_errorbar(position=position_dodge(width=0.5),aes(x=model,ymin=q5, ymax=q95), width=.2, alpha=0.4) +
  scale_color_manual(values=c("#F28E21", "#59A14F", "#59E01F", "#756BB1", "#54278F")   ) +
  geom_hline(yintercept=0.016, lty=2, col="darkgrey") +
  theme_minimal() +
  ylab('Estimate') +
  ggtitle("Real data: no sensitivity and specificity") 

## intercept only ### 
tab_nop_samp |> 
  pivot_longer(c(indPred_mean, mrp_mean, int_mean), names_to = c("type","measure"),
               names_sep = "_",
               values_to = "mean") %>%
  select(-measure) |> 
  left_join(t3, by=c('model', "type")) |> 
  select(model, type:q95) |> 
  filter(type == "int") |> 
  ggplot(aes(x=model, y = mean, col=type)) + 
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_line(position=position_dodge(width=0.5)) +
  geom_errorbar(position=position_dodge(width=0.5),aes(x=model,ymin=q5, ymax=q95), width=.2, alpha=0.4) +
  scale_color_manual(values=c("#000000", "#59A14F", "#59E01F", "#756BB1", "#54278F")   ) +
  theme_minimal() +
  ylab('Intercept estimate') +
  ggtitle("Real data: no sensitivity and specificity") +
  ylim(-5, -2)

