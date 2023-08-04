library(tidyverse)


# Load simulated data
res_df_MC <- readRDS(file = "simulations/data/sim_crossDGP_20230525.rds")


# Set options
model_choice <- "lin_6p"
q1 <- 3
diag_entries <- which(MTS::Vech(diag(2*q1))==1)

# Estimated Covariance Matrix and standard deviations (diagonal covariance entries only) as data frames
res_MC_cov <- res_df_MC %>%
  filter(model==model_choice, type %in% c("cov_est_asy")) %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure, type) %>%
  summarize(theta_cov_mean = mean(value, na.rm = TRUE),
            theta_cov_median = median(value, na.rm = TRUE))

res_MC_sd <- res_MC_cov %>%
  filter(theta_index %in% diag_entries) %>%
  group_by(model, n, alpha, beta, risk_measure) %>%
  mutate(theta_index = 1:length(theta_index)) %>% # Reset theta index after selecting diagonal entries only
  group_by(theta_index) %>%
  mutate(sd_asy_mean = sqrt(theta_cov_mean),
         sd_asy_median = sqrt(theta_cov_median))

# True parameter values
res_MC_true <- res_df_MC %>%
  filter(model==model_choice, type=="true_value") %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure) %>%
  summarize(theta_true = mean(value,na.rm = TRUE))

# Extract the mean parameter estimates
res_MC_est <- res_df_MC %>%
  dplyr::filter(model==model_choice, type=="param_est") %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure) %>%
  summarize(theta_mean = mean(value),
            theta_median = median(value),
            theta_sd_emp = sd(value)) %>%
  ungroup()

# Join data frames together
res_MC <- res_MC_est %>%
  left_join(res_MC_true,
            by=c("theta_index", "model", "n", "alpha", "beta", "risk_measure")) %>%
  left_join(res_MC_sd,
            by=c("theta_index", "model", "n", "alpha", "beta", "risk_measure")) %>%
  dplyr::select(theta_index, model, risk_measure, beta, alpha, n, theta_true, theta_mean, theta_median, theta_sd_emp, sd_asy_mean, sd_asy_median)



## Confidence interval coverage
res_MC_CI <- res_df_MC %>%
  filter(model==model_choice, type %in% c("cov_est_asy")) %>%
  group_by(i_MC, model, n, alpha, beta, risk_measure, type) %>%
  filter(theta_index %in% diag_entries) %>%
  mutate(theta_index=1:length(diag_entries)) %>% # This line is dangerous!
  reshape2::dcast(i_MC + theta_index + model + n + alpha + beta + risk_measure ~ type) %>%
  left_join(res_df_MC %>% dplyr::filter(type=="param_est") %>% dplyr::select(-type) %>% rename(param_est = value),
            by=c("i_MC", "theta_index", "model", "n", "alpha", "beta", "risk_measure")) %>%
  mutate(CI_asy_lower = param_est - qnorm(0.975)*sqrt(cov_est_asy),
         CI_asy_upper = param_est + qnorm(0.975)*sqrt(cov_est_asy)) %>%
  left_join(res_df_MC %>% filter(type=="true_value") %>% dplyr::select(-type) %>% rename(true_value = value),
            by=c("i_MC", "theta_index", "model", "n", "alpha", "beta", "risk_measure")) %>%
  as_tibble()


# Compute CI coverage
res_MC_CIcoverage <- res_MC_CI %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure) %>%
  summarize(CI_asy_coverage = mean( (true_value >= CI_asy_lower) & (true_value <= CI_asy_upper)),
            CI_asy_length = median(CI_asy_upper - CI_asy_lower))





# Look at specific results
res_MC %>% filter(n==4000, beta==0.975) 

res_MC_CIcoverage %>% filter(n==4000) %>% arrange(beta) %>% print(n=50)












##########
##   Automatically print results to LaTeX tables
##########



# Joint table displaying both, consistency and covariance estimation results
res_MC %>%
  dplyr::mutate(bias_mean=(theta_mean-theta_true),
                bias_median=(theta_median-theta_true)) %>%
  dplyr::select(theta_index, beta, n, bias_mean, theta_sd_emp, sd_asy_mean) %>%
  full_join(res_MC_CIcoverage %>% dplyr::select(c(theta_index, n, beta, CI_asy_coverage)),
            by=c("theta_index", "beta", "n")) %>%
  dplyr::filter(n >= 500, beta <= 0.99) %>%
  mutate(theta_type=ifelse(theta_index<=q1, "VaR", "MES"),
         theta_index_partial=ifelse(theta_index<=q1, theta_index, theta_index-q1))  %>%   # This appears weird, transform MES parameters with number (3,4) to (1,2)...
  tidyr::pivot_wider(id_cols = c("beta", "n", "theta_type"),
                     names_from = theta_index_partial,
                     values_from = c(bias_mean, theta_sd_emp, sd_asy_mean, CI_asy_coverage)) %>%
  arrange(desc(theta_type), beta, n) %>%
  mutate(empty1=NA, empty2=NA, empty3=NA) %>%
  dplyr::select(beta, n, empty1, bias_mean_1, theta_sd_emp_1, sd_asy_mean_1, CI_asy_coverage_1,
                empty2, bias_mean_2, theta_sd_emp_2, sd_asy_mean_2, CI_asy_coverage_2,
                empty3, bias_mean_3, theta_sd_emp_3, sd_asy_mean_3, CI_asy_coverage_3) %>%
  xtable::xtable(digits=c(0, 3, 0, 0, c(3,3,3,2), 0, c(3,3,3,2), 0, c(3,3,3,2) )) %>%
  print(file="simulations/output/MES_AllResults.txt", include.rownames=FALSE, booktabs=TRUE)




######## Compare with the regressions from the Brunnermeier et al (BRS) papers
# Extract the BRS parameter estimates
res_MC_BRSest <- res_df_MC %>%
  dplyr::filter(model==model_choice, type=="BRSreg_est") %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure) %>%
  summarize(theta_mean = mean(value),
            theta_median = median(value),
            theta_sd_emp = sd(value)) %>%
  ungroup()


# Plot histograms
beta_choice <- 0.95
res_df_plot <- res_df_MC %>% 
  dplyr::filter(n==4000, 
                model==model_choice, 
                theta_index %in% 4:6, 
                type %in% c("BRSreg_est", "param_est", "true_value"), 
                beta==beta_choice) %>% as_tibble()
  
res_df_plot$theta_index <- factor(res_df_plot$theta_index,
                                  labels=c(expression(theta[1]^m), expression(theta[2]^m), expression(theta[3]^m)))
 

# res_df_plot$type <- factor(res_df_plot$type, 
#        labels = c("MES regression", expression(paste("Mean regression of ", Y^ast)), "true_value"))



ggplot(res_df_plot %>% dplyr::filter(type!="true_value")) + 
  geom_histogram(aes(x=value, y = after_stat(density), fill=type), bins=100) + 
  geom_vline(data=res_df_plot %>% dplyr::filter(type=="true_value"), 
             aes(xintercept=value), color="black") +
  facet_wrap(~theta_index, scales="free", labeller = label_parsed) +
  scale_fill_manual(values = c('red','blue'),
                    name = '', 
                    labels = expression(paste(paste("Mean regression of ", Y[t],"*     ")), "MES regression    ")) +
  xlab("Parameter value") +
  theme_bw() +
  theme(legend.position="bottom")

ggsave("simulations/output/ParamHist.pdf", width = 8, height = 3, units = "in")



