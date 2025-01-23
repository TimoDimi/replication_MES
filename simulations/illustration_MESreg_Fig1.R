
library(MASS)
library(dplyr)
library(tibble)
library(ggplot2)
library(doParallel)
library(MTS)
library(tsibble)
library(SystemicRisk)
library(slider)

# Load DGP function
source("simulations/SRM_DGP.R")

gamma <- c(1, 1.5, 0, 0.25, 1, 0)
sigma1 <- 2
sigma2 <- 1
rho <- 0.5
Sigma2 <- matrix(c(sigma1^2, rho*sigma1*sigma2, rho*sigma1*sigma2, sigma2^2), nrow=2) # Covariance matrix
Sigma <- t(chol(Sigma2)) # Square-root of Covariance matrix
phi <- c(0.5, 0.8)
nu <- 8

beta <- 0.9
alpha <- 0.95

sim_help <- sim_SRM_DGP(n=1000, gamma=gamma, Sigma=Sigma, nu=nu, phi=phi, beta=beta, alpha=alpha)


df_plot <- sim_help$data %>%
  mutate(Expectation = gamma[1] + gamma[2] * z2,
         Quantile = sim_help$theta_true$VaR[1] + sim_help$theta_true$VaR[2] * z2,
         QViol = (x >= Quantile),
         ES = 1.2*sim_help$theta_true$VaR[1] + 1.3*sim_help$theta_true$VaR[2] * z2,
         MES = sim_help$theta_true$MES[1] + sim_help$theta_true$MES[2] * z2,
         CoVaR = sim_help$theta_true$CoVaR[1] + sim_help$theta_true$CoVaR[2] * z2) %>%
  na.omit()



# Conditional MES
df_plot_MES <- bind_rows(df_plot %>% dplyr::select(z2, QViol, x, Quantile) %>% mutate(Panel="X / VaR") %>% rename(Outcome=x, RiskMeasure=Quantile),
                          df_plot %>% dplyr::select(z2, QViol, y, MES)  %>% mutate(Panel="Y / MES") %>% rename(Outcome=y, RiskMeasure=MES)) %>%
  arrange(QViol) %>%
  rename("VaR exceedance" = "QViol")


ggplot(df_plot_MES) +
  geom_point(aes(x=z2, y=Outcome, col=`VaR exceedance`)) +
  ggplot2::scale_colour_manual(values = c("grey", "black")) +
  ggnewscale::new_scale_color() + # For two color scales!!!
  geom_line(aes(x=z2, y=RiskMeasure, col=Panel)) +
  ggplot2::scale_colour_manual(values = c("blue", "red")) +
  facet_wrap(~Panel, nrow=2, scales="free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Covariate")

ggsave("simulations/output/illustration_MESreg.pdf", 
       width = 8, height = 5.5, units = "in")

