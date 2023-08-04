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


# Preliminary settings
M <- 5000 # 5000 replications took 12 hours on 10 kernels.
optim_replications <- c(1,3)
n_set <- c(500, 1000, 2000, 4000)

prob_level_list <- list(c(beta=0.9, alpha=0.5),
                        c(beta=0.95, alpha=0.5),
                        c(beta=0.975, alpha=0.5),
                        c(beta=0.99, alpha=0.5))

risk_measure_set <- c("MES")
model_set <- c("lin_6p")
gamma_list <- list(c(1, 1.5, 0, 0.25, 0.5, 0),
                   c(1, 1.5, 2, 0.25, 1, 0))

sigma1 <- 2
sigma2 <- 1
rho <- 0.6
Sigma2 <- matrix(c(sigma1^2, rho*sigma1*sigma2, rho*sigma1*sigma2, sigma2^2), nrow=2) # Covariance matrix
Sigma <- t(chol(Sigma2)) # Square-root of Covariance matrix
phi <- c(0.6, 0.75)
nu <- 6


# Cluster Settings
core.max <- 80
cl <- makeCluster(min(parallel::detectCores()-1, M, core.max) )
registerDoParallel(cl)
start.time <- Sys.time()
res_df_MC <- foreach(
  i_MC = 1:M,
  .combine=rbind,
  .packages=c("dplyr", "tibble", "MTS", "MASS", "mvtnorm", "cubature", "lubridate", "abind", "tsibble", "slider", "SystemicRisk"),
  .errorhandling="pass"
)%dopar%{
  set.seed(i_MC) # set seed for reproducibility
  res_df <- tibble()
  
  # Load DGP function
  source("simulations/SRM_DGP.R")
  
  for (model in model_set){
    for (n in n_set){
      for (prob_level in prob_level_list){
        beta <- prob_level[[1]]
        alpha <- prob_level[[2]]
        
        # Simulate data
        gamma <- switch(model,
                        lin_4p = {gamma_list[[1]]},
                        lin_6p = {gamma_list[[2]]})
        sim_help <- sim_SRM_DGP(n=n, gamma=gamma, Sigma=Sigma, nu=nu, phi=phi, beta=beta, alpha=alpha)
        
        # Transform data to tsibble
        data <- sim_help$data %>%
          mutate(Date=1:n()) %>%
          as_tsibble(index=Date)
        
        for (risk_measure in risk_measure_set){
          
          # True parameter, depends on the risk_measure
          theta_true <- unlist(sim_help$theta_true[c("VaR",risk_measure)])
          
          # Parameter estimate
          if (model=="lin_4p"){
            theta0 <- theta_true[-c(3,6)]
            est_obj <- SRM(data = data %>% dplyr::select(-z3),
                           model="joint_linear", risk_measure=risk_measure, beta=beta, alpha=alpha, theta0=theta0, optim_replications=optim_replications)
          } else if (model=="lin_6p"){
            theta0 <- theta_true
            est_obj <- SRM(data=data,
                           model="joint_linear", risk_measure=risk_measure, beta=beta, alpha=alpha, theta0=theta0, optim_replications=optim_replications)
          } else {
            stop("enter a correct model name!")
          }
          
          ### BRS-type MES regression
          S <- 250
          sim_help_BRSreg <- sim_SRM_DGP(n=n+S+2, gamma=gamma, Sigma=Sigma, nu=nu, phi=phi, beta=beta, alpha=alpha)
          
          data_BRSreg <- sim_help_BRSreg$data %>%
            mutate(Date=1:n(),
                   Qx_rolling = slider::slide_dbl(x, quantile, .before = S, .after = 0, beta, na.rm=TRUE),
                   y_x_viol = ifelse((x >= Qx_rolling),y,NA),
                   y_MESrolling = slider::slide_dbl(y_x_viol, mean, .before = S, .after = 0, na.rm=TRUE)) %>%
            tail(n+1)
          
          sum_BRSreg <- summary(lm(y_MESrolling ~ z2+z3, data_BRSreg))
          
          
          df_BRSreg_est <- data.frame(i_MC = i_MC,
                                       model=model,
                                       n=n,
                                       alpha=alpha,
                                       beta=beta,
                                       risk_measure = risk_measure,
                                       type = "BRSreg_est",
                                       theta_index = (1:length(sum_BRSreg$coefficients[,1])) + length(sum_BRSreg$coefficients[,1]),
                                       value = sum_BRSreg$coefficients[,1])
          
          
          df_BRSreg_cov_asy <-  data.frame(i_MC = i_MC,
                                          model=model,
                                          n=n,
                                          alpha=alpha,
                                          beta=beta,
                                          risk_measure = risk_measure,
                                          type = "BRSreg_cov_est_asy",
                                          theta_index = 1:length(MTS::Vech(sum_BRSreg$cov.unscaled)),
                                          value = MTS::Vech(sum_BRSreg$cov.unscaled))
          ### End BRS-type MES regression
          
          
          
          df_theta_est <- data.frame(i_MC = i_MC,
                                     model=model,
                                     n=n,
                                     alpha=alpha,
                                     beta=beta,
                                     risk_measure = risk_measure,
                                     type = "param_est",
                                     theta_index = 1:length(est_obj$theta),
                                     value = est_obj$theta)
          
          df_theta_true <- data.frame(i_MC = i_MC,
                                      model=model,
                                      n=n,
                                      alpha=alpha,
                                      beta=beta,
                                      risk_measure = risk_measure,
                                      type = "true_value",
                                      theta_index = 1:length(theta0),
                                      value = theta0)
          
          
          # Covariance estimation
          sum_obj_asy <- summary(est_obj)
          
          df_theta_cov_asy <-  data.frame(i_MC = i_MC,
                                          model=model,
                                          n=n,
                                          alpha=alpha,
                                          beta=beta,
                                          risk_measure = risk_measure,
                                          type = "cov_est_asy",
                                          theta_index = 1:length(MTS::Vech(sum_obj_asy$cov)),
                                          value = MTS::Vech(sum_obj_asy$cov))
          
          res_df <- rbind(res_df,
                          df_theta_est,
                          df_theta_true,
                          df_theta_cov_asy,
                          df_BRSreg_est,
                          df_BRSreg_cov_asy)
        }
      }
    }
  }
  res_df
}
stopCluster(cl)
end.time <- Sys.time()
(run.time <- end.time-start.time)

head(res_df_MC)
saveRDS(res_df_MC, file = "simulations/data/sim_crossDGP_20230525.rds")


