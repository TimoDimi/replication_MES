
# This function generates the data according to the DGP used in Section 4
sim_SRM_DGP <- function(n, gamma, Sigma, nu=8, phi=c(0.5,0.8), beta, alpha){
  logZ1 <- 0.5 + 0.3*arima.sim(list(order=c(1,0,0), ar=phi[1]), n=n)
  Z2 <- arima.sim(list(order=c(1,0,0), ar=phi[2]), n=n)
  z <- rbind(1, as.numeric(exp(logZ1)), as.numeric(Z2)) %>% t()

  cov_matrix <- Sigma %*% t(Sigma)

  if (is.finite(nu)) {scale_matrix <- (nu - 2) / nu * cov_matrix}
  else {scale_matrix <- cov_matrix}

  eps <- cbind(gamma[4] + gamma[5]*z[,2] + gamma[6]*z[,3],
               gamma[4] + gamma[5]*z[,2] + gamma[6]*z[,3]) *
    mvtnorm::rmvt(n=n, sigma=scale_matrix, df = nu)
  eps_df <- data.frame(eps_x=eps[,1], eps_y=eps[,2])

  # Lag the covariates by one time point for the new implementation in the "Systemic Risk" package
  xy <- as.numeric(z %*% as.matrix(gamma[1:3])) + eps
  x <- c(NA,xy[,1])
  y <- c(NA,xy[,2])
  z <- rbind(z, NA)
  data <- tibble::tibble(x=x, y=y, z1=z[,1], z2=z[,2], z3=z[,3])

  # True VaR parameter values
  VaR_tdist <- qt(beta, df=nu) * sqrt(scale_matrix[1,1])
  theta_v_true <- c(gamma[1] + gamma[4]*VaR_tdist,
                    gamma[2] + gamma[5]*VaR_tdist,
                    gamma[3] + gamma[6]*VaR_tdist)

  # True CoVaR parameter values
  CoVaR_tdist <- CoVaR_tdist(alpha=alpha, beta=beta, nu=nu, H=cov_matrix)
  theta_c_true <- c(gamma[1] + gamma[4]*CoVaR_tdist,
                    gamma[2] + gamma[5]*CoVaR_tdist,
                    gamma[3] + gamma[6]*CoVaR_tdist)

  # True MES parameter values
  MES_tdist <- MES_tdist(beta=beta, nu=nu, H=cov_matrix)
  theta_m_true <- c(gamma[1] + gamma[4]*MES_tdist,
                    gamma[2] + gamma[5]*MES_tdist,
                    gamma[3] + gamma[6]*MES_tdist)

  theta_true <- list(VaR=theta_v_true, CoVaR=theta_c_true, MES=theta_m_true)
  return(list(data=data, theta_true=theta_true))
}
