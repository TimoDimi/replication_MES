library(rmgarch)
library(mvtnorm)
library(cubature)

#============================================================
# DCC_to_CoVaR: Function to compute non-parametric VaR, MES and CoVaR forecasts for DCC-GARCH models
# ===========================================================

# Inputs:
# Sigma_FC     = forecasted (2 x 2) conditional 'square-root' decomposition of covariance matrix (for time n+1)
# nu           = scalar degrees of freedom of multivariate t-distribution (= NULL, if data_IS is supplied)
# data_IS      = a data.frame containing the in-sample values of x,y, and the standardized in-sample residuals (= NULL, if nu is supplied)
# alpha        = risk/probability level for CoVaR
# beta         = risk/probability levels for VaR and MES

DCC_to_CoVaR <- function(Sigma_FC, nu, data_IS, alpha, beta){
  
  if( is.null(data_IS) ){
    ### (A) Forecasts for t-distributed innovations
    # 1. Calculate VaR, ES and MES of eps
    VaR.hat.b <- MES.hat.b <- ES.hat.b <- numeric(length(beta))   # Note that CoVaR depends on $beta$ through the beta-quantile of eps_X
    sd.wt <- sqrt(nu / (nu-2) )
    ScaleMatrix_tdist <- (nu - 2) / nu * matrix(c(1, 0, 0, 1), nrow=2) # Correct for scale != Variance
    f.MES <- function(x, nu) { x[1] * dmvt(x, sigma = ScaleMatrix_tdist, df = nu, log = FALSE) } # "x" is vector
    
    for(b in 1:length(beta)){
      VaR.hat.b[b] <- qt(p=beta[b], df = nu) / sd.wt
      MES.hat.b[b] <- 1/(1-beta[b]) * adaptIntegrate(f.MES, lowerLimit = c(-Inf, VaR.hat.b[b]), upperLimit = c(Inf, Inf), nu)$integral
    }
    ES.hat.b <- 1/sd.wt * (dt(qt(beta, df=nu) , df=nu) / (1 - beta)) * ((nu + (qt(beta, df=nu))^2) / (nu - 1))  # see [MFE15, Example 2.15]
    
    # 2. Compute conditional VaR, CoVaR and MES
    # The following formula ONLY apply to the Cholesky decomposition! (This could be implemented easier given we have Sigma_FC already, right?)
    H_FC <- Sigma_FC %*% t(Sigma_FC) # Covariance forecast
    rho.np1         <- H_FC[2,1] / sqrt( H_FC[1,1] * H_FC[2,2] )
    theta.hat.n.b.1 <- sqrt(H_FC[2,2]) * sqrt(1 - rho.np1^2) * MES.hat.b      # MES part
    theta.hat.n.b.2 <- sqrt(H_FC[2,2]) * rho.np1 *  ES.hat.b                  # ES part
    MES_FC     <- theta.hat.n.b.1 + theta.hat.n.b.2                      # this is the MES forecast
    VaR_FC     <- sqrt(H_FC[1,1]) * VaR.hat.b                            # this is the VaR forecast
    CoVaR_FC   <- numeric( length(beta))
    for(b in 1:length(beta)){
      CoVaR_FC[b] <- CoVaR.true.t(alpha, beta[b], nu, H_FC)
    }
    
  } else {
    ### (B) Forecasts based on nonparametrically computing the (systemic) risk measures from past DCC-innovations
    # 1. Compute VaR forecasts
    n <- nrow(data_IS)
    
    # u_hat := Sigma_Forecast * eps_InSample, a combination of the in sample residuals and forecastsed covariance matrix
    u_hat <- Sigma_FC %*% (data_IS %>% dplyr::select(eps_x, eps_y) %>% as.matrix() %>% t())
    data_IS$u_hat_x <- u_hat[1,]
    data_IS$u_hat_y <- u_hat[2,]
    
    VaR_FC <- quantile(data_IS$u_hat_x, probs = beta)
    
    # 2. Compute conditional risk measures
    MES_FC <- CoVaR_FC <- numeric( length(beta) ) # Note that CoVaR depends on $beta$ through the beta-quantile of eps_X
    
    for(b in 1:length(beta)){
      u_hat_Y_VaRviolation <- data_IS %>%
        dplyr::filter(u_hat_x > VaR_FC[b]) %>%
        dplyr::pull(u_hat_y)
      CoVaR_FC[b]  <- quantile(u_hat_Y_VaRviolation, probs=alpha)
      MES_FC[b]    <-     mean(u_hat_Y_VaRviolation)
    }
  }
  
  return(list(VaR = VaR_FC, MES = MES_FC, CoVaR = CoVaR_FC))
}



# ===========================================================
# Function for CoVaR calculation for bivariate, zero-mean t-distribution
# ===========================================================

# Inputs:
# alpha = confidence level of CoVaR
# beta  = confidence level of VaR
# nu    = degrees of freedom
# H     = variance-covariance matrix

CoVaR.true.t <- function(alpha, beta, nu, H){
  root.CoVaR.t <- function(x){
    sigma <-  (nu - 2) / nu * H
    sd.wt <- sqrt(nu / (nu-2) )
    VaR   <- qt(p = beta, df = nu) / sd.wt * sqrt( H[1,1] )
    prob  <- mvtnorm::pmvt(lower=c(VaR, x), upper=c(Inf, Inf), df=nu, sigma=sigma)
    
    return( prob - (1-alpha) * (1-beta) )
  }
  uniroot(root.CoVaR.t, interval = c(0, 10), extendInt = "yes")$root
}


# ===========================================================
# Function for MES calculation for bivariate t-distribution
# ===========================================================

# Inputs:
# nu   = degrees of freedom
# rho  = correlation
# beta = risk level

MES.true.t <- function(nu, rho, beta){
  if(rho==0){
    MES <- 0
  }
  else{
    sd.wt <- sqrt(nu / (nu-2) )
    f.MES <- function(x, nu, rho) { x[2] * mvtnorm::dmvt(x, sigma = (nu - 2) / nu * matrix(c(1, rho, rho, 1), nrow=2), df = nu, log = FALSE) } # "x" is vector
    VaR   <- qt(p = beta, df = nu) / sd.wt
    MES   <- 1/(1 - beta) * cubature::adaptIntegrate(f.MES, lowerLimit = c(VaR, -Inf), upperLimit = c(Inf, Inf), nu, rho)$integral
  }
  return( MES )
}






