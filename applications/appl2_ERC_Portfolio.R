library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(tsibble)
library(rmgarch)
library(mvtnorm)
library(cubature)
library(lubridate)
library(doParallel)
library(xts)
library(abind)
library(MTS)
library(SystemicRisk)


# Settings
start_date <- "2006-12-25" # Monday of the last trading week of 2006
asset_names <- c("BA", "JPM", "LLY", "MSFT", "XOM")
covariate_names <- c("Intercept", "VIX")
m <- length(asset_names)


# Load VIX and transform to a weekly variable
data_VIX <- readRDS(file = "applications/data/data_VIX.rds") %>%
  mutate(Week=yearweek(Date),
         Month=yearmonth(Date)) %>%
  group_by(Week) %>%
  summarize(Date=min(Date), VIX=mean(VIX, na.rm=TRUE)) 


# Load portfolio assets and transform to a weekly variable
data_Assets_Portfolio <- readRDS(file = "applications/data/data_Assets.rds") %>%
  pivot_wider(names_from = Asset, values_from = NegReturn) %>%
  dplyr::select(Date, any_of(asset_names)) %>%
  mutate(Week=yearweek(Date),
         Month=yearmonth(Date)) %>%
  group_by(Week) %>%
  summarize(Date=min(Date), across(where(is.numeric), sum, na.rm=TRUE))


# Select covariates and assets for the portfolio construction
data_mod <- data_Assets_Portfolio %>%
  full_join(data_VIX, by=c("Week")) %>%
  dplyr::select(Week,
                any_of(asset_names),
                any_of(covariate_names)) %>%
  mutate(Intercept=1,
         Date=as_date(Week)) %>% # get Mondays!
  arrange(Date) %>%
  na.omit()


# Dates for refitting (indicating the ENDs of each estimation period in the rolling window)
refit_dates <- data_mod %>%
  dplyr::filter(Date >= start_date) %>%
  pull(Date)

M_dates <- length(refit_dates)



################################################################################
### GMVP Portfolio: DCC GARCH
################################################################################

# DCC specifications
xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)), 
                   variance.model = list(garchOrder = c(1,1), model = 'sGARCH'), 
                   distribution.model = 'std')
uspec = multispec(replicate(length(asset_names), xspec))
spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvt')

# DCC estimation
roll1 <- dccroll(spec = spec1, 
                 data = as.xts(data_Assets_Portfolio %>% dplyr::select(Date, asset_names)),
                 forecast.length = M_dates,
                 refit.every = 22,
                 fit.control = list(eval.se = TRUE))
Sigma_FC <- rcov(roll1)


# GMVP portfolio weights
iota <- rep(1,length(asset_names))
w_GMVPc_tbl <- tibble()
for (j in 1:(M_dates-1)){
  Sigma_inv <- solve(Sigma_FC[,,j])  
  w_GMVPc <- Sigma_inv %*% iota / as.numeric(iota %*% Sigma_inv %*% iota)
  w_GMVPc_tbl <- bind_rows(w_GMVPc_tbl,
                           as_tibble(t(w_GMVPc)) %>% 
                             mutate(Date = refit_dates[j+1])) # Shift by one week (unit) to match with the realizations!!!
}

# Convert to long data frame
res_appl_GMVP <- pivot_longer(w_GMVPc_tbl, cols=-Date, names_to="Asset", values_to="w_GMVP_DCC") %>%
  group_by(Date) %>%
  mutate(w_GMVP_DCC_normalized = pmin(1,pmax(0,w_GMVP_DCC,0)),
         w_GMVP_DCC_normalized = w_GMVP_DCC_normalized/sum(w_GMVP_DCC_normalized)) 





################################################################################
### ERC Portfolio:   Loop over all dates in parallel
################################################################################
core.max <- 80
cl <- makeCluster(min(parallel::detectCores()-1, M_dates, core.max) )
registerDoParallel(cl)
start.time <- Sys.time()
res_appl_ERC <- foreach(
  i_date = 1:M_dates,
  .combine=rbind,
  .packages=c("dplyr", "tibble", "tidyr", "MTS", "MASS", "mvtnorm", "cubature", "lubridate", "abind", "tsibble", "SystemicRisk"),
  .errorhandling="pass"
)%dopar%{
  
  set.seed(i_date) # set seed for reproducibility
  
  ###### Settings!
  beta <- 0.975
  max_counter <- 250
  eps <- 0.01
  
  date <- refit_dates[i_date]
  data_fit <- data_mod %>%
    dplyr::filter(Date <= date) %>%
    dplyr::select(-c("Week"))
  
  # extract matrix of asset returns
  asset_matrix <- data_fit %>%
    as_tibble() %>%
    dplyr::select(all_of(asset_names)) %>%
    as.matrix()
  
  
  # Initial weights!
  w <- rep(1/m,m)
  w_matrix <- matrix(NA,max_counter+1,m)
  MES_matrix <- matrix(NA,max_counter,m)
  MES_rel <- c(0,1) # Initialize the while loop
  w_matrix[1,] <- w
  
  # Initial starting values:
  theta0_list <- sapply(asset_names, function(x) NULL) # List with NULL entries
  
  counter <- 0
  while (counter <= max_counter-1 & (max(MES_rel) > 1/m+eps | min(MES_rel) < 1/m-eps) ){
    counter <- counter + 1
    # Filter with new, flexible weights!
    data_fit <- data_fit %>%
      mutate(Portfolio = as.numeric(asset_matrix %*% w)) %>%
      as_tsibble()
    
    # Loop over all asset (names) in the portfolio
    SRM_fits <- list()
    for (asset in asset_names) {
      SRM_fits[[asset]] <- SRM(data=data_fit %>%
                                 dplyr::select(-asset_names[!asset_names %in% asset]) %>%
                                 rename(Date=Date, x=Portfolio, y=!!asset) %>%
                                 as_tsibble(),
                               model = "joint_linear",
                               risk_measure="MES",
                               beta=beta,
                               theta0=theta0_list[[asset]],
                               optim_replications=c(1,3))
      # New starting values
      theta0_list[[asset]] <- SRM_fits[[asset]]$theta
    }
    
    # Vector of MES contributions
    MES_vec <- lapply(SRM_fits, FUN=function(SRM_obj) forecast(SRM_obj)[2]) %>%
      unlist() * w
    
    # Provide a lower bound for stability of the algorithm
    if (any(MES_vec < mean(pmax(0,MES_vec))/100)){
      warning("Hit the lower bound!")
      MES_vec <- pmax(MES_vec, mean(pmax(0,MES_vec))/100)
    }

    MES_matrix[counter,] <- MES_vec
    MES_rel <- MES_vec/sum(MES_vec)
    
    # Weights
    w_old <- w
    w_rel <- w_old/MES_vec
    w_rel <- w_rel/sum(w_rel)
    
    # Average between a new/MES-relative weight and the old one to converge "slower" 
    a <- 0.5
    w <- a*w_rel + (1-a)*w_old
    
    w <- w/sum(w)
    w_matrix[counter+1,] <- w
  }
  
  
  # Select the best weight in terms of equalizing MES
  MES_rel_matrix <- MES_matrix/rowSums(MES_matrix)
  index_best <- which.min(rowMeans((MES_rel_matrix - 1/m)^2))
  
  w_choice <- w_matrix[index_best,]
  MES_vec <- MES_matrix[index_best,]
  MES_EW_vec <- MES_matrix[1,]
  
  df <- tibble(Date=date, Asset = asset_names, w=w_choice, MES=MES_vec, MES_EW=MES_EW_vec, iterations=counter, index_best=index_best)
  df
}

stopCluster(cl)
end.time <- Sys.time()
(run.time <- end.time-start.time)

################################################################################
### End parallel loop
################################################################################



saveRDS(data_mod, file = "applications/data/data_mod_appl_ERC_20230802.rds")
saveRDS(res_appl_ERC, file = "applications/data/res_appl_ERC_20230802.rds")
saveRDS(res_appl_GMVP, file = "applications/data/res_appl_GMVP_20230802.rds")





