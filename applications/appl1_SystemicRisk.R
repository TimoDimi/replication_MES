rm(list = ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library(rmgarch)
library(mvtnorm)
library(cubature)
library(lubridate)
library(tsibble)
library(SystemicRisk)
library(rlist)
library(pracma)
library(patchwork)


# Load function to convert DCC modes to systemic risk forecasts
source("applications/DCC_to_SystemicRisk.R")


####################################################################################
###
###     NOTICE:
###     The following code requires the file "data_bubbles.rds" which we cannot make
###     available publicly as we do not have the licence to publish the raw data
###
####################################################################################

# Load data sets
data_Bubbles <- readRDS(file = "applications/data/data_bubbles.rds")

data_Spreads <- readRDS(file = "applications/data/data_Spreads.rds")%>%
  arrange(Date)

data_VIX <- readRDS(file = "applications/data/data_VIX_RR.rds") 

data_Assets_long <- readRDS(file = "applications/data/data_Assets_RR.rds") 
# data_Assets_long <- readRDS(file = "applications/data/data_Assets.rds") 

data_Assets_wide <- data_Assets_long %>% 
  pivot_wider(names_from = Asset, values_from = NegReturn) %>%
  mutate(GS = case_when(Date <= "1999-05-05" ~ -1000, .default = GS)) # A placeholder before GS started trading


data_joint <- full_join(data_Spreads, data_Assets_wide, by="Date") %>%
  full_join(data_VIX, by="Date") %>%
  full_join(data_Bubbles, by="Date") %>%
  dplyr::select(Date, ChangeSpread, TEDSpread, VIX, st_boom_bsadf, st_bust_bsadf,
         BAC, C, JPM, SPF, MS, GS, WFC, BK, STT) %>%
  dplyr::filter(Date >= date("1993-01-01") & Date <= date(" 2015-12-31")) %>%
  arrange(Date) %>%
  fill(c(st_boom_bsadf, st_bust_bsadf), .direction = "down") %>%
  dplyr::filter(Date >= date("1993-05-05")) %>%
  na.omit()


x_set <-  c("BAC", "C", "JPM", "MS", "GS", "WFC", "BK", "STT")
x_name <- c("BAC", "C", "JPM", "MS", "GS", "WFC", "BK", "STT")
beta <- 0.95

MES_obj_list <- list()
BRSreg_obj_list <- list()
DCC_df_list <- list()

MES_dfsummary_list <- list()
BRSreg_dfsummary_list <- list()

BRSreg_df_list <- list()

for (i in 1:length(x_set)){
  # Transform to a tsibble, where x is SPF, and y the individual banks
  df <- data_joint %>%
    reframe(Date=Date,
            y=!!sym(x_set[i]),
            x=SPF,
            Intercept=1,
            ChangeSpread = ChangeSpread,
            TEDSpread = TEDSpread,
            VIX = VIX,
            ST_boom=lead(st_boom_bsadf),  
            ST_bust=lead(st_bust_bsadf)
    ) %>%
    na.omit() %>%
    as_tsibble(index=Date)
  # Note that we need "lead" for the boom/bust indices as the joint linear (SRM) model automatically shifts the 
  # variables in the SystemicRisk package, but BRS(2020) use contemporaneous indices. 
  
  
  # Filter non-trading days of GS
  if (x_set[i] == "GS"){ df <- df %>% filter(y!=-1000) }
  
  set.seed(2023)
  MESreg_obj <- SRM(data=df,
                    model="joint_linear",
                    risk_measure="MES",
                    beta=beta,
                    optim_replications=c(3,10))
  
  MES_obj_list[[x_name[i]]] <- MESreg_obj
  
  MES_summary <- summary(MESreg_obj)
  MES_dfsummary_list[[x_name[i]]] <- full_join(as_tibble(MES_summary$coef_mat_VaR, rownames="Covariate") %>%
                                                 rename(Estimate=2, Error=3, tval=4, pval=5),
                                               as_tibble(MES_summary$coef_mat_risk_measure, rownames="Covariate") %>%
                                                 rename(Estimate=2, Error=3, tval=4, pval=5), 
                                               by="Covariate",
                                               suffix = c("_VaR", "_MES"))
  
  
  #### Brunnermeier et al (RFS, 2020) "MES" regression
  S <- 250
    
  data_BRSreg <- df %>%
    mutate(index_Date=1:n(),
           Qx_rolling = slider::slide_dbl(x, quantile, .before = S, .after = 0, beta, na.rm=TRUE),
           y_x_viol = ifelse((x >= Qx_rolling),y,NA),
           y_MESrolling = slider::slide_dbl(y_x_viol, mean, .before = S, .after = 0, na.rm=TRUE),
           ChangeSpread=lag(ChangeSpread),
           TEDSpread=lag(TEDSpread),
           VIX=lag(VIX),
           ST_boom=lag(ST_boom),
           ST_bust=lag(ST_bust)) %>%
    fill(y_MESrolling, .direction = 'down') 
  # Note: Now we lag the covariates as we the "lm" function below does not use lags automatically (in contast to the SRM function)!
  # We also lag the ST_boom and ST_bust variables as they were "leaded" in constructing df above.
           
  
  lmBRSreg <- lm(y_MESrolling ~ ChangeSpread + TEDSpread + VIX  + ST_boom + ST_bust,
                 data_BRSreg)
  BRSreg_obj_list[[x_name[i]]] <- lmBRSreg
  
  sum_BRSreg <- summary(lmBRSreg)
  BRSreg_dfsummary_list[[x_name[i]]] <- as_tibble(sum_BRSreg$coefficients, rownames="Covariate") %>%
    rename(Estimate=2, Error=3, tval=4, pval=5)
  
  BRSreg_df_list[[x_name[i]]] <- data_BRSreg
  
  #### DCC GARCH model
  # DCC GARCH specifications
  GARCH11_norm_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean=FALSE),
                                  variance.model = list(garchOrder = c(1,1), model = "sGARCH"),
                                  distribution.model = "norm")
  DCC_spec <- dccspec(uspec = multispec(replicate(2, GARCH11_norm_spec)), dccOrder = c(1,1), distribution = "mvnorm")
  
  # DCC fit
  DCC_fit <- dccfit(DCC_spec, 
                    data = as_tibble(df) %>% select(x,y), 
                    fit.control = list(eval.se = FALSE, stationarity = TRUE, scale = FALSE))
  
  # MES from DCC
  H <- rcov(DCC_fit)
  Sigma <- sapply(1:dim(H)[3], function(i) pracma::sqrtm(H[,,i])$B, simplify = "array")
  risk_measure_matrix <- sapply(1:dim(H)[3], 
                                function(i) DCC_to_CoVaR(Sigma_FC=Sigma[,,i], nu=10^4, data_IS=NULL, alpha=0.5, beta=beta),
                                simplify = "array")
  
  # Collect (in-sample) model predictions
  DCC_df_list[[x_name[i]]] <- df %>%
    mutate(VaR = risk_measure_matrix[1,] %>% as.numeric(),
           MES = risk_measure_matrix[2,] %>% as.numeric())
  
}






# Look at the parameter estimates
MES_dfsummary_list
BRSreg_dfsummary_list


# Print main(BAC, C, JPM) results to a LaTeX table
full_join(bind_rows(MES_dfsummary_list[c("BAC", "C", "JPM")], .id = "Asset"),
          bind_rows(BRSreg_dfsummary_list[c("BAC", "C", "JPM")], .id = "Asset") %>% 
            mutate(Covariate = ifelse(Covariate=="(Intercept)", "Intercept", Covariate)),
          by=c("Asset", "Covariate")) %>%
  mutate(empty1=NA, empty2=NA, empty3=NA, empty4=NA) %>%
  dplyr::select(Asset, empty1, Covariate, 
                empty2, Estimate_VaR, Error_VaR,  pval_VaR,
                empty3, Estimate_MES, Error_MES, pval_MES,
                empty4, Estimate, Error, pval) %>%
  xtable::xtable(digits=c(0, 0, 0, 0, 0, c(3,3,3), 0, c(3,3,3), 0, c(3,3,3))) %>%
  print(file="applications/output/MES_AB_regressions_RR.txt", include.rownames=FALSE, booktabs=TRUE)


# Print supplemental ("BK", "GS", "MS", "STT", "WFC") results to a LaTeX table
full_join(bind_rows(MES_dfsummary_list[c("BK", "GS", "MS", "STT", "WFC")], .id = "Asset"),
          bind_rows(BRSreg_dfsummary_list[c("BK", "GS", "MS", "STT", "WFC")], .id = "Asset") %>% 
            mutate(Covariate = ifelse(Covariate=="(Intercept)", "Intercept", Covariate)),
          by=c("Asset", "Covariate")) %>%
  mutate(empty1=NA, empty2=NA, empty3=NA, empty4=NA) %>%
  dplyr::select(Asset, empty1, Covariate, 
                empty2, Estimate_VaR, Error_VaR,  pval_VaR,
                empty3, Estimate_MES, Error_MES, pval_MES,
                empty4, Estimate, Error, pval) %>%
  xtable::xtable(digits=c(0, 0, 0, 0, 0, c(3,3,3), 0, c(3,3,3), 0, c(3,3,3))) %>%
  print(file="applications/output/MES_AB_regressions_AddBanks_RR.txt", include.rownames=FALSE, booktabs=TRUE)





### Plot predictions exemplary for BAC
asset_plot <- "BAC"

# Generate a data frame for plotting
df_MESBRSreg_plot <- MES_obj_list[[asset_plot]]$data %>% 
  na.omit() %>%
  dplyr::select(Date, x, y, VaR, risk_measure) %>%
  rename(MES_reg = risk_measure) %>%
  mutate(VaR_viol = x >= VaR,
         BRS_reg = BRSreg_obj_list[[asset_plot]] %>% predict()) %>%
  pivot_longer(cols=c(MES_reg, BRS_reg), names_to="Method") %>%
  na.omit() %>%
  dplyr::rename("VaR Exceedance" = "VaR_viol")

# Assign factor for nicer names and ordering
df_MESBRSreg_plot$Method <- factor(df_MESBRSreg_plot$Method,
                                   levels=c("MES_reg", "BRS_reg"),
                                   labels=c("MES Regression", "BRS Regression"))


# Plot the BRS and MES regression (in-sample) predictions
ggplot(df_MESBRSreg_plot) + 
  geom_point(data=df_MESBRSreg_plot %>% arrange(`VaR Exceedance`), 
             aes(x=Date, y=y, col=`VaR Exceedance`), size=1) +
  ggplot2::scale_colour_manual(values = c("grey", "black")) +
  ggnewscale::new_scale_color() + # For two color scales!!!
  geom_line(aes(x=Date, y=value, col=Method), linewidth=0.5) + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  ylab("Negative BAC Return") +
  coord_cartesian(ylim=c(-10,30))

ggsave("applications/output/MES_BRS_regressions.pdf", width = 8, height = 4.5, units = "in")






### Plot predictions including GARCH-DCC
asset_plot <- "BAC"

# Generate a data frame for plotting
df_MESBRSreg_plot <- MES_obj_list[[asset_plot]]$data %>% 
  dplyr::select(Date, x, y, VaR, risk_measure) %>%
  rename(MES_reg = risk_measure) %>%
  na.omit() %>%
  mutate(VaR_viol = (x >= VaR),
         BRS_reg = BRSreg_obj_list[[asset_plot]] %>% predict()) %>%
  full_join(DCC_df_list[[asset_plot]] %>% select(Date,y,x,MES) %>% rename(MES_DCC=MES), 
            by=c("Date","x","y")) %>%
  pivot_longer(cols=c(MES_reg, BRS_reg, MES_DCC), names_to="Method") %>%
  na.omit() %>%
  dplyr::rename("VaR Exceedance" = "VaR_viol")

# Assign factor for nicer names and ordering
df_MESBRSreg_plot$Method <- factor(df_MESBRSreg_plot$Method,
                                   levels=c("MES_reg", "BRS_reg", "MES_DCC"),
                                   labels=c("MES Regression", "BRS Regression", "DCC-GARCH"))


# Plot the BRS and MES regression (in-sample) predictions
ggplot(df_MESBRSreg_plot) + 
  geom_point(data=df_MESBRSreg_plot %>% arrange(`VaR Exceedance`), 
             aes(x=Date, y=y, col=`VaR Exceedance`), size=1) +
  ggplot2::scale_colour_manual(values = c("grey", "black")) +
  ggnewscale::new_scale_color() + # For two color scales!!!
  geom_line(aes(x=Date, y=value, col=Method), linewidth=0.5) + 
  theme_bw() +
  theme(legend.position = "bottom") + 
  ylab("Negative BAC Return") +
  coord_cartesian(ylim=c(-10,30))

ggsave("applications/output/MES_BRS_DCC_regressions.pdf", width = 8, height = 4.5, units = "in")








###### In-Sample Model Diagnostics

beta <- 0.95
stock_choice <- "BAC"

### Extract data and residuals 

# MES regression generalized residuals
df_ResidDiag_MESreg <- MES_obj_list[[stock_choice]]$data %>%
  select(Date,x,y,VaR,risk_measure,ST_boom, ST_bust) %>%
  mutate(ST_boom=lag(ST_boom),    # Lag covariates by one time period for MES reg as MES reg regresses on previous rows in the data frame! (see above)
         ST_bust=lag(ST_bust)) %>%
  tail(-1) %>%
  rename(MES=risk_measure) %>%
  mutate(GenResid_VaR = (x<=VaR) - beta,
         GenResid_MES_zeros = (x>VaR)*(MES-y),
         GenResid_MES = na_if(GenResid_MES_zeros, (x>VaR)),
         model="MESreg") %>%
  as_tibble


# Extract BRS regression data
df_ResidDiag_BRS <- BRSreg_df_list[[stock_choice]]  %>%
  tail(-1) %>%
  mutate(BRSreg = predict(BRSreg_obj_list[[stock_choice]])) %>%
  select(Date,x,y,y_x_viol,Qx_rolling, BRSreg,ST_boom, ST_bust) %>%
  rename(MES=BRSreg) %>%
  as_tibble() %>% 
  left_join(df_ResidDiag_MESreg %>% select(Date, VaR, GenResid_VaR), by="Date") %>% # Include QuantReg infomration from SRM to have a conditional quantile!
  # Option 1: Use rolling Qx as VaR 
  mutate(GenResid_VaR = (x<=Qx_rolling) - beta,
         GenResid_MES_zeros = (x>Qx_rolling)*(MES-y),
         GenResid_MES = na_if(GenResid_MES_zeros, (x>Qx_rolling)),
         model="BRSreg")


# Second option with SRM VaR as BRS has no VaR for model diagnostics!
df_ResidDiag_BRS_QR <- BRSreg_df_list[[stock_choice]]  %>%
  tail(-1) %>%
  mutate(BRSreg = predict(BRSreg_obj_list[[stock_choice]])) %>%
  select(Date,x,y,y_x_viol,Qx_rolling, BRSreg,ST_boom, ST_bust) %>%
  rename(MES=BRSreg) %>%
  as_tibble() %>% 
  left_join(df_ResidDiag_MESreg %>% select(Date, VaR, GenResid_VaR), by="Date") %>% # Include QuantReg infomration from SRM to have a conditional quantile!
  # Option 2: We use SRM VaR as BRS has no VaR for model diagnostics!
  mutate(GenResid_MES_zeros = (x>VaR)*(MES-y),
         GenResid_MES = na_if(GenResid_MES_zeros, (x>VaR)),
         model="BRSreg_QR")



# DCC generalized residuals
df_ResidDiag_DCC <- DCC_df_list[[stock_choice]] %>%
  tail(-1) %>%
  select(Date,x,y,VaR,MES,ST_boom, ST_bust) %>%
  mutate(GenResid_VaR = (x<=VaR) - beta,
         GenResid_MES_zeros = (x>VaR)*(MES-y),
         GenResid_MES = na_if(GenResid_MES_zeros, (x>VaR)),
         model="DCC") %>%
  as_tibble()



### Join MESreg, DCC and BRS data frames  
df_ResidDiag <- bind_rows(df_ResidDiag_MESreg, df_ResidDiag_DCC, df_ResidDiag_BRS, df_ResidDiag_BRS_QR) %>%
  arrange(Date, model)


# Edit the data frame such that it fits in our plotting routine.
df_joint_plot <- bind_rows(df_ResidDiag %>% 
                             select(Date, model, MES, GenResid_MES) %>%
                             filter(model %in% c("BRSreg", "MESreg")) %>%
                             rename(GenResid = GenResid_MES, Prediction = MES) %>%
                             mutate(model=case_when(model=="BRSreg" ~ "MES_BRSreg",
                                                    model=="MESreg" ~ "MES_MESreg")),
                           df_ResidDiag %>% 
                             select(Date, model, VaR, GenResid_VaR) %>%
                             filter(model %in% c("MESreg")) %>%
                             rename(GenResid = GenResid_VaR, Prediction = VaR) %>%
                             mutate(model=case_when( model=="MESreg" ~ "VaR_MESreg"))
)  %>% 
  tidyr::drop_na(GenResid)




# Tune facet labels
df_joint_plot$model <- factor(df_joint_plot$model,
                              levels=c("VaR_MESreg", "MES_MESreg", "MES_BRSreg"),
                              labels=c("VaR Regression", "MES Regression", "BRS Regression"))

p_Date <- ggplot(df_joint_plot, aes(x=Date, y=GenResid)) +
  facet_wrap(~model, scales="free") +
  geom_point() + 
  geom_smooth(method="loess", fill="blue") + 
  geom_hline(yintercept=0) + 
  theme_bw() +
  # theme(aspect.ratio=1) +
  ylab("Generalized Residuals") +
  xlab("Date") 

ggsave("applications/output/SysRisk_Diagnostics_Date.pdf", 
       p_Date, width = 9, height = 3, units = "in")


p_FC <- ggplot(df_joint_plot, aes(x=Prediction, y=GenResid)) +
  facet_wrap(~model, scales="free") +
  geom_point() + 
  geom_smooth(method="loess", fill="blue") + 
  geom_hline(yintercept=0) + 
  theme_bw() +
  # theme(aspect.ratio=1) +
  ylab("Generalized Residuals") +
  xlab("Model Predictions") 

ggsave("applications/output/SysRisk_Diagnostics_Predictions.pdf", 
       p_FC, width = 9, height = 3, units = "in")






###### Compare the two quantile regression evaluations for the BRS regression
df_plot_BRScomp <- df_ResidDiag %>% 
  select(Date, model, MES, GenResid_MES) %>%
  filter(model %in% c("BRSreg", "BRSreg_QR")) %>%
  rename(GenResid = GenResid_MES, Prediction = MES) %>%
  mutate(model=case_when(model=="BRSreg" ~ "MES_BRSreg_Qrolling",
                         model=="BRSreg_QR" ~ "MES_BRSreg_Qregression")) %>% 
  tidyr::drop_na(GenResid)



# Tune facet labels
df_plot_BRScomp$model <- factor(df_joint_plot$model,
                              levels=c("MES_BRSreg_Qrolling", "MES_BRSreg_Qregression"),
                              labels=c("BRS Regression with the rolling quantile", "BRS Regression with quantile regression" ))

p_Date_BRScomp  <- ggplot(df_plot_BRScomp, aes(x=Date, y=GenResid)) +
  facet_wrap(~model, scales="free") +
  geom_point() + 
  geom_smooth(method="loess", fill="blue") + 
  geom_hline(yintercept=0) + 
  theme_bw() +
  # theme(aspect.ratio=1) +
  ylab("Generalized Residuals") +
  xlab("Date") 

ggsave("applications/output/SysRisk_BRScomp_Diagnostics_Date.pdf", 
       p_Date_BRScomp, width = 9, height = 3, units = "in")


p_FC_BRScomp <- ggplot(df_plot_BRScomp, aes(x=Prediction, y=GenResid)) +
  facet_wrap(~model, scales="free") +
  geom_point() + 
  geom_smooth(method="loess", fill="blue") + 
  geom_hline(yintercept=0) + 
  theme_bw() +
  # theme(aspect.ratio=1) +
  ylab("Generalized Residuals") +
  xlab("Model Predictions") 

ggsave("applications/output/SysRisk_BRScomp_Diagnostics_Predictions.pdf", 
       p_FC_BRScomp, width = 9, height = 3, units = "in")









# Edit the data frame such that it fits in our plotting routine.
df_plot_DCC <- bind_rows(df_ResidDiag %>% 
                           select(Date, model, MES, GenResid_MES) %>%
                           filter(model %in% c("DCC", "MESreg")) %>%
                           rename(GenResid = GenResid_MES, Prediction = MES) %>%
                           mutate(model=case_when(model=="DCC" ~ "MES_DCC",
                                                  model=="MESreg" ~ "MES_MESreg")),
                         df_ResidDiag %>% 
                           select(Date, model, VaR, GenResid_VaR) %>%
                           filter(model %in% c("DCC", "MESreg")) %>%
                           rename(GenResid = GenResid_VaR, Prediction = VaR) %>%
                           mutate(model=case_when(model=="DCC" ~ "VaR_DCC",
                                                  model=="MESreg" ~ "VaR_MESreg"))
)  %>% 
  tidyr::drop_na(GenResid)




# Tune facet labels
df_plot_DCC$model <- factor(df_plot_DCC$model,
                              levels=c("VaR_MESreg", "MES_MESreg", "VaR_DCC", "MES_DCC"),
                              labels=c("VaR Regression", "MES Regression", "VaR DCC-GARCH", "MES DCC-GARCH"))

p_Date_DCC <- ggplot(df_plot_DCC, aes(x=Date, y=GenResid)) +
  facet_wrap(~model, scales="free", nrow=1) +
  geom_point() + 
  geom_smooth(method="loess", fill="blue") + 
  geom_hline(yintercept=0) + 
  theme_bw() +
  # theme(aspect.ratio=1) +
  ylab("Generalized Residuals") +
  xlab("Date") 

ggsave("applications/output/SysRisk_DCC_Diagnostics_Date.pdf", 
       p_Date_DCC, width = 12, height = 3, units = "in")


p_FC_DCC <- ggplot(df_plot_DCC, aes(x=Prediction, y=GenResid)) +
  facet_wrap(~model, scales="free", nrow=1) +
  geom_point() + 
  geom_smooth(method="loess", fill="blue") + 
  geom_hline(yintercept=0) + 
  theme_bw() +
  # theme(aspect.ratio=1) +
  ylab("Generalized Residuals") +
  xlab("Model Predictions") 

ggsave("applications/output/SysRisk_DCC_Diagnostics_Predictions.pdf", 
       p_FC_DCC, width = 12, height = 3, units = "in")


