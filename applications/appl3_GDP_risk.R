rm(list = ls())

library(dplyr)
library(tidyverse)
library(tsibble)
library(xts)
library(lubridate)
library(SystemicRisk)

# Load function to convert DCC modes to systemic risk forecasts
source("applications/DCC_to_SystemicRisk.R")


# Set some options
countries_list <-  c("DEU","FRA","GBR")
beta <- 0.9

# Load data frames
data_GaR <- readRDS("applications/data/data_GaR.rds")
data_GDPtot <- readRDS("applications/data/data_GaR_total.rds")
data_gFCI <- readRDS("applications/data/data_gFCI.rds") 


# Generate manual (individual country) growth rates such that we know the weights
data_GDP_hlp <- full_join(data_GaR %>% 
                            dplyr::filter(Date >= "1973-01-01",
                                          Country %in% countries_list),
                          data_GDPtot %>% 
                            dplyr::filter(Date >= "1973-01-01", 
                                          Country %in% countries_list) %>%
                            group_by(Date) %>%
                            mutate(GDPrel=GDPrel/sum(GDPrel)),
                          by=c("Date","Country")) %>%
  group_by(Country) %>%
  fill(GDPrel) %>%
  full_join(data_gFCI %>% dplyr::filter(Country %in% countries_list), by=c("Date","Country")) %>%
  ungroup()


 # Save the weights
data_GDP_weights <- data_GDP_hlp %>%
  dplyr::select(Date, Country, GDPrel) %>%
  rename(GDPweights = GDPrel)
  
# Display average weights:
data_GDP_weights %>% 
  group_by(Country) %>%
  summarize(mean(GDPweights))

# Generate data frame with absolute GDP growth rates
data_GaR_GDPrel_hlp <- data_GaR %>%
  filter(Country %in% c(countries_list)) %>%
  dplyr::select(Date, Country, GDP_Growth) %>% 
  inner_join(data_GDP_weights, by=c("Date", "Country")) %>% 
  na.omit() %>%
  mutate(GDP_loss_rel = -3*GDPweights*GDP_Growth,
         GDP_loss_abs = -GDP_Growth) %>%
  dplyr::select(Date, Country, GDP_loss_rel, GDP_loss_abs, GDPweights)


data_GaR_GDPrel <- bind_rows(data_GaR_GDPrel_hlp,
                             data_GaR_GDPrel_hlp %>%
                               group_by(Date) %>%
                               summarise(GDP_loss_abs=sum(GDP_loss_abs*GDPweights),
                                         Country="Region")) %>% 
  dplyr::select(Date, Country, GDP_loss_abs) %>%
  pivot_wider(names_from=Country, values_from=GDP_loss_abs)





# Filter and merge data that we need
data_fit <- full_join(data_GDP_hlp %>%
                        filter(Country=="DEU") %>% # This does not matter
                        dplyr::select(Date, ae_ew_fci),
                      data_GaR_GDPrel,
                      by="Date") %>%
  arrange(Date) %>%
  rename(FCI=ae_ew_fci) %>%
  na.omit() %>%
  dplyr::filter(Date <= date("2019-12-30"))



### In-sample fits
set.seed(2023)

SRM_est_list <- list()
SRM_df <- tibble()
for (country in countries_list){
  tsibble_fit <- data_fit %>% 
    mutate(Intercept=1, 
           Regionlag = Region) %>%
    dplyr::select(Date, Region, !!country, Intercept, FCI, Regionlag) %>%
    rename(Date=Date, x=Region, y=!!country) %>%
    as_tsibble() %>%
    dplyr::filter(Date <= date("2019-12-30"))
  
  SRM_est_list[[country]] <- SRM(data=tsibble_fit,
                                 model = "joint_linear",
                                 risk_measure="MES",
                                 beta=beta,
                                 optim_replications=c(3,10))
  
  SRM_sum <- summary(SRM_est_list[[country]])
  SRM_df <- bind_rows(SRM_df,
                      bind_rows(as_tibble(SRM_sum$coef_mat_VaR, rownames="Variable") %>% 
                                  mutate(Country=country, measure="VaR"), 
                                as_tibble(SRM_sum$coef_mat_risk_measure, rownames="Variable") %>% 
                                  mutate(Country=country, measure="MES")))
}


### Compare with a two-step ES regression
set.seed(2023)

# Unconstraint ES optimization
SRM_ES <- SRM(data=tsibble_fit %>% mutate(y=x),
              model = "joint_linear",
              risk_measure="MES",
              beta=0.9,
              optim_replications=c(3,10)) 

summary(SRM_ES)


# Constraint ES optimization

# Second step ES opjective function
ES_twostep_loss <- function(theta, df_SRM_ES){
  model <- theta[1] + lag(df_SRM_ES$FCI) * theta[2] + lag(df_SRM_ES$Regionlag) * theta[3]
  0.5 * (df_SRM_ES$x - model)^2 * (df_SRM_ES$x >= df_SRM_ES$VaR)
}

# Pertubed unconstraint QR estimates as starting values such that they are in the interior of the feasible parameter space
theta_init <- SRM_ES$theta[1:3] + c(0.05,0,0)  

# Linear ineqaulity constraints
ui_constOpt <- SRM_ES$data %>% as_tibble()  %>% select(Intercept, FCI, Regionlag) %>% as.matrix() %>% head(-1)
ci_constOpt <- tail(SRM_ES$data$VaR, -1)

# Constraint optimization
val_constOpt <-  constrOptim(theta=theta_init,
                             f=function(theta,...){
                               mean(ES_twostep_loss(theta,...), na.rm=TRUE)},
                             grad=NULL,
                             ui = ui_constOpt,
                             ci=ci_constOpt,
                             df_SRM_ES=SRM_ES$data)

ES_pred <- ui_constOpt %*% val_constOpt$par


# Print the MES coefficients
SRM_df %>%
  dplyr::filter(measure=="MES") %>%
  arrange(Variable) %>%
  mutate(CI_lower = Estimate - qnorm(0.95) * `Std. Error`,
         CI_upper = Estimate + qnorm(0.95) * `Std. Error`) 



# Joint plot of in-sample fits
df_ES <- SRM_ES$data %>%
  dplyr::mutate(VaR_violation=(x > VaR),
                Country="Joint Economic Region",
                ES=c(NA,ES_pred)) %>%
  dplyr::select(Date, Country, VaR_violation, x, VaR, risk_measure, ES) %>%
  dplyr::rename(ES_unconstr=risk_measure,
                return=x) %>%
  pivot_longer(cols=c(VaR,ES_unconstr,ES), names_to="risk_measure") %>%
  dplyr::rename("VaR Exceedance" = "VaR_violation") %>%
  as_tibble()


tol <- 1e-5 # Tolerance for rounding errors
df_MES <- sapply(SRM_est_list, 
                 function(SRM_obj){SRM_obj %>% 
                     .[["data"]] %>% 
                     dplyr::mutate(VaR_violation=(x > VaR - tol)) %>%
                     dplyr::select(Date, VaR_violation, y, risk_measure) %>%
                     dplyr::rename(value=risk_measure,
                                   return=y) %>%
                     mutate(risk_measure="MES") %>%
                     dplyr::rename("VaR Exceedance" = "VaR_violation") %>%
                     as_tibble()},
                 simplify=FALSE) %>%
  bind_rows(.id = "Country")
                   

                   
df_plot <- bind_rows(df_ES %>% filter(risk_measure %in% c("VaR", "ES")), df_MES) %>% 
  na.omit() %>%
  dplyr::rename("Risk Measure" = "risk_measure")


# Arrange facets and legend orders
df_plot$Country <- factor(df_plot$Country, 
                          levels=c("Joint Economic Region", "DEU", "FRA", "GBR"),
                          labels=c("Joint Economic Region", "Germany", "France", "United Kingdom"))

df_plot$`Risk Measure` <- factor(df_plot$`Risk Measure`, 
                          levels=c("VaR", "ES", "MES"))


# Set manual colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# ggplot2
p_ES <- ggplot(df_plot %>% arrange(Country, `VaR Exceedance`)) +
  ggplot2::geom_point(aes(x=Date, y=return, color=`VaR Exceedance`)) +
  ggplot2::scale_colour_manual(values = c("grey", "black")) +
  ggnewscale::new_scale_color() + # For two color scales!!!
  ggplot2::geom_line(aes(x=Date, y=value, color=`Risk Measure`)) +
  ggplot2::scale_colour_manual(values = c("VaR"=gg_color_hue(3)[3], "ES"=gg_color_hue(3)[2], "MES"=gg_color_hue(3)[1])) +
  ggplot2::facet_wrap(~Country, ncol=2) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position="bottom") +
  ggplot2::ylab("Negative GDP Growth")

p_ES

ggsave("applications/output/MES_GDPgrowth.pdf", width = 8, height = 6, units = "in")


# Show VaR exceedances
df_plot %>% 
  group_by(Country) %>% 
  summarize(sum(`VaR Exceedance`),
            mean(`VaR Exceedance`))












#### DCC GARCH model and Plots

# DCC GARCH specifications
GARCH11_norm_spec <- ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean=FALSE),
                                variance.model = list(garchOrder = c(1,1), model = "sGARCH"),
                                distribution.model = "norm")
DCC_spec <- dccspec(uspec = multispec(replicate(4, GARCH11_norm_spec)), dccOrder = c(1,1), distribution = "mvnorm")

# DCC fit
set.seed(2023)
DCC_fit_GDP <- dccfit(DCC_spec, 
                      data = data_fit %>% select(Date, Region, DEU, FRA, GBR) %>% as.xts,
                      fit.control = list(eval.se = FALSE, stationarity = TRUE, scale = FALSE))

# MES from DCC
H_GDP <- rcov(DCC_fit_GDP)
Sigma_GDP <- sapply(1:dim(H_GDP)[3], function(i) pracma::sqrtm(H_GDP[,,i])$B, simplify = "array")


# Preliminary save of DCC MES results...
DCC_risk_measure_df <- tibble(Date=data_fit$Date)

for (country in c("DEU", "FRA", "GBR")){
  vector_index_country <- c("DEU"=2, "FRA"=3, "GBR"=4)
  index_country <- vector_index_country[country]
  
  # Contribution to the MES without the ARMA part
  DCC_risk_measures_woARMA <- sapply(1:dim(H_GDP)[3], 
                                  function(i) DCC_to_CoVaR(Sigma_FC=Sigma_GDP[c(1,index_country),c(1,index_country),i], nu=10^4, data_IS=NULL, alpha=0.5, beta=beta),
                                  simplify = "array")
  #Manually add ES forecasts
  DCC_risk_measures_woARMA <- rbind(DCC_risk_measures_woARMA,
                                    as.numeric(DCC_risk_measures_woARMA[1,]) * dnorm(qnorm(beta))/(1-beta) / qnorm(beta))
  
  # Add to data frame and add ARMA part
  DCC_risk_measure_df <- DCC_risk_measure_df %>% 
    mutate(!!paste0("MES_",country) := as.numeric(DCC_risk_measures_woARMA[2,]) + as.numeric(fitted(DCC_fit_GDP)[,country]),
           VaR_Region = as.numeric(DCC_risk_measures_woARMA[1,]) + as.numeric(fitted(DCC_fit_GDP)[,"Region"]),
           ES_Region = as.numeric(DCC_risk_measures_woARMA[4,]) + as.numeric(fitted(DCC_fit_GDP)[,"Region"]))
}


DCC_risk_measure_df2 <- 
  DCC_risk_measure_df %>% 
  rename("MES_Germany"="MES_DEU",
         "MES_France"="MES_FRA",
         "MES_United Kingdom"="MES_GBR",
         "VaR_Joint Economic Region"="VaR_Region",
         "ES_Joint Economic Region"="ES_Region") %>%
  pivot_longer(cols=!Date,
               values_to="value",
               names_to=c("Risk Measure","Country"),
               names_pattern = "(.*)_(.*)") %>%
  filter(Date >= "1995-04-01")



# Join with MES reg data frame
df_plot_DCC <- bind_rows(df_plot %>% mutate(Model="MESreg"), 
                         DCC_risk_measure_df2 %>% mutate(Model="DCC"))


# Arrange facets and legend orders
df_plot_DCC$Country <- factor(df_plot_DCC$Country, 
                              levels=c("Joint Economic Region", "Germany", "France", "United Kingdom"),
                              labels=c("Joint Economic Region", "Germany", "France", "United Kingdom"))

df_plot_DCC$`Risk Measure` <- factor(df_plot_DCC$`Risk Measure`, 
                                     levels=c("VaR", "ES", "MES"))

df_plot_DCC$Model <- factor(df_plot_DCC$Model, 
                            levels=c("MESreg", "DCC"),
                            labels=c("MES regression", "DCC-GARCH"))


p_MES_DCC <- ggplot(df_plot_DCC %>% arrange(Country, `VaR Exceedance`)) +
  ggplot2::geom_point(aes(x=Date, y=return, color=`VaR Exceedance`)) +
  ggplot2::scale_colour_manual(values = c("grey", "black"), na.translate = F) +
  ggnewscale::new_scale_color() + # For two color scales!!!
  ggplot2::geom_line(aes(x=Date, y=value, color=`Risk Measure`, linetype=Model)) +
  ggplot2::scale_colour_manual(values = c("VaR"=gg_color_hue(3)[3], "ES"=gg_color_hue(3)[2], "MES"=gg_color_hue(3)[1])) +
  ggplot2::scale_linetype_manual(values = c("MES regression"="solid", "DCC-GARCH"="dashed")) +
  ggplot2::facet_wrap(~Country, ncol=2) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position="bottom") +
  ggplot2::ylab("Negative GDP Growth")

p_MES_DCC

ggsave("applications/output/MES_DCC_GDPgrowth.pdf", width = 10, height = 7, units = "in")







#### Model Diagnostics

for (country_select in c("DEU", "FRA", "GBR")){
    
  # MES regression generalized residuals
  df_ResidDiag_MESreg <- SRM_est_list[[country_select]]$data %>%
    tail(-1) %>%
    select(Date,x,y,VaR,risk_measure) %>%
    rename(MES=risk_measure) %>%
    mutate(GenResid_VaR = (x<=VaR) - beta,
           GenResid_MES_zeros = (x>VaR)*(MES-y),
           GenResid_MES = na_if(GenResid_MES_zeros, (x>VaR)),
           model="MESreg") %>%
    as_tibble
  
  
  # DCC generalized residuals
  df_ResidDiag_DCC <- DCC_risk_measure_df %>%
    full_join(data_fit %>% select(-FCI), by="Date") %>%
    tail(-1) %>%
    rename(all_of(c(x="Region", y=paste(country_select), VaR="VaR_Region", MES=paste0("MES_",country_select)))) %>%
    select(Date,x,y,VaR,MES) %>%
    mutate(GenResid_VaR = (x<=VaR) - beta,
           GenResid_MES_zeros = (x>VaR)*(MES-y),
           GenResid_MES = na_if(GenResid_MES_zeros, (x>VaR)),
           model="DCC") %>%
    as_tibble()
  
  
  
  ### Join MESreg, DCC and BRS data frames  
  df_ResidDiag <- bind_rows(df_ResidDiag_MESreg, df_ResidDiag_DCC) %>%
    arrange(Date, model)
  
  
  # Edit the data frame such that it fits in our plotting routine.
  df_joint_plot <- bind_rows(df_ResidDiag %>% 
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
  df_joint_plot$model <- factor(df_joint_plot$model,
                                levels=c("VaR_MESreg", "MES_MESreg", "VaR_DCC", "MES_DCC"),
                                labels=c("VaR Regression", "MES Regression", "VaR DCC-GARCH", "MES DCC-GARCH"))
  
  p_Date <- ggplot(df_joint_plot, aes(x=Date, y=GenResid)) +
    facet_wrap(~model, scales="free", nrow=1) +
    geom_point() + 
    geom_smooth(method="loess", fill="blue") + 
    geom_hline(yintercept=0) + 
    theme_bw() +
    # theme(aspect.ratio=1) +
    ylab("Generalized Residuals") +
    xlab("Date") 
  
  ggsave(paste0("applications/output/GDP_Diagnostics_Date_",country_select,".pdf"), 
         p_Date, width = 12, height = 3, units = "in")
  
  
  p_FC <- ggplot(df_joint_plot, aes(x=Prediction, y=GenResid)) +
    facet_wrap(~model, scales="free", nrow=1) +
    geom_point() + 
    geom_smooth(method="loess", fill="blue") + 
    geom_hline(yintercept=0) + 
    theme_bw() +
    # theme(aspect.ratio=1) +
    ylab("Generalized Residuals") +
    xlab("Model Predictions") 
  
  ggsave(paste0("applications/output/GDP_Diagnostics_Predictions_",country_select,".pdf"), 
         p_FC, width = 12, height = 3, units = "in")
  
}
  



