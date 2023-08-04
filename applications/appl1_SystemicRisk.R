library(tidyverse)
library(lubridate)
library(tsibble)
library(SystemicRisk)


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

data_VIX <- readRDS(file = "applications/data/data_VIX.rds") 

data_Assets_long <- readRDS(file = "applications/data/data_Assets.rds") 

data_Assets_wide <- data_Assets_long %>% 
  pivot_wider(names_from = Asset, values_from = NegReturn)
  
data_joint <- full_join(data_Spreads, data_Assets_wide, by="Date") %>%
  full_join(data_VIX, by="Date") %>%
  full_join(data_Bubbles, by="Date") %>%
  dplyr::select(Date, ChangeSpread, TEDSpread, VIX, st_boom_bsadf, st_bust_bsadf, 
         BAC, C, JPM, SPF) %>%  
  dplyr::filter(Date >= date("1990-01-01") & Date <= date(" 2015-12-31")) %>%
  arrange(Date) %>%
  fill(c(st_boom_bsadf, st_bust_bsadf), .direction = "down") %>%
  na.omit()

x_set <- c( "BAC", "C", "JPM")
x_name <- c("BAC", "C", "JPM")

MES_obj_list <- list()
MES_dfsummary_list <- list()
BRSreg_obj_list <- list()
BRSreg_dfsummary_list <- list()

for (i in 1:length(x_set)){
  # Transform to a tsibble, where x is SPF, and y the individual banks
  df <- data_joint %>%
    summarize(Date=Date,
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
  # Note that we need "lead" for the indices as the joint linear model automatically shifts the 
  # variables in the package, but BRS(2020) use contemporaneous indices. 
  
  set.seed(2023)
  MESreg_obj <- SRM(data=df,
                    model="joint_linear",
                    risk_measure="MES",
                    beta=0.95,
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
           VIX=lag(VIX))
  # Note: Now we lag the covariates as we the "lm" function below does not use lags automatically!
           
  
  lmBRSreg <- lm(y_MESrolling ~ ChangeSpread + TEDSpread + VIX  + ST_boom + ST_bust,
                 data_BRSreg)
  BRSreg_obj_list[[x_name[i]]] <- lmBRSreg
  
  sum_BRSreg <- summary(lmBRSreg)
  BRSreg_dfsummary_list[[x_name[i]]] <- as_tibble(sum_BRSreg$coefficients, rownames="Covariate") %>%
    rename(Estimate=2, Error=3, tval=4, pval=5)
  
}



# Look at the parameter estimates
MES_dfsummary_list

BRSreg_dfsummary_list


# Print results to a LaTeX table
full_join(bind_rows(MES_dfsummary_list, .id = "Asset"),
          bind_rows(BRSreg_dfsummary_list, .id = "Asset") %>% 
            mutate(Covariate = ifelse(Covariate=="(Intercept)", "Intercept", Covariate)),
          by=c("Asset", "Covariate")) %>%
  mutate(empty1=NA, empty2=NA, empty3=NA, empty4=NA) %>%
  dplyr::select(Asset, empty1, Covariate, 
                empty2, Estimate_VaR, Error_VaR,  pval_VaR,
                empty3, Estimate_MES, Error_MES, pval_MES,
                empty4, Estimate, Error, pval) %>%
  xtable::xtable(digits=c(0, 0, 0, 0, 0, c(3,3,3), 0, c(3,3,3), 0, c(3,3,3))) %>%
  print(file="applications/output/MES_AB_regressions.txt", include.rownames=FALSE, booktabs=TRUE)






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




### Look at the model residuals to show that the BRS regression residuals have some time-variation.
# ggplot(df_MESBRSreg_plot %>%
#          dplyr::filter(`VaR Exceedance`==T) %>%
#          mutate(residuals=y-value)) +
#   geom_point(aes(x=Date, y=residuals, col=Method)) +
#   geom_smooth(aes(x=Date, y=residuals)) +
#   facet_grid(~Method)







