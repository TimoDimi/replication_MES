library(dplyr)
library(tidyverse)
library(tsibble)
library(lubridate)
library(SystemicRisk)

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


# Compare with a two-step ES regression
set.seed(2023)
SRM_ES <- SRM(data=tsibble_fit %>% mutate(y=x),
              model = "joint_linear",
              risk_measure="MES",
              beta=0.9,
              optim_replications=c(3,10)) 

summary(SRM_ES)


# Print the MES coefficients
SRM_df %>%
  dplyr::filter(measure=="MES") %>%
  arrange(Variable) %>%
  mutate(CI_lower = Estimate - qnorm(0.95) * `Std. Error`,
         CI_upper = Estimate + qnorm(0.95) * `Std. Error`) 



# Joint plot of in-sample fits
df_ES <- SRM_ES$data %>%
  dplyr::mutate(VaR_violation=(x > VaR),
                Country="Joint Economic Region") %>%
  dplyr::select(Date, Country, VaR_violation, x, VaR, risk_measure) %>%
  dplyr::rename(ES=risk_measure,
                return=x) %>%
  pivot_longer(cols=c(VaR,ES), names_to="risk_measure") %>%
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
                   

                   
df_plot <- bind_rows(df_ES, df_MES) %>% 
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


