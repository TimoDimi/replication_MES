library(tidyverse)
library(tsibble)
library(reshape2)
library(rmgarch)
library(mvtnorm)
library(cubature)
library(lubridate)
library(doParallel)
library(SystemicRisk)
library(patchwork)



# Set manual colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


################################################################################
###    Manually set options
################################################################################

start_date <- "2006-12-25"
asset_names <- c("BA", "JPM", "LLY", "MSFT", "XOM") 

# Set manual colors
cols_manual <- c("BA"=gg_color_hue(6)[2],
                 "JPM"=gg_color_hue(6)[3],
                 "LLY"=gg_color_hue(6)[4],
                 "MSFT"=gg_color_hue(6)[5],
                 "XOM"=gg_color_hue(6)[6],
                 "ERC Portfolio"="red",
                 "EW Portfolio"="black",
                 "GMVP-DCC Portfolio"="blue")

cols_stocks_manual <- c("BA"=gg_color_hue(6)[2],
                        "JPM"=gg_color_hue(6)[3],
                        "LLY"=gg_color_hue(6)[4],
                        "MSFT"=gg_color_hue(6)[5],
                        "XOM"=gg_color_hue(6)[6])



# Load data and transform
data_mod <- readRDS(file = "applications/data/data_mod_appl_ERC_20230802.rds")
res_appl_GMVP <- readRDS(file = "applications/data/res_appl_GMVP_20230802.rds")

res_appl <- readRDS(file = "applications/data/res_appl_ERC_20230802.rds") %>%
  group_by(Asset) %>%
  mutate(Date = lead(Date)) %>% # Lead Date by one week to match with the realizations and GMVP weights
  ungroup() %>%
  group_by(Date) %>%
  mutate(MES_rel = MES/sum(MES),
         MES_EW_rel = MES_EW/sum(MES_EW),
         w_stacked=cumsum(w))



# Stretch the res_appl data set to all dates with constant weights in the case that we do not refit every corresponding day or week!
# merge by (date, asset) with the log-losses df; multiply together through mutate
# Careful: The "date" column corresponds to the beginning of the last week used for fitting, but as we forecast MES, I think we have to shift by one week!
df_merged <- full_join(data_mod %>%
                         dplyr::filter(Date >= start_date) %>%
                         dplyr::select(Date, all_of(asset_names)) %>%
                         pivot_longer(cols=asset_names, names_to="Asset", values_to="LogLosses"),
                       res_appl,
                       by=c("Date","Asset")) %>%
  full_join(res_appl_GMVP, by=c("Date","Asset")) %>% 
  arrange(Date, Asset) %>%
  group_by(Asset) %>%
  fill(w, MES) %>% # This is necessary for days where we do not have a convergence
  ungroup() %>%
  na.omit()


# Compute tibble with portfolio returns/and relative prices
df_P <- df_merged %>%
  group_by(Date) %>%
  summarize(Portf_ERC=sum(LogLosses*w),
            Portf_GMVP_DCC=sum(LogLosses*w_GMVP_DCC_normalized),
            # Portf_GMVP_DCC=sum(LogLosses*w_GMVP_DCC),
            Portf_EW=sum(LogLosses*1/m)) %>%
  pivot_longer(cols=c("Portf_ERC", "Portf_EW", "Portf_GMVP_DCC"), names_to="Asset", values_to="LogLosses")


# Append Portfolio and Stock data frames
df_performance <- bind_rows(df_P,
                            df_merged %>% dplyr::select(Date, Asset, LogLosses)) %>%
  mutate(type_Asset=ifelse((Asset %in% c("Portf_EW","Portf_ERC", "Portf_GMVP_DCC")), "Portfolio", "Stock")) %>%
  arrange(Date, Asset) %>%
  group_by(Asset) %>%
  mutate(Price_rel=exp(cumsum(-LogLosses/100)))


### Table 3: Portfolio performance
df_performance %>%
  dplyr::filter(Date > as_date("2007-01-01"), Date < as_date("2024-01-01")) %>%
  dplyr::filter(Asset %in% c("Portf_EW","Portf_ERC", "Portf_GMVP_DCC")) %>%
  group_by(Asset) %>%
  summarize(mean=mean(-LogLosses),
            sd=sd(-LogLosses),
            VaR=quantile(-LogLosses, 0.025),
            ES=mean(-LogLosses[-LogLosses < VaR]),
            SharpR=mean/sd,
            RORAC=mean/-ES) %>%
  arrange(desc(RORAC))





################################################################################
###   Plot Portfolio Performance 
################################################################################

df_performance_full <- df_performance %>%
  bind_rows(tibble(Date=date("2007-01-01"), Asset=unique(df_performance$Asset), type_Asset=c(rep("Stock",length(asset_names)),rep("Portfolio",3)), LogLosses=0, Price_rel=NA)) %>%
  arrange(Date, Asset) %>%
  group_by(Asset) %>%
  mutate(Price_rel=exp(cumsum(-LogLosses/100)),
         Asset=recode(Asset, "Portf_ERC"="ERC Portfolio", "Portf_EW"="EW Portfolio", "Portf_GMVP_DCC"="GMVP-DCC Portfolio"), ##, "Portf_GMVP_uc"="GMVP-uc Portfolio"),
         Asset=factor(Asset))


df_performance_full$Asset <- factor(df_performance_full$Asset, 
                                    levels = c(asset_names, "EW Portfolio", "ERC Portfolio", "GMVP-DCC Portfolio")) ##, "GMVP-uc Portfolio"))



# Plot performance
ggplot(df_performance_full %>% arrange(Asset)) +
  geom_line(aes(x=Date, y=Price_rel, col=Asset, linetype=type_Asset, linewidth=type_Asset)) +
  scale_color_manual(values = cols_manual) +
  ylab("Normalized Value") +
  scale_x_date(date_labels = "%Y",
               date_minor_breaks = "1 year",
               breaks = function(x) seq.Date(from = date("2007-01-01"),
                                             to = date("2023-01-01"),
                                             by = "3 years")) +
  scale_fill_manual(values = cols_manual) +
  scale_linetype_discrete(guide = "none") +
  scale_discrete_manual("linewidth", values = c(1, 0.5), guide = "none") +
  theme_bw() +
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 1)) + 
  coord_cartesian(ylim=c(0.5,4))



################################################################################
###  Plot ABSOLUTE MES forecasts of the two portfolio
################################################################################
res_MES_abs <- res_appl %>%
  rename(EW = MES_EW, ERC=MES) %>%
  pivot_longer(cols=c("EW", "ERC"), names_to="Portfolio", values_to="MES") %>%
  arrange(Asset) %>%
  group_by(Date, Portfolio) %>%
  mutate(MES_stacked = cumsum(MES)) %>%
  arrange(Date, Portfolio, Asset) %>%
  dplyr::select(Date, Portfolio, Asset, MES, MES_stacked)

res_MES_rel <- res_appl %>%
  rename(EW = MES_EW_rel, ERC=MES_rel) %>%
  dplyr::select(Date, Asset, EW, ERC) %>%
  pivot_longer(cols=c("EW", "ERC"), names_to="Portfolio", values_to="MES") %>%
  arrange(Asset) %>%
  group_by(Date, Portfolio) %>%
  mutate(MES_stacked = cumsum(MES)) %>%
  arrange(Date, Portfolio, Asset) 


p1 <- ggplot(res_MES_abs %>% filter(Portfolio=="EW")) +
  geom_area(aes(x=Date, y=MES, fill=Asset)) +
  scale_fill_manual(values = cols_manual) +
  scale_y_continuous(name="Absolute Risk Contributions") +
  scale_x_date(date_labels = "%Y",
               date_minor_breaks = "1 year",
               breaks = function(x) seq.Date(from = date("2008-01-01"),
                                             to = date("2023-01-01"),
                                             by = "3 years")) +
  theme_bw() + 
  theme(legend.position = "bottom")

p2 <- ggplot(res_MES_rel %>% filter(Portfolio=="EW")) +
  geom_area(aes(x=Date, y=MES, fill=Asset)) +
  scale_fill_manual(values = cols_manual) +
  geom_hline(yintercept=c(seq(0,1,0.2)), linetype="dotted", alpha=0.8) +
  scale_y_continuous(name="Relative Risk Contributions",
                     breaks=seq(0,1,0.2)) +
  scale_x_date(date_labels = "%Y",
               date_minor_breaks = "1 year",
               breaks = function(x) seq.Date(from = date("2008-01-01"),
                                             to = date("2023-01-01"),
                                             by = "3 years")) +
  theme_bw() + 
  theme(legend.position = "bottom")

((p1 / p2) & theme(legend.position = "bottom")) + plot_layout(guides = "collect")

ggsave("applications/output/MES_EW.pdf", width = 8, height = 7, units = "in")


################################################################################
###  Area plot of ERC portfolio weights
################################################################################

ggplot(res_appl) +
  geom_area(aes(x=Date, y=w, fill=Asset)) +
  geom_hline(yintercept=c(seq(0,1,0.2)), linetype="dotted", alpha=0.8) +
  scale_fill_manual(values = cols_manual) +
  scale_y_continuous(name="ERC Portfolio Weights",
                     breaks=seq(0,1,0.2)) +
  scale_x_date(date_labels = "%Y",
               date_minor_breaks = "1 year",
               breaks = function(x) seq.Date(from = date("2008-01-01"),
                                             to = date("2023-01-01"),
                                             by = "3 years")) +
  theme_bw() + 
  theme(legend.position = "bottom")

# ggsave("applications/output/ERC_Portfolio_Weights.pdf", width = 8, height = 4, units = "in")




################################################################################
###  Area plot of GMVP portfolio weights
################################################################################

ggplot(res_appl_GMVP %>% ungroup) +
  geom_area(aes(x=Date, y=w_GMVP_DCC_normalized, fill=Asset)) +
  geom_hline(yintercept=c(seq(0,1,0.2)), linetype="dotted", alpha=0.8) +
  scale_fill_manual(values = cols_manual) +
  scale_y_continuous(name="GMVP Portfolio Weights",
                     breaks=seq(0,1,0.2)) +
  scale_x_date(date_labels = "%Y",
               date_minor_breaks = "1 year",
               breaks = function(x) seq.Date(from = date("2008-01-01"),
                                             to = date("2023-01-01"),
                                             by = "3 years")) +
  theme_bw() + 
  theme(legend.position = "bottom")

# ggsave("applications/output/GMVP_Portfolio_Weights.pdf", width = 8, height = 4, units = "in")



################################################################################
###  Joint area plot of the ERC and GMVP portfolio weights
################################################################################

df_weights <- full_join(res_appl %>% dplyr::select(Date, Asset, w),
                        res_appl_GMVP %>% ungroup() %>% dplyr::select(Date, Asset, w_GMVP_DCC_normalized),
                        by=c("Date", "Asset")) %>%
  rename(ERC=w, GMVP=w_GMVP_DCC_normalized) %>%
  pivot_longer(cols=c(ERC, GMVP), names_to="Portfolio", values_to="Weights")

ggplot(df_weights) +
  geom_area(aes(x=Date, y=Weights, fill=Asset)) +
  geom_hline(yintercept=c(seq(0,1,0.2)), linetype="dotted", alpha=0.8) +
  scale_fill_manual(values = cols_manual) +
  scale_y_continuous(name="Portfolio Weights",
                     breaks=seq(0,1,0.2)) +
  scale_x_date(date_labels = "%Y",
               date_minor_breaks = "1 year",
               breaks = function(x) seq.Date(from = date("2008-01-01"),
                                             to = date("2023-01-01"),
                                             by = "3 years")) +
  facet_wrap(~Portfolio, ncol=1) + 
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("applications/output/Portfolio_Weights_ERC_GMVP.pdf", width = 8, height = 7, units = "in")



