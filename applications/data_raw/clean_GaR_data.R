library(tidyverse)
library(lubridate)
library(readxl)
library(slider)

################################################################################
### (1) Download and clean absolute GDP data for G7 countries
################################################################################

G7countries_list <-  c("USA","CAN","DEU","FRA","GBR","JPN","ITA")

# Data downloaded from the OECD data archive; by manual selections here: 
# https://data.oecd.org/gdp/gross-domestic-product-gdp.htm#indicator-chart
data_GDPtot_raw <- read_csv(file = "applications/data_raw/data_GDPtotal_OECD.csv")

data_GDPtot <- data_GDPtot_raw %>%
  filter(LOCATION %in% G7countries_list) %>%
  mutate(Date = ymd(sprintf("%d-01-01", TIME))) %>%
  rename(Country=LOCATION, GDP=Value) %>%
  dplyr::select(Date, Country, GDP) %>%
  arrange(Date, Country) %>%
  group_by(Date) %>%
  mutate(GDPrel = GDP/sum(GDP))

saveRDS(data_GDPtot, "applications/data/data_GaR_total.rds")




################################################################################
### (2) Download and clean GDP growth data
################################################################################

# Data downloaded from the OECD data archive (by manual selections here:
# https://data.oecd.org/gdp/gross-domestic-product-gdp.htm#indicator-chart
data_GDP_raw <- read_csv(file = "applications/data_raw/data_GDPgrowth_OECD_PreviousYear.csv")

data_GDP <- data_GDP_raw %>%
  mutate(Date = yq(TIME)) %>%
  rename(Country=LOCATION, GDP_Growth=Value) %>%
  dplyr::select(Date, Country, GDP_Growth) %>%
  arrange(Date, Country)

# Have a look at the data
ggplot(data_GDP) +
  geom_line(aes(x=Date, y=GDP_Growth, col=Country)) +
  facet_wrap(~Country)

saveRDS(data_GDP, "applications/data/data_GaR.rds")


################################################################################
### (3) Download and clean global FCI data from Arrigoni et al (2022, IMF ER)
################################################################################

# Download from https://sites.google.com/view/simonearrigoni/data?authuser=0
data_gFCI_raw <- read_excel("applications/data_raw/FCIs_ABV_IMFER_2022.xlsx", sheet = 2)

# Select variables, and convert monthly to quarterly data through averaging
data_gFCI <- data_gFCI_raw %>%
  mutate(Date=date(date)) %>% 
  dplyr::select(Date, iso_alpha3, 
                ew_fci, pca_fci, tvp_fci, 
                global_ew_fci, global_tvp_fci, global_pca_fci,
                ae_ew_fci, ae_tvp_fci, ae_pca_fci, 
                ea_ew_fci, ea_tvp_fci, ea_pca_fci, 
                share) %>%
  rename(Country=iso_alpha3)

# Monthly to Quarterly frequency through averaging
data_gFCI <- data_gFCI  %>%
  mutate(Quarter=yq(quarter(Date, with_year=T))) %>%
  group_by(Country, Quarter) %>%
  summarize(across(!Date, mean)) %>% 
  rename(Date=Quarter)              

# Save for now as individual file, maybe generalize later and merge to data_GaR
saveRDS(data_gFCI, "applications/data/data_gFCI.rds")





