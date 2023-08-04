library(tidyverse)
library(quantmod)

# ===================================================
# Import FED Data
# ===================================================

data_ChangeSpread <- read_csv("./applications/data_raw/BAA10Y.csv",
                              col_types="Dd")

data_Spread <- read_csv("./applications/data_raw/T10Y3M.csv",
                        col_types="Dd")

data_TEDSpread <- read_csv("./applications/data_raw/TEDRATE.csv",
                           col_types="Dd")

# merge all data frames through a list
data_Spreads <- list(data_ChangeSpread, data_Spread, data_TEDSpread) %>%
  reduce(full_join, by='DATE') %>%
  rename(Date=DATE,
         ChangeSpread=BAA10Y,
         Spread=T10Y3M,
         TEDSpread=TEDRATE)

saveRDS(data_Spreads, file = "applications/data/data_Spreads.rds")


