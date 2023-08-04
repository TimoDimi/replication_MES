library(tidyverse)
library(lubridate)
library(readxl)
library(slider)



################################################################################
### (1) Load buddle indicator data
################################################################################

# Load data
data_bubbles_raw <- read_excel("applications/data_raw/BubbleEpisodes.xlsx", sheet = 2, skip=0)

data_bubbles <- data_bubbles_raw %>%
  mutate(Date=date(t)) %>% 
  dplyr::filter(country=="us") %>% 
  select(-t)

# Save for now as individual file
saveRDS(data_bubbles, "applications/data/data_bubbles.rds")


