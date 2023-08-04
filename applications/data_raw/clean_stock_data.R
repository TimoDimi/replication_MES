# install.packages('quantmod')
library(quantmod)
library(dplyr)
library(tibble)
library(lubridate)
library(reshape2)

Symbol_list <- c("AAPL", "MSFT", "IBM",  "JPM", "AMZN", "BAC", "C", "GS", "UNH", "XOM", "DIS", "NFLX", "BRK-A", "JNJ", "LLY", "BA", "^VIX","^SPX", "^SP500-40")
Symbol_names <- c("AAPL", "MSFT", "IBM", "JPM", "AMZN", "BAC", "C", "GS", "UNH", "XOM", "DIS", "NFLX", "BRK", "JNJ", "LLY", "BA", "VIX", "SP500", "SPF")

# Download assets for the symbol list
data_Symbols <- lapply(Symbol_list, function(x) {
  getSymbols(x,
             from = "1990/01/01",
             to = "2023/05/31",
             periodicity = "daily",
             auto.assign = FALSE) %>%
    data.frame(Date=index(.), check.names=FALSE) %>%
    tibble::as_tibble() %>%
    rename_all(~stringr::str_replace_all(., paste0(x,"."), ""))
})
names(data_Symbols) <- Symbol_names


# Correct a mistake in the data
data_Symbols$C <- data_Symbols$C %>% 
  rename(Close=ose)


# Compute neg returns for all assets but VIX (in long format)
data_Assets <- bind_rows(data_Symbols, .id = "Asset") %>%
  dplyr::filter(Asset != "VIX") %>%
  dplyr::group_by(Asset) %>%
  dplyr::mutate(Date=lubridate::as_date(Date),
                NegReturn= -100*(log(Close) - log(lag(Close)))) %>%
  dplyr::select(Date, Asset, NegReturn)

data_Assets %>% group_by(Asset) %>% summarize(sum_NA = sum(is.na(NegReturn)))


# Look at NA values in the assets
data_Assets_wide <- data_Assets %>% reshape2::dcast(Date~Asset) %>% as_tibble()
data_Assets_wide[rowSums(is.na(data_Assets_wide)) > 0, ]

# Optional: Delete NAs uniformly through all assets!
# data_Assets <- data_Assets %>%
#   reshape2::dcast(Date~Asset) %>%
#   na.omit() %>%
#   reshape2::melt(id.vars = c("Date"), variable.name="Asset", value.name="NegReturn") %>%
#   tibble::as_tibble()

# Look at some specifics
data_Assets %>% group_by(Asset) %>% summarize(n=n())

saveRDS(data_Assets, file = "applications/data/data_Assets.rds")



# Special case: The VIX index
saveRDS(data_Symbols$VIX %>%
          dplyr::select(Date, Close) %>%
          rename(VIX=Close),
        file = "applications/data/data_VIX.rds")


