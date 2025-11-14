# 1.- Event detection and Event detection for FT's mean signals  ----------
## ------------------------------------------------------------------------#
## @date 2025-05-20
## @paper Assessing the spatio-temporal cumulative drought effects
##        on two Mediterranean forest types in Central Chile.
## @R version R version 4.3.2 (2023-10-31 ucrt)
## @OS system Windows 10 x64
## @author Jose A. Lastra
## @email jose.lastramunoz@wur.nl | jose.lastra@pucv.cl
## ------------------------------------------------------------------------#
# 2.- Libraries --------
library(tidyverse)

## Clean environment ----
rm(list = ls(all = T))

## source functions -----
source("RFUN/aux_functions.R")

## ---------------------------------------------------------------------#
# 3.- Event detection ----
## tables -----
tb_list <- fs::dir_ls(
  path = "02_Tabular_data/A_mean_time_series/",
  glob = "*mean_anom*.csv"
)

## read table -----
variable <- "EVI" ## change accordingly
anom_accum_tb <- read_csv(str_subset(tb_list, variable))

## Detect events ----
### setup ----
max_events <- 10
max_gap_days <- 60
rfd_thr <- 90
min_tud <- 1 # to retain all the events this was set to 1 day
output <- "all"
rfd_vals <- rfd_polarity(
  rfd_vals = anom_accum_tb$dec_rfd,
  polarity = sign(anom_accum_tb$dec_anom)
)

### events deciduous ----
events <- fixed_event_detection(
  accum_anom = anom_accum_tb$dec_cusum,
  rfd_vals = rfd_vals,
  dates = anom_accum_tb$dates,
  max_events = max_events,
  max_gap_days = max_gap_days,
  output = output,
  min_tud = min_tud
)

start_dates <- events[names(events) == "start_date"] %>% ym()
end_dates <- events[names(events) == "end_date"] %>% ym()
auc_values <- events[names(events) == "auc"]
polarity <- events[names(events) == "polarity"]
tud_days <- events[names(events) == "tud_days"]
tud_months <- events[names(events) == "tud_months"]


### events sclerophyllous ----
rfd_vals <- rfd_polarity(
  rfd_vals = anom_accum_tb$scl_rfd,
  polarity = sign(anom_accum_tb$scl_anom)
)

events <- fixed_event_detection(
  accum_anom = anom_accum_tb$scl_cusum,
  rfd_vals = rfd_vals,
  dates = anom_accum_tb$dates,
  max_events = max_events,
  max_gap_days = max_gap_days,
  output = output,
  min_tud = min_tud
)

start_dates <- events[names(events) == "start_date"] %>% ym()
end_dates <- events[names(events) == "end_date"] %>% ym()
auc_values <- events[names(events) == "auc"]
polarity <- events[names(events) == "polarity"]
tud_days <- events[names(events) == "tud_days"]
tud_months <- events[names(events) == "tud_months"]
