# 1.- Event detection for raster data ----------
## ------------------------------------------------------------------------#
## @date 2025-03-26
## @paper Assessing the spatio-temporal cumulative drought effects
##        on two Mediterranean forest types in Central Chile.
## @R version R version 4.3.2 (2023-10-31 ucrt)
## @OS system Windows 10 x64
## @author Jose A. Lastra
## @email jose.lastramunoz@wur.nl | jose.lastra@pucv.cl
## ------------------------------------------------------------------------#
# 2.- Libraries --------
pacman::p_load(tidyverse, terra)

## Clean environment ----
rm(list = ls(all = T))

source("RFUN/aux_functions.R")
## ------------------------------------------------------------------------#
# 3.- Event detection  ----
## Spat raster data ----
### Cusums -----
variable <- "EVI" # change accordingly EVI, MSI, ET

path_acum <- switch(variable,
  "EVI" = "01_Raster_data/A_Cusums/EVI_cusum_2000_2024_MCD43A4.061_forest.tif",
  "MSI" = "01_Raster_data/A_Cusums/MSI_cusum_2000_2024_MCD43A4.061_forest.tif",
  "ET" = "01_Raster_data/A_Cusums/ET_cusum_2000_2024_MODIS16A2GF.061_forest.tif"
)

spat_accum <- rast(path_acum)

### RFD polarity -----
path_rfd_pol <- switch(variable,
  "EVI" = "01_Raster_data/B_RFD/EVI_RFD_polarity_2000_2024_MCD43A4.061_forest.tif",
  "MSI" = "01_Raster_data/B_RFD/MSI_RFD_polarity_2000_2024_MCD43A4.061_forest.tif",
  "ET" = "01_Raster_data/B_RFD/ET_RFD_polarity_2000_2024_MODIS16A2GF.061_forest.tif"
)

spat_rfd_pol <- rast(path_rfd_pol)


### combine data for detection ------
# 1:8764 ; 8765:17528 # for daily MODIS BRDF data
# 1:1104; 1105:2208 # For MODIS 8-day ET data

spat_detection <- c(spat_accum_crop, spat_rfd_crop)

## setup ----
dates <- names(spat_accum) %>%
  str_extract("\\d{4}-\\d{2}-\\d{2}") %>%
  ymd()

## Clean memory
rm(spat_accum, spat_rfd_pol)
gc()


max_events <- 20
max_gap_days <- 60
min_tud <- 30 # in days
output <- "negative" # change this accordingly negative, positive, all
outpath <- switch(output,
  "negative" = "negative",
  "positive" = "positive",
  "all" = "all_events"
)
rfd_thr <- 90 # extreme values characterization for events

outname <- paste0(
  "01_Raster_data/C_Events/", outpath, "/", variable, "_",
  output, "_events_", ifelse(variable == "ET", "MODIS16A2GF_v2.tif", "MCD43A4_v2.tif")
)

### event detection and TUD metrics ----
#### cluster setup -----
library(parallel)
ncores <- detectCores() - 1 # leave one core free
cl <- makeCluster(ncores)
clusterExport(cl, c(
  "fixed_event_detection", "dates",
  "max_events", "max_gap_days",
  "rfd_thr", "min_tud", "output"
))

# Cargar paquetes en los workers
clusterEvalQ(cl, {
  library(terra)
  library(dplyr)
  library(lubridate)
  library(magrittr)
  library(purrr)
})

system.time(spat_events <- app(spat_detection, fun = function(x) {
  # BRDF 1:8764 ; 8765:17528 | ET  1:1104; 1105:2208

  # Call vector function
  fixed_event_detection(
    accum_anom = as.numeric(x[1:1104]), # change subset values accordingly
    rfd_vals = as.numeric(x[1105:2208]), # change subset values accordingly
    output = output,
    rfd_thr = rfd_thr,
    dates = dates,
    max_events = max_events,
    max_gap_days = max_gap_days,
    min_tud = min_tud
  )
}, cores = cl, filename = outname, overwrite = T, wopt = list(datatype = "INT4S")))

stopCluster(cl)
