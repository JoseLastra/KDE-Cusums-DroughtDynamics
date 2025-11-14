# 1.- LM 2.6.	Environmental and topographic influence  ----------
## ------------------------------------------------------------------------#
## @date 2025-09-23
## @paper Assessing the spatio-temporal cumulative drought effects
##        on two Mediterranean forest types in Central Chile.
## @R version R version 4.3.2 (2023-10-31 ucrt)
## @OS system Windows 10 x64
## @author Jose A. Lastra
## @email jose.lastramunoz@wur.nl | jose.lastra@pucv.cl
## ------------------------------------------------------------------------#
# 2.- Libraries --------
pacman::p_load(tidyverse)
## Clean environment ----
rm(list = ls(all = T))
gc()

## ------------------------------------------------------------------------#
# 3.- Modelling setup -----
## load data ------
data_climate_terrain <- read_csv("02_Tabular_data/B_model_data/MODIS_forest_strata_sample_1000N.csv")

## add auxiliary fields and scale data -------
data_model <- data_climate_terrain %>%
  mutate(
    year = year(dates),
    month = month(dates),
    tmmx = tmmx * 0.1,
    tmmn = tmmn * 0.1,
    sm = sm * 0.1,
    aspect_sin = sin(aspect),
    aspect_cos = cos(aspect)
  )


### monthly climatology 2000-2010------
vars_standard <- c("pr", "tmmx", "tmmn", "sm")

climatology <- data_model %>%
  filter(between(year, 2000, 2010)) %>%
  group_by(fid, x, y, DN, month) %>%
  summarise(
    across(
      all_of(vars_standard), list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  )


## ------------------------------------------------------------------------#
# 4.- Modelling all Growing seasons (GS) -----
## Base setup ----
variable <- "EVI" # change accordingly EVI, MSI, ET
response <- paste0(tolower(variable), "_acum") # response variable
# DN: forest type
formula_lm <- as.formula(
  paste(response, "~ DN + pr_zacum + tmmx_zacum + sm_zacum +
         elev + slope + aspect_sin + aspect_cos + TPI + cwd +
         clay_0_5cm + total_awc")
)

### Z-score anomalies and Cusum calculation -----
gs_anom <- data_model %>%
  # join climatology
  left_join(climatology, by = c("fid", "month")) %>%
  # z-score for climate data
  mutate(
    pr_z = ifelse(pr_sd == 0, NA, (pr - pr_mean) / pr_sd),
    tmmx_z = ifelse(tmmx_sd == 0, NA, (tmmx - tmmx_mean) / tmmx_sd),
    tmmn_z = ifelse(tmmn_sd == 0, NA, (tmmn - tmmn_mean) / tmmn_sd),
    sm_z = ifelse(sm_sd == 0, NA, (sm - sm_mean) / sm_sd),
  ) %>%
  group_by(fid, x.x, y.x, DN.x) %>%
  # cusum climate anomalies
  mutate(
    pr_zacum = collapse::fcumsum(pr_z),
    tmmx_zacum = collapse::fcumsum(tmmx_z),
    sm_zacum = collapse::fcumsum(sm_z)
  ) %>%
  ungroup()


## Initiate GS modelling -----
yy_init <- 2000:2023 # define init years

model_list <- lapply(yy_init, function(yy) {
  # dates to filter
  gs_ini <- paste0(yy, "-07-01")
  gs_fin <- paste0(yy + 1, "-06-30")

  # prepare data set sum all temporal variables, and getting static variables values
  data_model_lpp <- gs_anom %>%
    rename(x = x.x, y = y.x, DN = DN.x) %>% # rename for convinience
    filter(between(dates, ymd(gs_ini), ymd(gs_fin))) %>% # filtering
    group_by(fid, x, y, DN) %>%
    # create data for model
    summarise(
      pr_zacum = sum(pr_zacum, na.rm = T),
      sm_zacum = sum(sm_zacum, na.rm = T),
      tmmx_zacum = sum(tmmx_zacum, na.rm = T),
      cwd = sum(cwd, na.rm = T),
      evi_acum = sum(evi_acum, na.rm = T),
      msi_acum = sum(msi_acum, na.rm = T),
      et_acum = sum(et_acum, na.rm = T),
      bulkd_0_5cm = first(bulkd_0_5cm),
      clay_0_5cm = first(clay_0_5cm),
      slope = first(slope),
      elev = first(elev),
      aspect_sin = sin(first(aspect)),
      aspect_cos = cos(first(aspect)),
      total_awc = first(total_awc),
      TPI = first(TPI),
      .groups = "drop"
    ) %>%
    drop_na()

  # add levels to forest type column
  data_model_lpp$DN <- factor(data_model_lpp$DN, levels = c(1, 2), labels = c("deciduous", "evergreen"))

  # build model
  m <- lm(
    formula = formula_lm,
    data = data_model_lpp
  )

  rm(data_model_lpp) # clean memory
  m
})


### Tidy up model output -----
names(model_list) <- yy_init

#### Broom data to tibble ------
summary_list <- lapply(model_list, broom::tidy)

#### R squared calculations for the global model using f-stat ----

rsquared <- lapply(model_list, function(x) {
  rsq <- 1 - sum(resid(x)^2) / sum((x$model[[1]] - mean(x$model[[1]]))^2)

  tibble(
    R2 = rsq %>% round(3),
    F_value = summary.lm(x)$fstatistic[1],
    df1 = summary.lm(x)$fstatistic[2],
    df2 = summary.lm(x)$fstatistic[3],
    p_global = pf(F_value,
      df1,
      df2,
      lower.tail = FALSE
    )
  )
}) %>%
  # bind all models
  bind_rows() %>%
  # add data for easy reading.
  mutate(
    year = as.character(yy_init),
    GS = paste0(year, "-", substr(as.numeric(year) + 1, 3, 4)),
    color_p_global = if_else(p_global < 0.05, "p-val < 0.05", "p-val > 0.05"),
    model = variable
  ) %>%
  select(GS, R2, p_global, year, color_p_global, model)

rsq_outname <- paste0("02_Tabular_data/B_model_data/results/RSquared_", variable, "_allGS_models.csv")

write_csv(rsquared, rsq_outname)

#### Combine results from all GS  with Rsquared -----
results <- bind_rows(summary_list, .id = "year")

results_full <- results %>%
  left_join(rsquared, by = "year") %>%
  mutate(
    signo = ifelse(estimate >= 0, "+", "–"),
    model = variable
  )


##### Filter based on significance ocurrence (>= 2) ------
## calculate ocurrence
vars_signif <- results_full %>%
  group_by(term) %>%
  summarise(signif_once = sum(p.value < 0.05, na.rm = T)) %>%
  filter(signif_once >= 2) %>%
  pull(term) %>%
  unique()

## filter
results_filter <- results_full %>%
  filter(term %in% vars_signif) %>%
  mutate(
    color_sig = if_else(p.value < 0.05, "p-val < 0.05", "p-val > 0.05"),
    color_p_global = if_else(p_global < 0.05, "p-val < 0.05", "p-val > 0.05")
  ) %>%
  ## Rename terms for easier interpretation
  mutate(estimates_names = case_when(
    term == "DNevergreen" ~ "Forest type (sclerophyllous)",
    term == "tmmx_zacum" ~ "Max. Temperature (z-cumsum)",
    term == "pr_zacum" ~ "Precipitation (z-cumsum)",
    term == "sm_zacum" ~ "Soil moisture (z-cumsum)",
    term == "elev" ~ "Elevation",
    term == "slope" ~ "Slope",
    term == "aspect_sin" ~ "Aspect sine (N-S)",
    term == "aspect_cos" ~ "Aspect cosine (E-W)",
    term == "TPI" ~ "Topographic Position Index",
    term == "cwd" ~ "Cumulative Water Deficit",
    term == "clay_0_5cm" ~ "Clay content (0-5cm)",
    term == "total_awc" ~ "Total Available Water Content"
  ))

model_estimates_outname <- paste0("02_Tabular_data/B_model_data/results/Estimates_Filtered_", variable, "_allGS_models.csv")
write_csv(results_filter, model_estimates_outname)
