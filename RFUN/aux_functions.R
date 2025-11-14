#' @title Cumulative Sum of Anomalies
#' @description Computes the cumulative sum of anomalies, optionally scaling and inverting the result.
#' @param anom A numeric vector of anomalies.
#' @param scale A numeric value to scale the anomalies considering x/scale. Example: EVI MODIS anomalies --> x/10000. Defaults = NULL, which means no scaling.
#' @param invert A logical value indicating whether to invert the cumulative sum. Example: MSI anomalies. Defaults = FALSE.
#' @returns A numeric vector of the cumulative sum of anomalies, scaled and/or inverted as specified.
#'

cumsum_anom <- function(anom, scale = NULL, invert = F) {
  if (all(is.na(anom))) {
    return(rep(NA, length(anom)))
  }
  if (is.null(scale)) scale <- 1

  anom_scale <- anom / scale

  anom_acum <- collapse::fcumsum(anom_scale)

  if (invert) {
    anom_acum <- -anom_acum
  }

  anom_acum
}

## -------------------------------------------------------------------------------------#
#' @title Polarity setting for RFD
#' @description Assign positive or negative signs to RFD values based on its anomalies properties.
#' @param rfd_vals A numeric vector of rfd values.
#' @param polarity A numeric vector of the specific signs `-1 = Negative`, `1 = Positive`
#' @returns A numeric vector of signed rfd values.
#'

rfd_polarity <- function(rfd_vals, polarity) {
  if (all(is.na(rfd_vals))) {
    return(rep(NA, length(rfd_vals)))
  }

  rfd_pol <- rfd_vals * polarity

  rfd_pol
}


## -------------------------------------------------------------------------------------#
#' @title Event detection from accumulated anomalies v1
#' @description Detects events based on accumulated anomalies and calculates metrics.
#' @param accum_anom A numeric vector of accumulated anomalies.
#' @param dates A character vector of dates corresponding to the anomalies.
#' @param max_events Maximum number of events to detect. Defaults to 3. Only works when `output = "all"`
#' @param max_gap_days Maximum gap in days allowed between anomalies to consider them part of the same event. Defaults to 60. If there is a gap bigger than 60 days, sign changes are evaluated to see if condition persist or if it can be a new event.
#' @param output Indicates the desired output for the events. Choices are: `"all"` to extract all the events, `"negative"` largest negative eventt, and `"positive"` largest positive.
#' @param rfd_thr Numeric valu used to characterize the proportion of extreme (positive/negative) anomalies inside the detected events. Default is 90.
#' @param min_tud Minimum time under disturbance in days to flag the event. Default is 30.
#' @returns A numeric vector containing the detected events and their metrics. For `output = "all"`, the vector includes:
#' - Total number of positive events
#' - Total number of negative events
#' - Start date of the event
#' - End date of the event
#' - Area under the curve (AUC) of the accumulated anomalies
#' - Max AUC for the specific event
#' - Months to reach the max AUC
#' - Months from max AUC to the change in condition
#' - Polarity of the event (1 for positive, -1 for negative)
#' - Time under disturbance in days
#' - Time under disturbance in months
#' - Id of detection for the specific event
#' # If fewer than `max_events` events are detected, the remaining values will be filled with NA.
#'
#' For `output = "positive"` or `output = "negative"`, the vector consider the biggest event in the specific polarity and its metrics, without the total number of events.


fixed_event_detection <- function(accum_anom,
                                  rfd_vals,
                                  dates,
                                  max_events = 3,
                                  max_gap_days = 60,
                                  rfd_thr = 90,
                                  output = NULL,
                                  min_tud = 30) {
  if (length(accum_anom) != length(dates) || length(accum_anom) != length(rfd_vals)) stop("Dates and anomalies must have the same length.")


  if (output == "all") {
    out_length <- max_events * 13 + 2
  }
  if (output %in% c("positive", "negative")) {
    out_length <- 13
  }


  if (all(is.na(accum_anom)) || (length(unique(accum_anom)) * 100) / length(na.omit(accum_anom)) <= 50) {
    return(as.numeric(rep(NA, out_length)))
  }

  # initial data arrangement
  df <- data.frame(dates = lubridate::ymd(dates), accum = as.numeric(accum_anom), rfd = as.numeric(rfd_vals)) %>%
    na.omit() %>%
    dplyr::arrange(dates) %>%
    mutate(
      day_diff = as.numeric(difftime(dates, dplyr::lag(dates), units = "days")),
      sign = sign(accum),
      sign_lag = dplyr::lag(sign),
      gap = ifelse(is.na(day_diff), 0, day_diff),
      rfd_extreme = ifelse(abs(rfd) >= rfd_thr, 1 * sign(rfd), 0),
      # condition 1: simple sign change
      change_sign = sign != sign_lag & !is.na(sign) & !is.na(sign_lag),

      # condition 2: Na gap + sign change
      change_sign_after_gap = gap > max_gap_days & change_sign,

      # New event grouping conditions
      new_event = change_sign | change_sign_after_gap,

      # group events
      grouped_anom = collapse::fcumsum(collapse::replace_na(new_event, FALSE))
    )

  # Filter values
  events_raw <- df %>%
    dplyr::group_by(grouped_anom) %>%
    dplyr::arrange(dates) %>%
    dplyr::group_split() %>%
    sapply(FUN = as.data.frame, simplify = F)


  if (length(events_raw) == 0) {
    return(as.numeric(rep(NA, out_length)))
  }

  # Calculate event metrics
  events <- purrr::map_dfr(events_raw, function(df_evt) {
    total <- df_evt %>%
      filter(dates >= (base::min(dates) + months(6)) & dates <= (base::max(dates) + months(-6))) %>%
      nrow()
    extreme_pos <- df_evt %>%
      filter(dates >= (base::min(dates) + months(6)) & dates <= (base::max(dates) + months(-6)) & rfd_extreme == 1) %>%
      summarise((n() * 100) / total) %>%
      as.numeric()
    extreme_neg <- df_evt %>%
      filter(dates >= (base::min(dates) + months(6)) & dates <= (base::max(dates) + months(-6)) & rfd_extreme == -1) %>%
      summarise((n() * 100) / total) %>%
      as.numeric()

    out_df <- data.frame(
      start_date = as.integer(base::format(base::min(df_evt$dates), "%Y%m")),
      end_date = as.integer(base::format(base::max(df_evt$dates), "%Y%m")),
      auc = base::sum(abs(df_evt$accum), na.rm = TRUE),
      auc_max = base::max(abs(df_evt$accum), na.rm = TRUE),
      aux_max_date = as.integer(base::format(df_evt$dates[which.max(abs(df_evt$accum))], "%Y%m")),
      months_to_max = lubridate::interval(base::min(df_evt$dates), df_evt$dates[which.max(abs(df_evt$accum))]) %/% months(1) + 1,
      months_since_max = lubridate::interval(df_evt$dates[which.max(abs(df_evt$accum))], base::max(df_evt$dates)) %/% months(1) + 1,
      polarity = sign(df_evt$accum[1]),
      tud_days = as.integer(base::max(df_evt$dates) - base::min(df_evt$dates)) + 1,
      tud_months = lubridate::interval(base::min(df_evt$dates), base::max(df_evt$dates)) %/% months(1) + 1,
      extreme_neg = extreme_neg,
      extreme_pos = extreme_pos
    )

    out_df
  }) %>%
    dplyr::arrange(start_date) %>% # sort
    dplyr::mutate(event_id = dplyr::row_number()) %>% # add detection ID
    dplyr::filter(tud_days >= min_tud) # filter based on duration



  # Take events
  if (output == "all") {
    if (nrow(events) == 0) {
      return(as.numeric(rep(NA, out_length)))
    }

    total_pos <- base::sum(events$polarity == 1)
    total_neg <- base::sum(events$polarity == -1)

    events <- events %>% dplyr::slice_head(n = max_events)
    n_rep <- base::nrow(events)

    # to vector
    vec_events <- events %>% base::t()
    name_vec <- vec_events %>%
      row.names() %>%
      base::rep(times = n_rep)
    vec_events2 <- vec_events %>%
      as.vector() %>%
      round(2)
    # fill missing events
    missing_val <- (max_events * 13) - length(vec_events2)

    if (missing_val > 0) vec_events2 <- c(vec_events2, as.vector(rep(NA, missing_val)))

    # add total events both pos and neg
    out_events <- c(total_pos, total_neg, vec_events2)
    names(out_events) <- c("total_pos", "total_neg", name_vec, rep("not_flagged", missing_val))
  }
  if (output == "positive") {
    events <- events %>%
      filter(polarity == 1) %>%
      dplyr::slice_max(order_by = auc, n = 1)

    if (nrow(events) == 0) {
      return(as.numeric(rep(NA, out_length)))
    }

    vec_events <- events %>% t()
    name_vec <- vec_events %>% row.names()
    out_events <- vec_events %>%
      as.vector() %>%
      round(2)
    names(out_events) <- name_vec
  }
  if (output == "negative") {
    events <- events %>%
      filter(polarity == -1) %>%
      dplyr::slice_max(order_by = auc, n = 1)

    if (nrow(events) == 0) {
      return(as.numeric(rep(NA, out_length)))
    }

    vec_events <- events %>% base::t()
    name_vec <- vec_events %>% row.names()
    out_events <- vec_events %>%
      as.vector() %>%
      round(2)
    names(out_events) <- name_vec
  }

  if (length(out_events) != out_length) {
    stop("different lengths for output")
  }

  return(out_events)
}
