
#' @export
ES_clean_data <- function(long_data,
                          outcomevar,
                          unit_var,
                          cal_time_var,
                          onset_time_var,
                          min_control_gap = 1,
                          max_control_gap = Inf,
                          omitted_event_time = -2,
                          treated_subset_var=NA,
                          treated_subset_event_time=NA,
                          control_subset_var=NA,
                          control_subset_event_time=NA) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  input_dt = copy(long_data)

  # Restriction based on supplied omitted_event_time

  min_eligible_cohort = min(input_dt[get(cal_time_var) - get(onset_time_var) == omitted_event_time][[onset_time_var]])
  input_dt = input_dt[get(onset_time_var) >= min_eligible_cohort]
  gc()

  # Setting up

  onset_times <- input_dt[, sort(unique(get(onset_time_var)))]
  cal_times <- input_dt[, sort(unique(get(cal_time_var)))]

  min_onset_time <- min(onset_times)
  max_onset_time <- max(onset_times)

  min_cal_time <- min(cal_times)
  max_cal_time <- max(cal_times)

  # Main code
  j <- 0

  stack_across_cohorts_balanced_treated_control <- list()

  # sequence of treated cohorts based on user choices
  # will start with min_onset_time

  # last possible treated cohort:
  # a) if there are cohorts treated after the end of the panel, determined by max_cal_time and min_control_gap
  # b) if last onset is before end of panel, determined by max_onset_time, max_cal_time, and min_control_gap

  if (max_onset_time > max_cal_time) {
    last_treat_grp_time <- max_cal_time - (min_control_gap - 1)
  } else if (max_onset_time <= max_cal_time) {
    last_treat_grp_time <- max_onset_time - min_control_gap
  }

  for (e in min_onset_time:last_treat_grp_time) {
    j <- j + 1


    # For a given treated cohort, possible_treated_control is the subset of possible treated and control observations
    possible_treated_control <- list()

    possible_treated_control[[1]] <- input_dt[get(onset_time_var) == e,
                                              na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var)),
                                              with = FALSE
                                              ]
    gc()
    possible_treated_control[[1]][, ref_onset_time := e]
    possible_treated_control[[1]][, treated := 1]

    possible_treated_control[[2]] <- input_dt[between(get(onset_time_var), e + min_control_gap, e + max_control_gap, incbounds = TRUE),
                                              na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var)),
                                              with = FALSE]
    gc()
    possible_treated_control[[2]][, ref_onset_time := e]
    possible_treated_control[[2]][, treated := 0]

    possible_treated_control <- rbindlist(possible_treated_control, use.names = TRUE)
    gc()

    possible_treated_control[, ref_event_time := get(cal_time_var) - ref_onset_time]

    # # Key step -- making sure to only use control groups pre-treatment
    # possible_treated_control <- possible_treated_control[treated == 1 | treated == 0 & get(cal_time_var) < get(onset_time_var)]
    # gc()

    # # Key step -- making sure to only use control groups pre-treatment
    # possible_treated_control <- possible_treated_control[treated == 1 | treated == 0 & get(cal_time_var) < get(onset_time_var)- (min_control_gap - 1)]
    # gc()

    # # Key step -- making sure to only use control groups pre-treatment
    # possible_treated_control <- possible_treated_control[treated == 1 & get(cal_time_var) <= min(max_cal_time, max_onset_time) - (min_control_gap - 1) | treated == 0 & get(cal_time_var) < get(onset_time_var)- (min_control_gap - 1)]
    # gc()

    # Key step -- making sure to only use control groups pre-treatment and treated groups where there are control observations
    max_control_cohort = max(possible_treated_control[[onset_time_var]])
    possible_treated_control <- possible_treated_control[treated == 1 & get(cal_time_var) < max_control_cohort - (min_control_gap - 1) | treated == 0 & get(cal_time_var) < get(onset_time_var)- (min_control_gap - 1)]
    gc()

    i <- 0

    # For a given cohort-specific ATT, balanced_treated_control is the subset of possible treated and control observations
    # that are valid in the sense that the same cohorts are used for both the pre and post differences that form the DiD.
    balanced_treated_control <- list()

    temp <- possible_treated_control[, .N, by = "ref_event_time"]
    years <- sort(unique(temp$ref_event_time))
    temp <- NULL
    gc()

    for (t in setdiff(years, omitted_event_time)) { # excluding focal cohort's omitted_event_time -- recall, panel such that all cohorts have omitted_event_time

      i <- i + 1

      if (t < 1) {
        balanced_treated_control[[i]] <- possible_treated_control[(get(onset_time_var) == e & ref_event_time %in% c(omitted_event_time, t)) | (get(onset_time_var) > e & ref_event_time %in% c(omitted_event_time, t))]
        gc()
      } else if (t >= 1) {
        balanced_treated_control[[i]] <- possible_treated_control[(get(onset_time_var) == e & ref_event_time %in% c(omitted_event_time, t)) | (get(onset_time_var) > (e + t) & ref_event_time %in% c(omitted_event_time, t))]
        gc()
      }
      # Code above ensures that I don't continue using the "omitted_event_time" observations for the control group after it
      # has become treated

      balanced_treated_control[[i]] <- na.omit(balanced_treated_control[[i]])
      gc()
      balanced_treated_control[[i]][, catt_specific_sample := i]
    }

    possible_treated_control <- NULL
    gc()

    # Now let's get all of these in one stacked data set specific to cohort es
    # Will build up from the pieces of the prior step together with a catt_specific_sample dummy

    balanced_treated_control <- rbindlist(balanced_treated_control, use.names = TRUE)
    balanced_treated_control <- na.omit(balanced_treated_control)
    gc()

    stack_across_cohorts_balanced_treated_control[[j]] <- balanced_treated_control[, na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var ,"ref_onset_time", "ref_event_time", "catt_specific_sample", "treated")), with = FALSE]
    gc()

    balanced_treated_control <- NULL
    gc()
  }

  # All years and cohorts at once
  stack_across_cohorts_balanced_treated_control <- rbindlist(stack_across_cohorts_balanced_treated_control, use.names = TRUE)
  gc()

  if(!is.na(treated_subset_var) & is.na(control_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := max(as.integer(get(treated_subset_var)*(ref_event_time==treated_subset_event_time))), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[treated==0 | valid_treated_group==1]
    gc()
  }

  if(!is.na(control_subset_var) & is.na(treated_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_control_group := max(as.integer(get(control_subset_var)*(ref_event_time==control_subset_event_time))), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[valid_control_group==1  | treated==1 ]
    gc()
  }

  if(!is.na(control_subset_var) & !is.na(treated_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := max(as.integer(get(treated_subset_var)*(ref_event_time==treated_subset_event_time))), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control[, valid_control_group := max(as.integer(get(control_subset_var)*(ref_event_time==control_subset_event_time))), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[valid_control_group==1  | (treated==1 & valid_treated_group==1) ]
    gc()
  }

  flog.info("Successfully produced a stacked dataset.")
  return(stack_across_cohorts_balanced_treated_control)
}

ES_parallelize_trends <- function(long_data,
                                  outcomevar,
                                  unit_var,
                                  cal_time_var,
                                  onset_time_var) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  input_dt = copy(long_data)

  cal_times <- input_dt[, sort(unique(get(cal_time_var)))]
  min_cal_time <- min(cal_times)

  start_cols <- copy(colnames(input_dt))

  lm_formula_input <- paste(c(sprintf("factor(%s)", cal_time_var), sprintf("factor(%s)*%s", onset_time_var, cal_time_var)), collapse = "+")

  est <- lm(as.formula(paste0(eval(outcomevar), " ~ ", lm_formula_input)),
    data = input_dt[get(cal_time_var) < get(onset_time_var)]
  )
  gc()

  results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
  est <- NULL
  gc()
  results <- results[grep("\\:", rn)]
  results[, rn := gsub("\\:tax\\_yr", "", rn)]
  results[, rn := gsub("factor\\(win\\_yr\\)", "", rn)]
  results[, rn := as.integer(rn)]
  results <- results[, list(rn, Estimate)]
  setnames(results, c("rn", "Estimate"), c(onset_time_var, "pre_slope"))

  input_dt <- merge(input_dt, results, by = onset_time_var, all.x = TRUE)
  input_dt[!is.na(pre_slope), (outcomevar) := get(outcomevar) - (get(cal_time_var) - min_cal_time) * pre_slope]

  all_added_cols <- setdiff(colnames(input_dt), start_cols)
  input_dt[, (all_added_cols) := NULL]
  gc()

  setorderv(input_dt, c(unit_var, cal_time_var))
  return(input_dt)
}

ES_residualize_time_varying_covar <- function(long_data,
                                              outcomevar,
                                              unit_var,
                                              cal_time_var,
                                              onset_time_var,
                                              time_vary_covar) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  input_dt = copy(long_data)

  cal_times <- input_dt[, sort(unique(get(cal_time_var)))]
  min_cal_time <- min(cal_times)

  start_cols <- copy(colnames(input_dt))

  est <- felm(as.formula(paste0(eval(outcomevar), " ~ -1 + ", sprintf("factor(%s)", time_vary_covar), " | 0 | 0 | 0")),
              data = input_dt[get(cal_time_var) < get(onset_time_var)],
              nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
  )
  gc()

  results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
  est <- NULL
  gc()
  results[, rn := gsub("\\_", "", rn)]
  results[, rn := gsub(sprintf("factor\\(%s\\)", gsub("\\_", "", time_vary_covar)), "", rn)]
  results[, rn := as.integer(rn)]
  results <- results[, list(rn, Estimate)]
  setnames(results, c("rn", "Estimate"), c(time_vary_covar, "pre_intercept"))

  input_dt <- merge(input_dt, results, by = time_vary_covar, all.x = TRUE)
  input_dt[!is.na(pre_intercept), (outcomevar) := get(outcomevar) - pre_intercept]

  all_added_cols <- setdiff(colnames(input_dt), start_cols)
  input_dt[, (all_added_cols) := NULL]
  gc()

  setorderv(input_dt, c(unit_var, cal_time_var))
  return(input_dt)
}
