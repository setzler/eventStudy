
#' @export
ES_clean_data <- function(long_data,
                          outcomevar,
                          unit_var,
                          cal_time_var,
                          onset_time_var,
                          cluster_vars,
                          discrete_covars = NULL,
                          cont_covars = NULL,
                          anticipation = 0,
                          min_control_gap = 1,
                          max_control_gap = Inf,
                          omitted_event_time = -2,
                          treated_subset_var=NA,
                          treated_subset_event_time=NA,
                          control_subset_var=NA,
                          control_subset_event_time=NA,
                          never_treat_action = 'none',
                          never_treat_val = NA,
                          reg_weights = NULL,
                          event_vs_noevent = FALSE) {

  # Restriction based on supplied omitted_event_time

  min_eligible_cohort <- long_data[get(cal_time_var) - get(onset_time_var) == omitted_event_time, min(get(onset_time_var))]
  long_data[, relevant_subset := as.integer(get(onset_time_var) >= min_eligible_cohort)]
  gc()

  # Setting up

  onset_times <- long_data[relevant_subset == 1, sort(unique(get(onset_time_var)))]
  cal_times <- long_data[relevant_subset == 1, sort(unique(get(cal_time_var)))]

  min_onset_time <- min(onset_times)
  max_onset_time <- max(onset_times)

  min_cal_time <- min(cal_times)
  max_cal_time <- max(cal_times)

  # Main code
  j <- 0

  stack_across_cohorts_balanced_treated_control <- list()

  # sequence of treated cohorts based on user choices
  # will start with min_onset_time

  # Last possible treated cohort:
  # a) If there are cohorts treated after the end of the panel, determined by max_cal_time and min_control_gap
  # b) If last onset is before end of panel, determined by max_onset_time, max_cal_time, and min_control_gap
  # c) If there are never-winners, by construction, max_onset_time > max_cal_time irrespective of a) or b) above;
  #     however, only want the min_control_gap to bind to the not-yet treated

  if (max_onset_time > max_cal_time & !(never_treat_action %in% c('keep', 'only'))) {
    last_treat_grp_time <- max_cal_time - (min_control_gap - 1)
  } else if (max_onset_time > max_cal_time & (never_treat_action %in% c('keep', 'only'))) {
    last_treat_grp_time <- max_cal_time
  } else if (max_onset_time <= max_cal_time) {
    last_treat_grp_time <- max_onset_time - min_control_gap
  }

  for (e in intersect(min_onset_time:last_treat_grp_time, onset_times)) {

    # Depending on the structure of the data and choices of min_control_gap / max_control_g, it possible that no control group exists for 'e'
    # E.g., if treated cohorts are 2001, 2002, 2004, 2005 and min_control_gap = max_control_gap = 1, the 2002 has no control group
    # Will skip such a treated cohort at this stage

    count_potential_controls <- dim(long_data[relevant_subset == 1 & between(get(onset_time_var), e + min_control_gap, e + max_control_gap, incbounds = TRUE)])[1]
    if(count_potential_controls <= 0){
      flog.info(sprintf("Given min_control_gap='%s' & max_control_gap='%s', no control units found for %s='%s'",
                        min_control_gap, max_control_gap, onset_time_var, e)
      )
    } else{

      j <- j + 1


      # For a given treated cohort, possible_treated_control is the subset of possible treated and control observations
      possible_treated_control <- list()

      possible_treated_control[[1]] <- long_data[relevant_subset == 1 & get(onset_time_var) == e,
                                                 unique(na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, cluster_vars, discrete_covars, cont_covars, reg_weights))),
                                                 with = FALSE
                                                 ]
      gc()
      possible_treated_control[[1]][, ref_onset_time := e]
      possible_treated_control[[1]][, treated := 1]

      possible_treated_control[[2]] <- long_data[relevant_subset == 1 & between(get(onset_time_var), e + min_control_gap, e + max_control_gap, incbounds = TRUE),
                                                 unique(na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, cluster_vars, discrete_covars, cont_covars, reg_weights))),
                                                 with = FALSE
                                                 ]

      if(never_treat_action == 'only'){
        # will only form control groups using the never-treated units (so excluding the not-yet treated)
        possible_treated_control[[2]] <- possible_treated_control[[2]][get(onset_time_var) == never_treat_val]
      }

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
      # max_control_cohort <- max(possible_treated_control[[onset_time_var]])
      possible_treated_control <- possible_treated_control[(treated == 1)  | ((treated == 0) & (get(cal_time_var) < get(onset_time_var) - anticipation))]

      max_control_year <- possible_treated_control[treated == 0, max(get(cal_time_var))]
      possible_treated_control <- possible_treated_control[get(cal_time_var) <= max_control_year]

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

        # note: the next selection step will have a selection for the treated group, then an "or" |, then the selection for the control group
        # treated group: (get(onset_time_var) == e & ref_event_time %in% c(omitted_event_time, t))
        # treated group: this says onset time must be the current value of e (e is the reference cohort currently under consideration) and then considers only reference event time and outcome event time "t"

        if (t < 1 | event_vs_noevent==TRUE) {
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

      # Now let's get all of these in one stacked data set specific to cohort e
      # Will build up from the pieces of the prior step together with a catt_specific_sample dummy

      balanced_treated_control <- rbindlist(balanced_treated_control, use.names = TRUE)
      balanced_treated_control <- na.omit(balanced_treated_control)
      gc()

      stack_across_cohorts_balanced_treated_control[[j]] <- balanced_treated_control[, unique(na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, cluster_vars, discrete_covars, cont_covars, reg_weights, "ref_onset_time", "ref_event_time", "catt_specific_sample", "treated"))), with = FALSE]
      gc()

      balanced_treated_control <- NULL
      gc()

    }

  }

  # All years and cohorts at once
  stack_across_cohorts_balanced_treated_control <- rbindlist(stack_across_cohorts_balanced_treated_control, use.names = TRUE)
  gc()

  if(!is.na(treated_subset_var) & is.na(control_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := max(as.integer(get(treated_subset_var)*(ref_event_time==treated_subset_event_time))), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[treated==0 | valid_treated_group==1]
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := NULL] # this variable is no longer needed
    gc()
  }

  if(!is.na(control_subset_var) & is.na(treated_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_control_group := max(as.integer(get(control_subset_var)*(ref_event_time==control_subset_event_time))), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[valid_control_group==1  | treated==1 ]
    stack_across_cohorts_balanced_treated_control[, valid_control_group := NULL] # this variable is no longer needed

    gc()
  }

  if(!is.na(control_subset_var) & !is.na(treated_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := max(as.integer(get(treated_subset_var)*(ref_event_time==treated_subset_event_time))), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control[, valid_control_group := max(as.integer(get(control_subset_var)*(ref_event_time==control_subset_event_time))), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[valid_control_group==1  | (treated==1 & valid_treated_group==1) ]
    stack_across_cohorts_balanced_treated_control[, c("valid_treated_group", "valid_control_group") := NULL] # these variables are no longer needed
    gc()
  }

  if(never_treat_action == 'only'){
    rm(never_treat_val)
  }

  long_data[, relevant_subset := NULL]
  gc()

  flog.info(sprintf("Successfully produced a stacked dataset with %s rows.",
                    format(dim(stack_across_cohorts_balanced_treated_control)[1], scientific = FALSE, big.mark = ","))
            )

  return(stack_across_cohorts_balanced_treated_control)
}





ES_parallelize_trends <- function(long_data,
                                  outcomevar,
                                  unit_var,
                                  cal_time_var,
                                  onset_time_var,
                                  anticipation = 0,
                                  reg_weights = NULL) {

  cal_times <- long_data[, sort(unique(get(cal_time_var)))]
  min_cal_time <- min(cal_times)

  start_cols <- copy(colnames(long_data))

  lm_formula_input <- paste(c(sprintf("factor(%s)", cal_time_var), sprintf("factor(%s)*%s", onset_time_var, cal_time_var)), collapse = "+")

  if(!(is.null(reg_weights))){
    est <- lm(as.formula(paste0(eval(outcomevar), " ~ ", lm_formula_input)),
              data = long_data[get(cal_time_var) < (get(onset_time_var) - anticipation)],
              model = FALSE, weights = long_data[get(cal_time_var) < (get(onset_time_var) - anticipation)][[reg_weights]]
    )
  } else{
    est <- lm(as.formula(paste0(eval(outcomevar), " ~ ", lm_formula_input)),
              data = long_data[get(cal_time_var) < (get(onset_time_var) - anticipation)],
              model = FALSE
    )
  }
  gc()

  results <- data.table(rn = names(est$coefficients), estimate = est$coefficients)
  est <- NULL
  gc()
  results <- results[grep("\\:", rn)]
  results[, rn := gsub("\\_", "", rn)]
  results[, rn := gsub(sprintf("\\:%s", gsub("\\_", "", cal_time_var)), "", rn)]
  results[, rn := gsub(sprintf("factor\\(%s\\)", gsub("\\_", "", onset_time_var)), "", rn)]
  results[, rn := as.integer(rn)]
  setnames(results, c("rn", "estimate"), c(onset_time_var, "pre_slope"))

  long_data <- merge(long_data, results, by = onset_time_var, all.x = TRUE)
  long_data[!is.na(pre_slope), (outcomevar) := get(outcomevar) - ((get(cal_time_var) - min_cal_time) * pre_slope)]

  all_added_cols <- setdiff(colnames(long_data), start_cols)
  long_data[, (all_added_cols) := NULL]
  gc()

  setorderv(long_data, c(unit_var, cal_time_var))

  # flog.info("Residualized out cohort-specific pre-trends using only pre-treatment data.")

  return(long_data)
}

ES_residualize_covariates <- function(long_data,
                                      outcomevar,
                                      unit_var,
                                      cal_time_var,
                                      onset_time_var,
                                      discrete_covars = NULL,
                                      cont_covars = NULL,
                                      anticipation = 0,
                                      reg_weights = NULL) {

  start_cols <- copy(colnames(long_data))

  if(!is.null(discrete_covars)){

    # For lm(), want to tell which factor combination will be sent to the intercept
    # Let's set it to the min

    i <- 0
    reference_lookup <- list()
    discrete_covar_formula_input <- c()
    for(var in discrete_covars){
      i <- i + 1
      min_pre_value <- long_data[get(cal_time_var) < (get(onset_time_var) - anticipation), min(get(var))]
      reference_lookup[[i]] <- data.table(varname = var, reference_level = as.character(min_pre_value))

      discrete_covar_formula_input <- c(discrete_covar_formula_input,
                                        sprintf("relevel(as.factor(%s), ref = '%s')", var, min_pre_value))

      var <- NULL
      min_pre_value <- NULL
    }
    reference_lookup <- rbindlist(reference_lookup, use.names = TRUE)

  } else{
    discrete_covar_formula_input = NA
  }

  if(!is.null(cont_covars) & !is.null(discrete_covars)){

    formula_input <- paste(paste0(na.omit(discrete_covar_formula_input), collapse = "+"),
                           paste0(na.omit(cont_covars), collapse = "+"),
                           sep = " + "
                           )

  } else if(is.null(cont_covars) & !is.null(discrete_covars)){

    formula_input <- paste0(paste0(na.omit(discrete_covar_formula_input), collapse = "+"),
                           paste0(na.omit(cont_covars), collapse = "+")
                           )

  } else{

    formula_input <- paste0(paste0(na.omit(cont_covars), collapse = "+"),
                           paste0(na.omit(discrete_covar_formula_input), collapse = "+")
                           )

  }

  if(!(is.null(reg_weights))){
    est <- lm(as.formula(paste0(eval(outcomevar), "~", formula_input)),
              data = long_data[get(cal_time_var) < (get(onset_time_var) - anticipation)],
              model = FALSE, weights = long_data[get(cal_time_var) < (get(onset_time_var) - anticipation)][[reg_weights]]
    )
  } else{
    est <- lm(as.formula(paste0(eval(outcomevar), "~", formula_input)),
              data = long_data[get(cal_time_var) < (get(onset_time_var) - anticipation)],
              model = FALSE
    )
  }
  gc()

  results <- data.table(rn = names(est$coefficients), estimate = est$coefficients)
  est <- NULL
  gc()

  # remove the contribution of any continuous covariates, if any
  if(!is.null(cont_covars)){

    for(var in cont_covars){
      loading = results[rn == var]$estimate
      long_data[, (outcomevar) := get(outcomevar) - (loading*get(var))]
    }
  }

  if(!is.null(discrete_covars)){

    for(var in discrete_covars){
      results[rn == "(Intercept)", (var) := as.integer(reference_lookup[varname == var]$reference_level)]
      append = data.table(estimate = 0)
      append[, (var) := as.integer(reference_lookup[varname == var]$reference_level)]
      results = rbindlist(list(results, append), use.names = TRUE, fill = TRUE)
      results[grep(sprintf("relevel\\(as\\.factor\\(%s\\), ref =", var), rn), (var) := as.integer(gsub('.*[)]', "", rn))]
    }

    # at this stage, want to break apart the intercept and other factor levels
    intercept <- results[rn == "(Intercept)"]$estimate
    results <- results[is.na(rn) | rn != "(Intercept)", c("estimate", discrete_covars), with = FALSE]

    pre_intercepts = c()
    for(var in discrete_covars){
      results_for_merge <- results[, c("estimate", var), with = FALSE]
      setnames(results_for_merge, c("estimate"), sprintf("pre_%s", var))
      results_for_merge <- na.omit(results_for_merge)
      long_data <- merge(long_data, results_for_merge, by = var, all.x = TRUE)
      pre_intercepts = c(pre_intercepts, sprintf("pre_%s", var))
    }

    long_data[, pre_intercept := rowSums(.SD), .SDcols = pre_intercepts]
    long_data <- long_data[!is.na(pre_intercept)]
    # ^ Note: prior step above will restrict the data to only values of discrete_covars found in the pre period

    long_data[, (outcomevar) := get(outcomevar) - pre_intercept - intercept]

  }

  all_added_cols <- setdiff(colnames(long_data), start_cols)
  if(length(all_added_cols) > 0){
    long_data[, (all_added_cols) := NULL]
  }
  gc()

  setorderv(long_data, c(unit_var, cal_time_var))

  # flog.info("Residualized out contribution of covariates using only pre-treatment data.")

  return(long_data)
}

ES_expand_to_balance <- function(long_data,
                                 vars_to_fill = NULL,
                                 unit_var,
                                 cal_time_var,
                                 onset_time_var,
                                 time_invar_covars = NULL
                                 ) {

  # Expand long_data to maximal balanced panel, introducing NAs where needed
  long_data <- setDT(long_data, key = c(unit_var, cal_time_var))[CJ(get(unit_var), get(cal_time_var), unique=TRUE)]

  if(!is.null(vars_to_fill)){
    # replace NAs of provided vars_to_fill columns with zeros
    for (j in vars_to_fill){
      set(long_data, which(is.na(long_data[[j]])), j, 0)
    }
  }

  setorderv(long_data, c(unit_var, cal_time_var))
  gc()

  long_data[, (onset_time_var) := unique(get(onset_time_var)[!is.na(get(onset_time_var))]), list(get(unit_var))]

  if(!is.null(time_invar_covars)){
    # fill in time-invariant values
    long_data[, (time_invar_covars) := lapply(.SD, max, na.rm = TRUE), by = unit_var, .SDcols = time_invar_covars]
  }


  #flog.info("Expanded provided data to balanced panel and filled missings with 0 for vars_to_fill.")

  return(long_data)
}

ES_subset_to_balance <- function(long_data,
                                 unit_var,
                                 cal_time_var
                                 ) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  input_dt <- copy(long_data)

  # find the minimal and maximal calendar times and keep only units that have them all
  years <- length(sort(unique(input_dt[[cal_time_var]])))
  keep_units = input_dt[, .N, by = unit_var][N == years][[unit_var]]
  input_dt <- input_dt[get(unit_var) %in% keep_units]

  setorderv(input_dt, c(unit_var, cal_time_var))
  gc()

  flog.info("Reduced provided data to balanced panel.")

  return(input_dt)
}
