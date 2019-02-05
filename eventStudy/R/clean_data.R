
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
                          control_subset_event_time=NA,
                          ipw = FALSE,
                          ipw_model = "linear",
                          ipw_covars_discrete = NA,
                          ipw_covars_cont = NA,
                          ipw_composition_change = FALSE,
                          never_treat_action,
                          never_treat_val = NA) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  input_dt <- copy(long_data)

  # Restriction based on supplied omitted_event_time

  min_eligible_cohort <- min(input_dt[get(cal_time_var) - get(onset_time_var) == omitted_event_time][[onset_time_var]])
  input_dt <- input_dt[get(onset_time_var) >= min_eligible_cohort]
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
                                              na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, ipw_covars_discrete, ipw_covars_cont)),
                                              with = FALSE
                                              ]
    gc()
    possible_treated_control[[1]][, ref_onset_time := e]
    possible_treated_control[[1]][, treated := 1]

    possible_treated_control[[2]] <- input_dt[between(get(onset_time_var), e + min_control_gap, e + max_control_gap, incbounds = TRUE),
                                              na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, ipw_covars_discrete, ipw_covars_cont)),
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
    max_control_cohort <- max(possible_treated_control[[onset_time_var]])
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

      if(ipw == TRUE){

        balanced_treated_control[[i]] <- ES_make_ipw_dt(did_dt = balanced_treated_control[[i]],
                                                          ipw_composition_change = ipw_composition_change,
                                                          ipw_covars_discrete = ipw_covars_discrete,
                                                          ipw_covars_cont = ipw_covars_cont
                                                          )

      } else{
        pr = NA
        weight = NA
        # just to make sure the next na.omit() below doesn't break
      }

    }

    possible_treated_control <- NULL
    gc()

    # Now let's get all of these in one stacked data set specific to cohort e
    # Will build up from the pieces of the prior step together with a catt_specific_sample dummy

    balanced_treated_control <- rbindlist(balanced_treated_control, use.names = TRUE)
    balanced_treated_control <- na.omit(balanced_treated_control)
    gc()

    stack_across_cohorts_balanced_treated_control[[j]] <- balanced_treated_control[, na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, ipw_covars_discrete, ipw_covars_cont, pr, weight, "ref_onset_time", "ref_event_time", "catt_specific_sample", "treated")), with = FALSE]
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

  if(ipw == TRUE & ipw_composition_change == TRUE){
    flog.info("Estimated three seprate propensity score models per (ref_onset_time, CATT) with Treated-Post as the target population.")
  } else if(ipw == TRUE & ipw_composition_change == FALSE){
    flog.info("Estimated one propensity score model per (ref_onset_time, CATT) with Treated as the target population.")
  }

  if(never_treat_action == 'only'){
    rm(never_treat_val)
  }

  return(stack_across_cohorts_balanced_treated_control)
}

ES_parallelize_trends <- function(long_data,
                                  outcomevar,
                                  unit_var,
                                  cal_time_var,
                                  onset_time_var) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  input_dt <- copy(long_data)

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
  results[, rn := gsub("\\_", "", rn)]
  results[, rn := gsub(sprintf("\\:%s", gsub("\\_", "", cal_time_var)), "", rn)]
  results[, rn := gsub(sprintf("factor\\(%s\\)", gsub("\\_", "", onset_time_var)), "", rn)]
  results[, rn := as.integer(rn)]
  results <- results[, list(rn, Estimate)]
  setnames(results, c("rn", "Estimate"), c(onset_time_var, "pre_slope"))

  input_dt <- merge(input_dt, results, by = onset_time_var, all.x = TRUE)
  input_dt[!is.na(pre_slope), (outcomevar) := get(outcomevar) - (get(cal_time_var) - min_cal_time) * pre_slope]

  all_added_cols <- setdiff(colnames(input_dt), start_cols)
  input_dt[, (all_added_cols) := NULL]
  gc()

  setorderv(input_dt, c(unit_var, cal_time_var))

  # flog.info("Residualized out cohort-specific pre-trends using only pre-treatment data.")

  return(input_dt)
}

ES_residualize_covariates <- function(long_data,
                                              outcomevar,
                                              unit_var,
                                              cal_time_var,
                                              onset_time_var,
                                              discrete_covars = NULL,
                                              cont_covars = NULL) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  input_dt <- copy(long_data)

  start_cols <- copy(colnames(input_dt))

  if(!is.null(discrete_covars)){

    # For lm(), want to tell which factor combination will be sent to the intercept
    # Let's set it to the min

    i <- 0
    reference_lookup <- list()
    discrete_covar_formula_input <- c()
    for(var in discrete_covars){
      i <- i + 1
      min_pre_value <- min(input_dt[get(cal_time_var) < get(onset_time_var)][[var]])
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

    formula_input <- paste(paste0(na.omit(discrete_covar_formula_input), collapse = "+"),
                           paste0(na.omit(cont_covars), collapse = "+"),
                           collapse = " + "
    )

  } else{

    formula_input <- paste(paste0(na.omit(cont_covars), collapse = "+"),
                           paste0(na.omit(discrete_covar_formula_input), collapse = "+"),
                           collapse = " + "
    )

  }

  est <- lm(as.formula(paste0(eval(outcomevar), "~", formula_input)),
            data = input_dt[get(cal_time_var) < get(onset_time_var)],
            model = FALSE
  )
  gc()

  results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
  est <- NULL
  gc()

  # remove the contribution of any continuous covariates, if any
  if(!is.null(cont_covars)){

    for(var in cont_covars){
      loading = results[rn == var]$Estimate
      input_dt[, outcome := outcome - loading*get(var)]
    }
  }

  if(!is.null(discrete_covars)){

    for(var in discrete_covars){
      results[rn == "(Intercept)", (var) := as.integer(reference_lookup[varname == var]$reference_level)]
      append = data.table(Estimate = 0)
      append[, (var) := as.integer(reference_lookup[varname == var]$reference_level)]
      results = rbindlist(list(results, append), use.names = TRUE, fill = TRUE)
      results[grep(sprintf("relevel\\(as\\.factor\\(%s\\), ref =", var), rn), (var) := as.integer(gsub('.*[)]', "", rn))]
    }

    # at this stage, want to break apart the intercept and other factor levels
    intercept <- results[rn == "(Intercept)"]$Estimate
    results <- results[is.na(rn) | rn != "(Intercept)", c("Estimate", discrete_covars), with = FALSE]

    pre_intercepts = c()
    for(var in discrete_covars){
      results_for_merge <- results[, c("Estimate", var), with = FALSE]
      setnames(results_for_merge, c("Estimate"), sprintf("pre_%s", var))
      results_for_merge <- na.omit(results_for_merge)
      input_dt <- merge(input_dt, results_for_merge, by = var, all.x = TRUE)
      pre_intercepts = c(pre_intercepts, sprintf("pre_%s", var))
    }

    input_dt[, pre_intercept := rowSums(.SD), .SDcols = pre_intercepts]
    input_dt <- input_dt[!is.na(pre_intercept)]
    # ^ Note: prior step above will restrict the data to only values of discrete_covars found in the pre period

    input_dt[, (outcomevar) := get(outcomevar) - pre_intercept - intercept]

  }

  all_added_cols <- setdiff(colnames(input_dt), start_cols)
  if(length(all_added_cols) > 0){
    input_dt[, (all_added_cols) := NULL]
  }
  gc()

  setorderv(input_dt, c(unit_var, cal_time_var))

  # flog.info("Residualized out contribution of covariates using only pre-treatment data.")

  return(input_dt)
}

ES_expand_to_balance <- function(long_data,
                                 vars_to_fill = NULL,
                                 unit_var,
                                 cal_time_var,
                                 onset_time_var,
                                 time_invar_covars = NULL
                                 ) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  input_dt <- copy(long_data)

  # Expand input_dt to maximal balanced panel, introducing NAs where needed
  input_dt <- setDT(input_dt, key = c(unit_var, cal_time_var))[CJ(get(unit_var), get(cal_time_var), unique=TRUE)]

  if(!is.null(vars_to_fill)){
    # replace NAs of provided vars_to_fill columns with zeros
    for (j in vars_to_fill){
      set(input_dt, which(is.na(input_dt[[j]])), j, 0)
    }
  }

  setorderv(input_dt, c(unit_var, cal_time_var))
  gc()

  input_dt[, (onset_time_var) := unique(get(onset_time_var)[!is.na(get(onset_time_var))]), list(get(unit_var))]

  if(!is.null(time_invar_covars)){
    # fill in time-invariant values
    input_dt[, (time_invar_covars) := lapply(.SD, max, na.rm = TRUE), by = unit_var, .SDcols = time_invar_covars]
  }


  #flog.info("Expanded provided data to balanced panel and filled missings with 0 for vars_to_fill.")

  return(input_dt)
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

ES_make_ipw_dt <- function(did_dt,
                           ipw_model = "linear",
                           ipw_covars_discrete = NULL,
                           ipw_covars_cont = NULL,
                           ipw_composition_change = FALSE){

  # Just in case
  input_dt <- copy(did_dt)

  if(ipw_composition_change == TRUE){

    # We have four groups in each CATT-specific sample -- Treated-Pre, Treated-Post, Control-Pre, and Control-Post
    # Set reference group was treated == 1 and ref_event_time != omitted_event_time
    # Approach: calculate three propensity scores, and weight relative to reference propensity score

    # Reference group -- Treated-Post

    input_dt[, ref_group := as.integer(treated == 1 & ref_event_time != omitted_event_time)]
    input_dt[, treated_pre := as.integer(treated == 1 & ref_event_time == omitted_event_time)]
    input_dt[, control_post := as.integer(treated == 0 & ref_event_time != omitted_event_time)]
    input_dt[, control_pre := as.integer(treated == 0 & ref_event_time == omitted_event_time)]

    if(!is.null(ipw_covars_discrete)){

      levels <- sort(unique(input_dt[[ipw_covars_discrete]]))
      dropped_levels = NA
      dropped_target_levels = NA

      for(l in levels){

        # Pre-check to establish where common support fails, if at all
        check <- input_dt[, list(mean_age = mean(get(ipw_covars_discrete) == l)),
                          by = list(ref_onset_time, catt_specific_sample, ref_group, treated, ref_event_time)
                          ]
        gc()
        varname_mean <- sprintf("mean_age%s", l)
        setnames(check, c("mean_age"), c(varname_mean))

        check <- check[get(varname_mean) == 0]
        if(dim(check)[1] > 0){

          input_dt <- input_dt[get(ipw_covars_discrete) != l]
          dropped_levels = na.omit(c(dropped_levels, l))

          if(sum(check$ref_group) == 1){
            dropped_target_levels = na.omit(c(dropped_target_levels, l))
          }

        }

      }

      if(!is.na(dropped_levels[1]) & length(dropped_levels) > 0){
        flog.info(paste0(sprintf("Dropped levels of %s for the ref_onset_time %s and ref_event_times (%s, %s) due to lack of common support: %s",
                                 ipw_covars_discrete,
                                 e,
                                 omitted_event_time,
                                 t,
                                 paste0(dropped_levels, collapse = " , ")
        )
        )
        )
      }

      if(!is.na(dropped_target_levels[1]) & length(dropped_target_levels) > 0){
        flog.info(paste0(sprintf("Dropped levels of %s in target population for the ref_onset_time %s and ref_event_times (%s, %s): %s",
                                 ipw_covars_discrete,
                                 e,
                                 omitted_event_time,
                                 t,
                                 paste0(dropped_target_levels, collapse = " , ")
        )
        )
        )
      }
    }

    if(!is.null(cont_covars)){

      formula_input <- paste(paste0("factor(", na.omit(ipw_covars_discrete), ")", collapse = "+"),
                             paste(na.omit(ipw_covars_cont), collapse = "+"),
                             sep = " + "
      )

    } else{

      formula_input <- paste(paste0("factor(", na.omit(ipw_covars_discrete), ")", collapse = "+"),
                             paste(na.omit(ipw_covars_cont), collapse = "+"),
                             collapse = " + "
      )
    }

    input_dt[ref_group == 1 | treated_pre == 1, pr := predict(lm(as.formula(paste0("ref_group ~ ", formula_input)),data = input_dt[ref_group == 1 | treated_pre == 1]))]
    input_dt[ref_group == 1 | control_post == 1, pr := predict(lm(as.formula(paste0("ref_group ~ ", formula_input)),data = input_dt[ref_group == 1 | control_post == 1]))]
    input_dt[ref_group == 1 | control_pre == 1, pr := predict(lm(as.formula(paste0("ref_group ~ ", formula_input)),data = input_dt[ref_group == 1 | control_pre == 1]))]

    input_dt[ref_group != 1, weight_unadj := pr / (1 - pr)]
    input_dt[ref_group != 1, normalizing_constant := (1 / dim(input_dt[ref_group != 1])[1]) * sum(weight_unadj)]
    input_dt[ref_group != 1, weight := weight_unadj / normalizing_constant]
    input_dt[ref_group == 1, weight := 1]
    input_dt[, c("weight_unadj", "normalizing_constant") := NULL]

    pr = "pr" # for the eventual na.omit()
    weight = "weight" # for the eventual na.omit()

  } else if(ipw_composition_change == FALSE){

    formula_input <- paste0(paste0("factor(", na.omit(ipw_covars_discrete), ")", collapse = "+"),
                            paste(na.omit(ipw_covars_cont), collapse = "+"),
                            collapse = "+"
    )

    if(ipw_model == "linear"){
      input_dt[, pr := predict(lm(as.formula(paste0("treated ~ ", formula_input)),data = input_dt))]
    } else if(ipw_model == "logit"){
      input_dt[, pr := predict(glm(as.formula(paste0("treated ~ ", formula_input)),family = binomial(link = "logit"), data = input_dt))]
    } else if(ipw_model == "probit"){
      input_dt[, pr := predict(lm(as.formula(paste0("treated ~ ", formula_input)),family = binomial(link = "probit"), data = input_dt))]
    }

    pr = "pr" # for the eventual na.omit()
    weight = NA # just to make sure the next na.omit() below doesn't break
  }

}
