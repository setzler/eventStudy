
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
                          ipw_composition_change = FALSE) {

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

        pr = "pr" # for the eventual na.omit()

        if(ipw_composition_change == TRUE){

          # # First method -- manually weighting
          # # Calculate mean of covariate in the four relevant cells - treat/control and pre/post
          # # Then assign weights relative to a reference cell, separately for within-cell values of time_vary_var
          #
          # means0 = balanced_treated_control[[i]][, list(p = mean(get(ipw_covars_discrete))), by = list(treated, ref_event_time)]
          # means0[, (ipw_covars_discrete) := 0]
          # means0[treated == 1 & ref_event_time != omitted_event_time, weight_old := 1]
          # p_ref = means0[treated == 1 & ref_event_time != omitted_event_time]$p
          # means0[!(treated == 1 & ref_event_time != omitted_event_time), weight_old := (1) / ((p_ref * (1-p))/(p * (1-p_ref)) + 1)]
          # means0[, p := NULL]
          #
          # means1 = balanced_treated_control[[i]][, list(p = mean(get(ipw_covars_discrete))), by = list(treated, ref_event_time)]
          # means1[, (ipw_covars_discrete) := 1]
          # means1[treated == 1 & ref_event_time != omitted_event_time, weight_old := 1]
          # p_ref = means1[treated == 1 & ref_event_time != omitted_event_time]$p
          # means1[!(treated == 1 & ref_event_time != omitted_event_time), weight_old := (p_ref * (1-p))/(p * (1-p_ref)) / ((p_ref * (1-p))/(p * (1-p_ref)) + 1)]
          # means1[, p := NULL]
          #
          # means = rbindlist(list(means0, means1), use.names = TRUE)
          #
          # balanced_treated_control[[i]] = merge(balanced_treated_control[[i]], means, by = c("treated", "ref_event_time", ipw_covars_discrete))
          # gc()

          # Could instead do this through the propensity
          # We have four groups in each CATT-specific sample -- Treated-Pre, Treated-Post, Control-Pre, and Control-Post
          # Recall, reference group was treated == 1 and ref_event_time != omitted_event_time
          # Approach: calculate three propensity scores, and weight relative to reference propensity score

          # Reference group -- Treated-Post
          balanced_treated_control[[i]][, ref_group := as.integer(treated == 1 & ref_event_time != omitted_event_time)]
          balanced_treated_control[[i]][, treated_pre := as.integer(treated == 1 & ref_event_time == omitted_event_time)]
          balanced_treated_control[[i]][, control_post := as.integer(treated == 0 & ref_event_time != omitted_event_time)]
          balanced_treated_control[[i]][, control_pre := as.integer(treated == 0 & ref_event_time == omitted_event_time)]

          formula_input <- paste0(paste0("factor(", na.omit(ipw_covars_discrete), ")", collapse = "+"),
                                  paste(na.omit(ipw_covars_cont), collapse = "+"),
                                  collapse = "+"
                                  )

          balanced_treated_control[[i]][ref_group == 1 | treated_pre == 1, pr := predict(lm(as.formula(paste0("ref_group ~ ", formula_input)),data = balanced_treated_control[[i]][ref_group == 1 | treated_pre == 1]))]
          balanced_treated_control[[i]][ref_group == 1 | control_post == 1, pr := predict(lm(as.formula(paste0("ref_group ~ ", formula_input)),data = balanced_treated_control[[i]][ref_group == 1 | control_post == 1]))]
          balanced_treated_control[[i]][ref_group == 1 | control_pre == 1, pr := predict(lm(as.formula(paste0("ref_group ~ ", formula_input)),data = balanced_treated_control[[i]][ref_group == 1 | control_pre == 1]))]

          balanced_treated_control[[i]][ref_group != 1, weight_unadj := pr / (1 - pr)]
          balanced_treated_control[[i]][ref_group != 1, normalizing_constant := (1 / dim(balanced_treated_control[[i]][ref_group != 1])[1]) * sum(weight_unadj)]
          balanced_treated_control[[i]][ref_group != 1, weight := weight_unadj / normalizing_constant]
          balanced_treated_control[[i]][ref_group == 1, weight := 1]
          balanced_treated_control[[i]][, c("weight_unadj", "normalizing_constant") := NULL]

          # ## TO DELETE
          #
          # pr_check1 <- balanced_treated_control[[i]][, list(mean_age = mean(get(ipw_covars_discrete)),
          #                             w_mean_age = weighted.mean(get(ipw_covars_discrete), w = weight_old)),
          #                      by = list(treated, ref_event_time)
          #                      ]
          # gc()
          # pr_check1 <- dcast(pr_check1, "ref_event_time ~ treated", value.var = c("mean_age", "w_mean_age"))
          # setorderv(pr_check1, c("ref_event_time"))
          #
          # print(pr_check1)
          #
          # pr_check2 <- balanced_treated_control[[i]][, list(mean_age = mean(get(ipw_covars_discrete)),
          #                                                   w_mean_age = weighted.mean(get(ipw_covars_discrete), w = weight)),
          #                                            by = list(treated, ref_event_time)
          #                                            ]
          # gc()
          # pr_check2 <- dcast(pr_check2, "ref_event_time ~ treated", value.var = c("mean_age", "w_mean_age"))
          # setorderv(pr_check2, c("ref_event_time"))
          #
          # print(pr_check2)
          #
          # print(all.equal(pr_check1, pr_check2, tol = 10^-12))
          #
          # check = dim(pr_check2[abs(w_mean_age_0 - w_mean_age_1) > 10^-12])[1]
          # print(ifelse(check > 0, sprintf("Weighted means of %s aren't identical - check propensity scores for errors.", ipw_covars_discrete), sprintf("Weighted means of %s are identical across 4 cells.", ipw_covars_discrete)))
          #
          # ## END TO DELETE

          weight = "weight" # for the eventual na.omit()


        } else if(ipw_composition_change == FALSE){

          formula_input <- paste0(paste0("factor(", na.omit(ipw_covars_discrete), ")", collapse = "+"),
                                  paste(na.omit(ipw_covars_cont), collapse = "+"),
                                  collapse = "+"
          )

          if(ipw_model == "linear"){
            balanced_treated_control[[i]][, pr := predict(lm(as.formula(paste0("treated ~ ", formula_input)),data = balanced_treated_control[[i]]))]
          } else if(ipw_model == "logit"){
            balanced_treated_control[[i]][, pr := predict(glm(as.formula(paste0("treated ~ ", formula_input)),family = binomial(link = "logit"), data = balanced_treated_control[[i]]))]
          } else if(ipw_model == "probit"){
            balanced_treated_control[[i]][, pr := predict(lm(as.formula(paste0("treated ~ ", formula_input)),family = binomial(link = "probit"), data = balanced_treated_control[[i]]))]
          }

          weight = NA # just to make sure the next na.omit() below doesn't break
        }


      } else{
        pr = NA
        weight = NA
      } # just to make sure the next na.omit() below doesn't break

    }

    possible_treated_control <- NULL
    gc()

    # Now let's get all of these in one stacked data set specific to cohort es
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

  flog.info("Residualized out cohort-specific pre-trends using only pre-treatment data.")

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

  input_dt <- copy(long_data)

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

  input_dt <- merge(input_dt, results, by = time_vary_covar)
  input_dt[, (outcomevar) := get(outcomevar) - pre_intercept]
  # ^ Note: merge() above will restrict the data to only values of time_vary_covar found in the pre period

  all_added_cols <- setdiff(colnames(input_dt), start_cols)
  input_dt[, (all_added_cols) := NULL]
  gc()

  setorderv(input_dt, c(unit_var, cal_time_var))

  flog.info("Residualized out contribution of covariates using only pre-treatment data.")

  return(input_dt)
}
