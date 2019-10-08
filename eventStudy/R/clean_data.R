#' @export
ES_clean_data <- function(long_data,
                          outcomevar,
                          unit_var,
                          cal_time_var,
                          onset_time_var,
                          cluster_vars,
                          discrete_covars = NULL,
                          cont_covars = NULL,
                          ref_discrete_covars = NULL,
                          ref_cont_covars = NULL,
                          anticipation = 0,
                          min_control_gap = 1,
                          max_control_gap = Inf,
                          omitted_event_time = -2,
                          treated_subset_var=NA,
                          treated_subset_event_time=NA,
                          control_subset_var=NA,
                          control_subset_event_time=NA,
                          treated_subset_var2=NA,
                          treated_subset_event_time2=NA,
                          control_subset_var2=NA,
                          control_subset_event_time2=NA,
                          never_treat_action = 'none',
                          never_treat_val = NA,
                          reg_weights = NULL,
                          ref_reg_weights= NULL,
                          event_vs_noevent = FALSE,
                          ntile_var = NULL,
                          ntile_event_time = -2,
                          ntiles = NA,
                          ntile_var_value = NA,
                          ntile_avg = FALSE,
                          endog_var = NULL,
                          linearDiD_treat_var = NULL) {

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

  for (event_cohort in intersect(min_onset_time:last_treat_grp_time, onset_times)) {

    # Depending on the structure of the data and choices of min_control_gap / max_control_g, it possible that no control group exists for 'e'
    # E.g., if treated cohorts are 2001, 2002, 2004, 2005 and min_control_gap = max_control_gap = 1, the 2002 has no control group
    # Will skip such a treated cohort at this stage

    count_potential_controls <- dim(long_data[relevant_subset == 1 & between(get(onset_time_var), event_cohort + min_control_gap, event_cohort + max_control_gap, incbounds = TRUE)])[1]
    if(count_potential_controls <= 0){
      flog.info(sprintf("Given min_control_gap='%s' & max_control_gap='%s', no control units found for %s='%s'",
                        min_control_gap, max_control_gap, onset_time_var, event_cohort)
      )
    } else{

      j <- j + 1


      # For a given treated cohort, possible_treated_control is the subset of possible treated and control observations
      possible_treated_control <- list()

      possible_treated_control[[1]] <- long_data[relevant_subset == 1 & get(onset_time_var) == event_cohort,
                                                 unique(na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, treated_subset_var2, control_subset_var2, cluster_vars, discrete_covars, cont_covars, reg_weights, ref_reg_weights, ref_discrete_covars, ref_cont_covars, ntile_var, endog_var, linearDiD_treat_var))),
                                                 with = FALSE
                                                 ]
      gc()
      possible_treated_control[[1]][, ref_onset_time := event_cohort]
      possible_treated_control[[1]][, treated := 1]

      possible_treated_control[[2]] <- long_data[relevant_subset == 1 & between(get(onset_time_var), event_cohort + min_control_gap, event_cohort + max_control_gap, incbounds = TRUE),
                                                 unique(na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, treated_subset_var2, control_subset_var2, cluster_vars, discrete_covars, cont_covars, reg_weights, ref_reg_weights, ref_discrete_covars, ref_cont_covars, ntile_var, endog_var, linearDiD_treat_var))),
                                                 with = FALSE
                                                 ]

      if(never_treat_action == 'only'){
        # will only form control groups using the never-treated units (so excluding the not-yet treated)
        possible_treated_control[[2]] <- possible_treated_control[[2]][get(onset_time_var) == never_treat_val]
      }

      gc()
      possible_treated_control[[2]][, ref_onset_time := event_cohort]
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
      if(event_vs_noevent==FALSE){
        possible_treated_control <- possible_treated_control[(treated == 1)  | ((treated == 0) & (get(cal_time_var) < get(onset_time_var) - anticipation))]
      }

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

      if(!(is.null(ntile_var)) & ntile_avg == TRUE){

        # Will be using the two-year mean of ntile_var instead of just omitted_event_time
        # 1) use sum2() instead of max2() in case 'ntile_var' takes on negative values (as max2() would then return 0)
        # na.rm = TRUE will mean the missings of ntile_var will be treated as 0
        possible_treated_control[, pre1 := sum2((ref_event_time == omitted_event_time)*get(ntile_var), na.rm = TRUE), by = list(get(unit_var))]
        if((omitted_event_time + anticipation < -1)){
          # case where omitted_event_time is at least 2 years before the start of anticipation/treatment
          # will use the year AFTER omitted_event_time (which stills satisfies the design as its before anticipation)
          possible_treated_control[, pre2 := sum2((ref_event_time == (omitted_event_time + 1))*get(ntile_var), na.rm = TRUE), by = list(get(unit_var))]
        } else if((omitted_event_time + anticipation == -1)){
          # case where omitted_event_time is 1 year before the start of anticipation/treatment
          # will use the year BEFORE omitted_event_time
          possible_treated_control[, pre2 := sum2((ref_event_time == (omitted_event_time - 1))*get(ntile_var), na.rm = TRUE), by = list(get(unit_var))]
        }

        possible_treated_control[, avg_ntile_var := ((pre1 + pre2) / 2)]
        possible_treated_control[, (ntile_var) := avg_ntile_var]
        possible_treated_control[, avg_ntile_var := NULL]
        gc()

      }

      for (time_since_event in setdiff(years, omitted_event_time)) { # excluding focal cohort's omitted_event_time -- recall, panel such that all cohorts have omitted_event_time

        i <- i + 1

        # note: the next selection step will have a selection for the treated group, then an "or" |, then the selection for the control group
        # treated group: (get(onset_time_var) == e & ref_event_time %in% c(omitted_event_time, t))
        # treated group: this says onset time must be the current value of e (e is the reference cohort currently under consideration) and then considers only reference event time and outcome event time "t"

        if (time_since_event < 1 | event_vs_noevent==TRUE) {
          balanced_treated_control[[i]] <- possible_treated_control[(get(onset_time_var) == event_cohort & ref_event_time %in% c(omitted_event_time, time_since_event)) | (get(onset_time_var) > event_cohort & ref_event_time %in% c(omitted_event_time, time_since_event))]
          gc()
        } else if (time_since_event >= 1) {
          balanced_treated_control[[i]] <- possible_treated_control[(get(onset_time_var) == event_cohort & ref_event_time %in% c(omitted_event_time, time_since_event)) | (get(onset_time_var) > (event_cohort + time_since_event) & ref_event_time %in% c(omitted_event_time, time_since_event))]
          gc()
        }
        # Code above ensures that I don't continue using the "omitted_event_time" observations for the control group after it
        # has become treated

        balanced_treated_control[[i]] <- na.omit(balanced_treated_control[[i]],
                                                 cols = unique(na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars, discrete_covars, cont_covars, reg_weights, ntile_var, endog_var, linearDiD_treat_var, "ref_onset_time", "ref_event_time", "treated")))
                                                 )
        gc()
        balanced_treated_control[[i]][, catt_specific_sample := i]

        if(!(is.null(ntile_var))){

          # constructing ntiles using the treated group
          q = quantile(balanced_treated_control[[i]][treated == 1 & ref_event_time == ntile_event_time][[ntile_var]],
                       probs = (seq(1,(ntiles-1), by =1) / ntiles)
          )
          gc()

          # using sum() instead of max() in case ntile_var has negative values
          # na.omit() above means that the na.rm = TRUE should have no bite below
          balanced_treated_control[[i]][, relevant_level := sum(get(ntile_var) * (ref_event_time == ntile_event_time), na.rm = TRUE), by = list(get(unit_var))]

          # loop through a mapping of values of ntile_var to catt_ntile
          # note: order of the loop is deliberate
          balanced_treated_control[[i]][relevant_level > q[(ntiles-1)], catt_ntile := ntiles]
          for(qt in ((ntiles-1):1)){
            balanced_treated_control[[i]][relevant_level <= q[qt], catt_ntile := qt]
          }

          if(!(is.na(ntile_var_value))){
            balanced_treated_control[[i]] <- balanced_treated_control[[i]][catt_ntile == ntile_var_value]
            balanced_treated_control[[i]][, catt_ntile := NULL]
          }
          balanced_treated_control[[i]][, relevant_level := NULL]
          rm(q)
          gc()

        }

      }

      possible_treated_control <- NULL
      gc()

      # Now let's get all of these in one stacked data set specific to cohort e
      # Will build up from the pieces of the prior step together with a catt_specific_sample dummy

      balanced_treated_control <- rbindlist(balanced_treated_control, use.names = TRUE)
      balanced_treated_control <- na.omit(balanced_treated_control,
                                          cols = unique(na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars, discrete_covars, cont_covars, reg_weights, ntile_var, endog_var, linearDiD_treat_var, "ref_onset_time", "ref_event_time", "catt_specific_sample", "treated")))
                                         )
      gc()

      stack_across_cohorts_balanced_treated_control[[j]] <- balanced_treated_control[, unique(na.omit(c(outcomevar, unit_var, cal_time_var, onset_time_var, treated_subset_var, control_subset_var, treated_subset_var2, control_subset_var2, cluster_vars, discrete_covars, cont_covars, reg_weights, ref_reg_weights, ref_discrete_covars, ref_cont_covars, ntile_var, endog_var, linearDiD_treat_var, "ref_onset_time", "ref_event_time", "catt_specific_sample", "treated"))), with = FALSE]
      gc()

      balanced_treated_control <- NULL
      gc()

    }

  }

  # All years and cohorts at once
  stack_across_cohorts_balanced_treated_control <- rbindlist(stack_across_cohorts_balanced_treated_control, use.names = TRUE)
  gc()

  # Restricting stacked data to selection of interest

  if(!is.na(treated_subset_var) & is.na(control_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[treated==0 | valid_treated_group==1]
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := NULL] # this variable is no longer needed
    gc()
  }

  if(!is.na(control_subset_var) & is.na(treated_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[valid_control_group==1  | treated==1 ]
    stack_across_cohorts_balanced_treated_control[, valid_control_group := NULL] # this variable is no longer needed
    gc()
  }

  if(!is.na(control_subset_var) & !is.na(treated_subset_var)){
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[(treated==0 & valid_control_group==1)  | (treated==1 & valid_treated_group==1) ]
    stack_across_cohorts_balanced_treated_control[, c("valid_treated_group", "valid_control_group") := NULL] # these variables are no longer needed
    gc()
  }

  if(!is.na(treated_subset_var2) & is.na(control_subset_var2)){
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := as.integer(max2((get(treated_subset_var2)*(ref_event_time==treated_subset_event_time2)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[treated==0 | valid_treated_group==1]
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := NULL] # this variable is no longer needed
    gc()
  }

  if(!is.na(control_subset_var2) & is.na(treated_subset_var2)){
    stack_across_cohorts_balanced_treated_control[, valid_control_group := as.integer(max2((get(control_subset_var2)*(ref_event_time==control_subset_event_time2)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[valid_control_group==1  | treated==1 ]
    stack_across_cohorts_balanced_treated_control[, valid_control_group := NULL] # this variable is no longer needed
    gc()
  }

  if(!is.na(control_subset_var2) & !is.na(treated_subset_var2)){
    stack_across_cohorts_balanced_treated_control[, valid_treated_group := as.integer(max2((get(treated_subset_var2)*(ref_event_time==treated_subset_event_time2)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control[, valid_control_group := as.integer(max2((get(control_subset_var2)*(ref_event_time==control_subset_event_time2)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
    stack_across_cohorts_balanced_treated_control <- stack_across_cohorts_balanced_treated_control[(treated==0 & valid_control_group==1)  | (treated==1 & valid_treated_group==1) ]
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
