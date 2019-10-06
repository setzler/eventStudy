#' @export
by_cohort_ES <- function(long_data, outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars,
               omitted_event_time= -2, anticipation = 0, min_control_gap=1, max_control_gap=Inf,
               linearize_pretrends=FALSE, fill_zeros=FALSE, residualize_covariates = FALSE,
               control_subset_var=NA, control_subset_event_time=0,
               treated_subset_var=NA, treated_subset_event_time=0,
               control_subset_var2=NA, control_subset_event_time2=0,
               treated_subset_var2=NA, treated_subset_event_time2=0,
               discrete_covars = NULL, cont_covars = NULL, never_treat_action = 'none',
               homogeneous_ATT = FALSE, reg_weights = NULL, add_unit_fes = FALSE,
               bootstrapES = FALSE, bootstrap_iters = 1, bootstrap_num_cores = 1,
               ipw = FALSE, ipw_model = 'linear', ipw_composition_change = FALSE, ipw_keep_data = FALSE, ipw_ps_lower_bound = 0, ipw_ps_upper_bound = 1,
               event_vs_noevent = FALSE,
               ref_discrete_covars = NULL, ref_discrete_covar_event_time=0, ref_cont_covars = NULL, ref_cont_covar_event_time=0,
               calculate_collapse_estimates = FALSE, collapse_inputs = NULL,
               ref_reg_weights = NULL, ref_reg_weights_event_time=0,
               ntile_var = NULL, ntile_event_time = -2, ntiles = NA, ntile_var_value = NA, ntile_avg = FALSE,
               endog_var = NULL,
               linearDiD = FALSE, linearDiD_treat_var = NULL,
               cohort_by_cohort = FALSE, cohort_by_cohort_num_cores = 1){

  flog.info("Beginning ES.")

  # run all preliminary checks
  # TO-DO: all arguments of ES_check_inputs() match ES() inputs, so would be nicer to just pass them all to ES_check_inputs in a less manual way
  ES_check_inputs(long_data = long_data, outcomevar = outcomevar, unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var, cluster_vars = cluster_vars,
                  omitted_event_time = omitted_event_time, anticipation = anticipation, min_control_gap = min_control_gap, max_control_gap = max_control_gap,
                  linearize_pretrends = linearize_pretrends, fill_zeros = fill_zeros, residualize_covariates = residualize_covariates,
                  control_subset_var = control_subset_var, control_subset_event_time = control_subset_event_time,
                  treated_subset_var = treated_subset_var, treated_subset_event_time = treated_subset_event_time,
                  control_subset_var2 = control_subset_var2, control_subset_event_time2 = control_subset_event_time2,
                  treated_subset_var2 = treated_subset_var2, treated_subset_event_time2 = treated_subset_event_time2,
                  discrete_covars = discrete_covars, cont_covars = cont_covars, never_treat_action = never_treat_action,
                  homogeneous_ATT = homogeneous_ATT, reg_weights = reg_weights, add_unit_fes = add_unit_fes,
                  bootstrapES = bootstrapES, bootstrap_iters = bootstrap_iters, bootstrap_num_cores = bootstrap_num_cores,
                  ipw = ipw, ipw_model = ipw_model, ipw_composition_change = ipw_composition_change, ipw_keep_data = ipw_keep_data, ipw_ps_lower_bound = ipw_ps_lower_bound, ipw_ps_upper_bound = ipw_ps_upper_bound,
                  event_vs_noevent = event_vs_noevent,
                  ref_discrete_covars = ref_discrete_covars, ref_discrete_covar_event_time = ref_discrete_covar_event_time,
                  ref_cont_covars = ref_cont_covars, ref_cont_covar_event_time = ref_cont_covar_event_time,
                  calculate_collapse_estimates = calculate_collapse_estimates, collapse_inputs = collapse_inputs,
                  ref_reg_weights = ref_reg_weights, ref_reg_weights_event_time = ref_reg_weights_event_time,
                  ntile_var = ntile_var, ntile_event_time = ntile_event_time, ntiles = ntiles, ntile_var_value = ntile_var_value, ntile_avg = ntile_avg,
                  endog_var = endog_var,
                  linearDiD = linearDiD, linearDiD_treat_var = linearDiD_treat_var,
                  cohort_by_cohort = cohort_by_cohort, cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)

  # if user will be producing bootstrap SEs, want to preserve a copy of the data as it is now;
  # otherwise, potential for changes to underlying sample in subsequent steps
  if(bootstrapES == TRUE){original_sample <- copy(long_data)}

  # edit long_data in line with supplied never_treat_action option
  # NOTE: this will alter supplied long_data
  if(never_treat_action == 'exclude'){
    never_treat_val = NA
    long_data <- long_data[!is.na(get(onset_time_var))]
    gc()
  } else if(never_treat_action %in% c('keep', 'only')){
    # assign the never-treated a unique onset time that ensures they are always part of the control group
    never_treat_val <- max(max(long_data[[onset_time_var]], na.rm = TRUE),
                           max(long_data[[cal_time_var]], na.rm = TRUE)
    ) + min_control_gap + anticipation + 1
    long_data[is.na(get(onset_time_var)), (onset_time_var) := never_treat_val]
  }

  # fill with zeros
  # NOTE: this will alter supplied long_data
  if(fill_zeros){
    flog.info("Filling in zeros.")
    long_data <- ES_expand_to_balance(long_data = long_data,
                                      vars_to_fill = na.omit(unique(c(outcomevar, endog_var))),
                                      unit_var = unit_var,
                                      cal_time_var = cal_time_var,
                                      onset_time_var = onset_time_var)
  }

  # Check that there exist cohorts with observations at omitted_event_time
  if(is.infinite(suppressWarnings(long_data[get(cal_time_var) - get(onset_time_var) == omitted_event_time, min(get(onset_time_var))]))){
    stop(sprintf("Variable onset_time_var='%s' has no treated groups with observations at pre-treatment event time %s.",onset_time_var, omitted_event_time))
  }

  # remove any cases with reg_weights <= 0, -Inf, Inf, or NA (if relevant)
  # NOTE: this will alter supplied long_data
  if(!(is.null(reg_weights))){

    count_long_data_initial <- dim(long_data)[1]
    count_reg_weights_isna <- dim(long_data[is.na(get(reg_weights))])[1]
    count_reg_weights_isinf <- dim(long_data[is.infinite(get(reg_weights))])[1]
    count_reg_weights_nonpositive <- dim(long_data[get(reg_weights) <= 0])[1]

    long_data <- long_data[(!is.na(get(reg_weights))) & (!is.infinite(get(reg_weights))) & (get(reg_weights) > 0)]
    gc()

    count_dropped <- (count_long_data_initial-dim(long_data)[1])
    rm(count_long_data_initial)

    if(count_dropped != 0){
      flog.info((sprintf("\n Warning: Droppped %s from long_data due to missing or extreme reg_weights. \n Of those dropped, the breakdown is: \n 1) %s%% had NA reg_weights \n 2) %s%% had Inf/-Inf reg_weights \n 3) %s%% had reg_weights <= 0",
                         format(count_dropped, scientific = FALSE, big.mark = ","),
                         round(((count_reg_weights_isna / count_dropped) * 100), digits = 4),
                         round(((count_reg_weights_isinf / count_dropped) * 100), digits = 4),
                         round(((count_reg_weights_nonpositive / count_dropped) * 100), digits = 4)
      )
      )
      )
    }
    rm(count_reg_weights_isna, count_reg_weights_isinf, count_reg_weights_nonpositive, count_dropped)
    gc()

  }

  # linearize pre-trends; never-treated will be treated as a single cohort
  # NOTE: this will alter supplied long_data
  if(linearize_pretrends){
    flog.info("Linearizing pre-trends.")
    long_data <- ES_parallelize_trends(long_data = long_data, outcomevar = outcomevar,
                                       unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                                       anticipation = anticipation, reg_weights = reg_weights)
  }

  # extrapolate contribution of covariates using pre-treated period
  # NOTE: this will alter supplied long_data
  if(residualize_covariates){
    flog.info("Residualizing on covariates.")
    long_data <- ES_residualize_covariates(long_data = long_data, outcomevar = outcomevar,
                                           unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                                           anticipation = anticipation, discrete_covars = discrete_covars, cont_covars = cont_covars,
                                           reg_weights = reg_weights)
  }

  # process data
  flog.info("Beginning data stacking.")
  ES_data <- ES_clean_data(long_data = long_data, outcomevar = outcomevar,
                           unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                           anticipation = anticipation, min_control_gap = min_control_gap, max_control_gap = max_control_gap, omitted_event_time = omitted_event_time,
                           control_subset_var = control_subset_var, control_subset_event_time = control_subset_event_time,
                           treated_subset_var = treated_subset_var, treated_subset_event_time = treated_subset_event_time,
                           control_subset_var2 = control_subset_var2, control_subset_event_time2 = control_subset_event_time2,
                           treated_subset_var2 = treated_subset_var2, treated_subset_event_time2 = treated_subset_event_time2,
                           never_treat_action = never_treat_action, never_treat_val = never_treat_val,
                           cluster_vars = cluster_vars, discrete_covars = discrete_covars, cont_covars = cont_covars, reg_weights = reg_weights, event_vs_noevent = event_vs_noevent,
                           ref_discrete_covars = ref_discrete_covars, ref_cont_covars = ref_cont_covars, ref_reg_weights = ref_reg_weights,
                           ntile_var = ntile_var, ntile_event_time = ntile_event_time, ntiles = ntiles, ntile_var_value = ntile_var_value, ntile_avg = ntile_avg,
                           endog_var = endog_var, linearDiD_treat_var = linearDiD_treat_var)

  # construct discrete covariates specific to a ref_event_time
  # will be time invariant for a given ref_onset_time == ref_discrete_covar_event_time, but time-varying across ref_onset_times

  if(!is.null(ref_discrete_covars)){

    start_cols <- copy(colnames(ES_data))

    for(var in ref_discrete_covars){

      # If, for some reason, there is overlap (in variable name) between discrete_covars and ref_discrete_covars, want to make sure not to overwrite
      if(var %in% discrete_covars){

        varname <- sprintf("%s_dyn", var)

        # For a variable that is a character, will need to generate an integer
        # Note: if there are missings, .GRP will assign them all to a single level, so need to turn to them back NAs
        if(class(ES_data[[var]]) == "character"){
          ES_data[, sortorder := .I]
          ES_data[, (varname) := .GRP, keyby = var]
          ES_data[(is.na(get(var)) | get(var) == ""), (varname) := NA]
          setorderv(ES_data, "sortorder")
          ES_data[, sortorder := NULL]
        } else{
          ES_data[, (varname) := get(var)]
        }

      } else{

        varname <- var

        # For a variable that is a character, will need to generate a numeric
        # Note: if there are missings, .GRP will assign them all to a single level, so need to turn to them back NAs
        if(class(ES_data[[var]]) == "character"){
          ES_data[, sortorder := .I]
          ES_data[, tempvar := .GRP, keyby = var]
          ES_data[(is.na(get(var)) | get(var) == ""), tempvar := NA]
          ES_data[, (varname) := NULL]
          setnames(ES_data, "tempvar", varname)
          # ^need to do it this way because (varname) is character and redefining and assigning in one step is a pain
          setorderv(ES_data, "sortorder")
          ES_data[, sortorder := NULL]
        }

      }

      # Define the time-invariant outcome for all units within a ref_onset_time
      # 1) use sum2() instead of max2() in case 'varname' takes on negative values (as max2() would then return 0)
      # 2) if ref_discrete_covar_event_time == omitted_event_time, we need to do sum2() by by=c(unit_var, "ref_onset_time", "catt_specific_sample")
      # ==> otherwise, end up summing undesirably over the repeated rows (as omitted_event_time rows are reused many times)
      # 3) replace values of non-target ref_event_time with NA
      # ==> this ensures that when we sum2() with a by(), we either return the value from the desired ref_event_time or NA
      # ==> otherwise, in case where varname is NA in desired ref_event_time, but not NA at some other times, sum2() would return 0
      ES_data[ref_event_time != ref_discrete_covar_event_time, (varname) := NA]
      if(ref_discrete_covar_event_time == omitted_event_time){
        ES_data[, (varname) := as.integer(sum2(get(varname)*(ref_event_time==ref_discrete_covar_event_time), na.rm = TRUE)), by=c(unit_var, "ref_onset_time", "catt_specific_sample")]
      } else{
        ES_data[, (varname) := as.integer(sum2(get(varname)*(ref_event_time==ref_discrete_covar_event_time), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
      }

    }

    # If any columns were generated above, we can add them to ref_discrete_covars
    # When these are used later in estimation, the loop is over unique values in the intersection of discrete_covars and ref_discrete_covars
    ref_discrete_covars <- unique(na.omit(c(ref_discrete_covars, setdiff(colnames(ES_data), start_cols))))
    rm(start_cols, var, varname)
  }

  # construct continuous covariates specific to a ref_event_time
  # will be time invariant for a given ref_onset_time, but time-varying across ref_onset_times

  if(!is.null(ref_cont_covars)){

    # Will construct this control variable which will be time-invariant, but defined as of the ref_cont_covar_event_time relative to each ref_onset_time
    # Shouldn't have any missings, but in the code, will just treat this as another category

    start_cols <- copy(colnames(ES_data))

    for(var in ref_cont_covars){

      # If, for some reason, there is overlap (in variable name) between cont_covars and ref_cont_covars, want to make sure not to overwrite
      if(var %in% cont_covars){
        varname <- sprintf("%s_dyn", var)
        ES_data[, (varname) := get(var)]
      } else{
        varname <- var
      }

      # Define the time-invariant outcome for all units within a ref_onset_time
      # 1) use sum2() instead of max2() in case 'varname' takes on negative values (as max2() would then return 0)
      # 2) if ref_cont_covar_event_time == omitted_event_time, we need to do sum2() by by=c(unit_var, "ref_onset_time", "catt_specific_sample")
      # ==> otherwise, end up summing undesirably over the repeated rows (as omitted_event_time rows are reused many times)
      # 3) replace values of non-target ref_event_time with NA
      # ==> this ensures that when we sum2() with a by(), we either return the value from the desired ref_event_time or NA
      # ==> otherwise, in case where varname is NA in desired ref_event_time, but not NA at some other times, sum2() would return 0
      ES_data[ref_event_time != ref_cont_covar_event_time, (varname) := NA]
      if(class(ES_data[[var]]) == "integer"){
        if(ref_cont_covar_event_time == omitted_event_time){
          ES_data[, (varname) := as.integer(sum2(get(varname)*(ref_event_time==ref_cont_covar_event_time), na.rm = TRUE)), by=c(unit_var, "ref_onset_time", "catt_specific_sample")]
        } else{
          ES_data[, (varname) := as.integer(sum2(get(varname)*(ref_event_time==ref_cont_covar_event_time), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
        }
      } else{
        if(ref_cont_covar_event_time == omitted_event_time){
          ES_data[, (varname) := as.numeric(sum2(get(varname)*(ref_event_time==ref_cont_covar_event_time), na.rm = TRUE)), by=c(unit_var, "ref_onset_time", "catt_specific_sample")]
        } else{
          ES_data[, (varname) := as.numeric(sum2(get(varname)*(ref_event_time==ref_cont_covar_event_time), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
        }
      }

    }

    # If any columns were generated above, we can add them to ref_cont_covars
    # When these are used later in estimation, the loop is over unique values in the intersection of cont_covars and ref_cont_covars
    ref_cont_covars <- unique(na.omit(c(ref_cont_covars, setdiff(colnames(ES_data), start_cols))))
    rm(start_cols)
  }

  # construct regression weights specific to a ref_event_time
  # will be time invariant for a given ref_onset_time, but time-varying across ref_onset_times

  if(!is.null(ref_reg_weights)){

    # Will construct this weight variable which will be time-invariant, but defined as of the ref_reg_weights_event_time relative to each ref_onset_time
    # Shouldn't have any missings, but in the code, will just treat this as another category

    # Unlike other two types of ref_* vars, no need to check for variable name overlap here with reg_weights
    # as we don't allow using both reg_weight and ref_reg_weights

    # Define the time-invariant outcome for all units within a ref_onset_time
    if(class(ES_data[[ref_reg_weights]]) == "integer"){
      ES_data[, (ref_reg_weights) := as.numeric(get(ref_reg_weights))]
      # needed to include type checks as sum() often converts from integer to numeric
      # use sum() instead of max() just to parallel earlier code; wouldn't anticipate negative values in a weight variable
    }
    ES_data[, (ref_reg_weights) := unique(get(ref_reg_weights)[ref_event_time==ref_reg_weights_event_time], na.rm = TRUE)[1], by=c(unit_var, "ref_onset_time")]


    # remove any cases with ref_reg_weights <= 0, -Inf, Inf, or NA
    count_ES_initial <- dim(ES_data)[1]
    count_ref_reg_weights_isna <- dim(ES_data[is.na(get(ref_reg_weights))])[1]
    count_ref_reg_weights_isinf <- dim(ES_data[is.infinite(get(ref_reg_weights))])[1]
    count_ref_reg_weights_nonpositive <- dim(ES_data[get(ref_reg_weights) <= 0])[1]

    ES_data <- ES_data[(!is.na(get(ref_reg_weights))) & (!is.infinite(get(ref_reg_weights))) & (get(ref_reg_weights) > 0)]
    gc()

    count_dropped <- (count_ES_initial-dim(ES_data)[1])
    rm(count_ES_initial)

    if(count_dropped != 0){
      flog.info((sprintf("\n Warning: Droppped %s from stacked data due to missing or extreme ref_reg_weights. \n Of those dropped, the breakdown is: \n 1) %s%% had NA ref_reg_weights \n 2) %s%% had Inf/-Inf ref_reg_weights \n 3) %s%% had ref_reg_weights <= 0",
                         format(count_dropped, scientific = FALSE, big.mark = ","),
                         round(((count_ref_reg_weights_isna / count_dropped) * 100), digits = 4),
                         round(((count_ref_reg_weights_isinf / count_dropped) * 100), digits = 4),
                         round(((count_ref_reg_weights_nonpositive / count_dropped) * 100), digits = 4)
      )
      )
      )
    }
    rm(count_ref_reg_weights_isna, count_ref_reg_weights_isinf, count_ref_reg_weights_nonpositive, count_dropped)
    gc()
  }

  # estimate inverse probability weights, if relevant
  if(ipw == TRUE){

    # Within each DiD sample, estimate weights to balance provided covariates
    ES_data[, did_id := .GRP, by = list(ref_onset_time, catt_specific_sample)]
    did_ids <- ES_data[, sort(unique(did_id))]

    if(ipw_composition_change == FALSE){

      ES_data[, pr_temp := as.numeric(NA)]

      for(i in did_ids){

        ipw_dt <- ES_make_ipw_dt(did_dt = copy(ES_data[did_id == i]),
                                 unit_var = unit_var,
                                 cal_time_var = cal_time_var,
                                 discrete_covars = discrete_covars,
                                 cont_covars = cont_covars,
                                 ref_discrete_covars = ref_discrete_covars,
                                 ref_cont_covars = ref_cont_covars,
                                 omitted_event_time = omitted_event_time,
                                 ipw_model = ipw_model,
                                 reg_weights = reg_weights,
                                 ref_reg_weights = ref_reg_weights,
                                 ipw_composition_change = ipw_composition_change
        )
        ipw_dt[, did_id := i]

        ES_data <- merge(ES_data, ipw_dt, by = c(unit_var, cal_time_var, "did_id"), all.x = TRUE, sort = FALSE)
        ES_data[is.na(pr_temp) & !is.na(pr), pr_temp := pr]
        ES_data[, pr := NULL]
        ipw_dt <- NULL
        gc()
      }

      setnames(ES_data, "pr_temp", "pr")

      # Restrict to propensity scores strictly in [a,b], a>0 and b<1 (default: a=0, b=1)
      # Also grab the total count, share NA pr, share infinite pr, share pr < a, and share pr > b
      count_ES_initial <- dim(ES_data)[1]
      count_pr_isna <- dim(ES_data[is.na(pr)])[1]
      count_pr_isinf <- dim(ES_data[is.infinite(pr)])[1]
      count_pr_extremelow <- dim(ES_data[pr <= ipw_ps_lower_bound])[1]
      count_pr_extremehigh <- dim(ES_data[pr >= ipw_ps_upper_bound])[1]

      ES_data <- ES_data[!is.na(pr) & !is.infinite(pr) & (between(pr, ipw_ps_lower_bound, ipw_ps_upper_bound, incbounds = FALSE))]
      gc()

      count_dropped <- (count_ES_initial-dim(ES_data)[1])
      rm(count_ES_initial)

      if(count_dropped != 0){
        flog.info((sprintf("\n Warning: Droppped %s from stacked data due to missing or extreme estimated propensity scores. \n Of those dropped, the breakdown is: \n 1) %s%% had a NA propensity score \n 2) %s%% had a Inf/-Inf propensity score \n 3) %s%% had a propensity score <= %s \n 4) %s%% had a propensity score >= %s",
                           format(count_dropped, scientific = FALSE, big.mark = ","),
                           round(((count_pr_isna / count_dropped) * 100), digits = 4),
                           round(((count_pr_isinf / count_dropped) * 100), digits = 4),
                           round(((count_pr_extremelow / count_dropped) * 100), digits = 4),
                           ipw_ps_lower_bound,
                           round(((count_pr_extremehigh / count_dropped) * 100), digits = 4),
                           ipw_ps_upper_bound
        )
        )
        )
      }

      rm(count_pr_isna, count_pr_isinf, count_pr_extremelow, count_pr_extremehigh, count_dropped)

      ES_data[treated == 1, pw := 1]
      ES_data[treated == 0, pw := (pr / (1 - pr))]

      if(ipw_keep_data == TRUE){
        ipw_dt <- ES_data[, list(get(unit_var), ref_onset_time, ref_event_time, catt_specific_sample, treated, pr)]
        setnames(ipw_dt, "V1", unit_var)
      }
    }
  }

  # collect ATT estimates
  if(homogeneous_ATT == FALSE){
    ES_results_hetero <- by_cohort_ES_estimate_ATT(ES_data = copy(ES_data),
                                                   outcomevar=outcomevar,
                                                   unit_var = unit_var,
                                                   onset_time_var = onset_time_var,
                                                   cluster_vars = cluster_vars,
                                                   homogeneous_ATT = homogeneous_ATT,
                                                   omitted_event_time = omitted_event_time,
                                                   discrete_covars = discrete_covars,
                                                   cont_covars = cont_covars,
                                                   ref_discrete_covars = ref_discrete_covars,
                                                   ref_cont_covars = ref_cont_covars,
                                                   residualize_covariates = residualize_covariates,
                                                   reg_weights = reg_weights,
                                                   ref_reg_weights = ref_reg_weights,
                                                   ipw = ipw,
                                                   ipw_composition_change = ipw_composition_change,
                                                   add_unit_fes = add_unit_fes,
                                                   linearDiD = linearDiD,
                                                   linearDiD_treat_var = linearDiD_treat_var,
                                                   cohort_by_cohort = cohort_by_cohort,
                                                   cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)

    gc()

    catt_coefs <- ES_results_hetero[[2]]
    catt_vcov <- ES_results_hetero[[3]]

    ES_results_hetero <- ES_results_hetero[[1]]
    gc()

    setnames(ES_results_hetero,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))
  } else{
    ES_results_hetero <- NULL
  }
  ES_results_homo <- by_cohort_ES_estimate_ATT(ES_data = copy(ES_data),
                                               outcomevar=outcomevar,
                                               unit_var = unit_var,
                                               onset_time_var = onset_time_var,
                                               cluster_vars = cluster_vars,
                                               homogeneous_ATT = TRUE,
                                               omitted_event_time = omitted_event_time,
                                               discrete_covars = discrete_covars,
                                               cont_covars = cont_covars,
                                               ref_discrete_covars = ref_discrete_covars,
                                               ref_cont_covars = ref_cont_covars,
                                               residualize_covariates = residualize_covariates,
                                               reg_weights = reg_weights,
                                               ref_reg_weights = ref_reg_weights,
                                               ipw = ipw,
                                               ipw_composition_change = ipw_composition_change,
                                               add_unit_fes = add_unit_fes,
                                               linearDiD = linearDiD,
                                               linearDiD_treat_var = linearDiD_treat_var,
                                               cohort_by_cohort = FALSE)[[1]]
  gc()
  setnames(ES_results_homo,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

  # collect levels by treatment/control
  if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
    # STILL NEED TO FIX THE SEs to be WEIGHTED; for now use unweighted
    ES_treatcontrol_means <- ES_data[,list(rn="treatment_means", estimate = weighted.mean(get(outcomevar), get(na.omit(unique(c(reg_weights, ref_reg_weights))))), cluster_se = sd(get(outcomevar))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]
  } else{
    ES_treatcontrol_means <- ES_data[,list(rn="treatment_means", estimate = mean(get(outcomevar)), cluster_se = sd(get(outcomevar))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]
  }

  # collect count of treated units by each (ref_onset_time, ref_event_time) for V1 of population-weighted ATTs
  # as we won't have an estimate for the omitted_event_time, exclude it below
  if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
    catt_treated_unique_units <- ES_data[treated == 1 & (ref_event_time != omitted_event_time),list(catt_treated_unique_units = sum(get(na.omit(unique(c(reg_weights, ref_reg_weights)))))), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  } else{
    catt_treated_unique_units <- ES_data[treated == 1 & (ref_event_time != omitted_event_time),list(catt_treated_unique_units = .N), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  }

  # collect count of treated+control units by each (ref_onset_time, ref_event_time) for V2 of population-weighted ATTs
  # as we won't have an estimate for the omitted_event_time, exclude it below
  if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
    catt_total_unique_units <- ES_data[ref_event_time != omitted_event_time,list(catt_total_unique_units = sum(get(na.omit(unique(c(reg_weights, ref_reg_weights)))))), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  } else{
    catt_total_unique_units <- ES_data[ref_event_time != omitted_event_time,list(catt_total_unique_units = .N), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  }

  # collect count of unique units that will be (implicitly) used to estimate cohort- and equally-weighted avgs for each event time
  # as we won't have an estimate for the omitted_event_time, exclude it below
  # will merge these onto figdata at the very end
  event_time_total_unique_units <- ES_data[ref_event_time != omitted_event_time, list(event_time_total_unique_units = uniqueN(get(unit_var), na.rm = TRUE)), by = list(ref_event_time)][order(ref_event_time)]

  # collect count of unique units that will be (implicitly) used to estimate collapsed estimates
  # will merge these onto figdata at the very end (if relevant)
  if(calculate_collapse_estimates == TRUE & homogeneous_ATT == FALSE){

    # quick copy to prevent from editing the supplied collapse_inputs directly
    collapse_input_dt <- copy(collapse_inputs)

    # internal rename columns of collapse_input_dt
    setnames(collapse_input_dt, c("name", "event_times"))

    # for each grouping, get a unique count separately, and then just make a small table
    # as we won't have an estimate for the omitted_event_time, exclude it
    collapsed_estimate_total_unique_units <- list()
    i <- 0
    for(g in unique(na.omit(collapse_input_dt[["name"]]))){
      i <- i + 1
      count <- ES_data[(ref_event_time != omitted_event_time) & (ref_event_time %in% unique(na.omit(unlist(collapse_input_dt[name == g][[2]])))),
                       uniqueN(get(unit_var), na.rm = TRUE)
                       ]
      collapsed_estimate_total_unique_units[[i]] <- data.table(grouping = g, collapsed_estimate_total_unique_units = count)
    }
    rm(i)

    collapsed_estimate_total_unique_units <- rbindlist(collapsed_estimate_total_unique_units, use.names = TRUE)

  }

  # finally, outermost count of unique units in full ES_data
  unique_units <- ES_data[, uniqueN(get(unit_var), na.rm = TRUE)]
  gc()

  ES_data <- NULL
  gc()

  figdata <- rbindlist(list(ES_results_hetero, ES_results_homo, ES_treatcontrol_means), use.names = TRUE, fill=TRUE)

  # merge various *_unique_units into figdata
  # exclude catt_treated_unique_units/catt_total_unique_units as these will be added as part of (un)weighted means below

  # calculate unweighted and V1/V2 weighted means
  catt_treated_unique_units[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  catt_total_unique_units[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  figdata <- merge(figdata, catt_treated_unique_units, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
  figdata <- merge(figdata, catt_total_unique_units, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)

  subsets_for_avgs <- figdata[rn %in% c("catt")]
  subsets_for_avgs[, unweighted_estimate := mean(estimate, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, cohort_weight_V1 := catt_treated_unique_units / sum(catt_treated_unique_units, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, cohort_weight_V2 := catt_total_unique_units / sum(catt_total_unique_units, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weighted_estimate_V1 := weighted.mean(x = estimate, w = cohort_weight_V1, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weighted_estimate_V2 := weighted.mean(x = estimate, w = cohort_weight_V2, na.rm = TRUE), by = list(ref_event_time)]

  # merge weights into figdata
  weights <- subsets_for_avgs[, list(ref_event_time, ref_onset_time, cohort_weight_V1, cohort_weight_V2)]
  figdata <- merge(figdata, weights, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)

  # rbind the (un)weighted avg estimates to figdata much like the "Pooled" data
  subsets_for_avgs[, rowid := seq_len(.N), by = list(ref_event_time)]
  subsets_for_avgs <- subsets_for_avgs[rowid == 1 | is.na(rowid)]

  unweighted <- subsets_for_avgs[, list(ref_event_time, unweighted_estimate)]
  unweighted[, ref_onset_time := "Equally-Weighted"]
  unweighted[, rn := "att"]
  setnames(unweighted, c("unweighted_estimate"), c("estimate"))
  setorderv(unweighted, c("ref_event_time"))

  weighted_V1 <- subsets_for_avgs[, list(ref_event_time, weighted_estimate_V1)]
  weighted_V1[, ref_onset_time := "Cohort-Weighted"]
  weighted_V1[, rn := "att"]
  setnames(weighted_V1, c("weighted_estimate_V1"), c("estimate"))
  setorderv(weighted_V1, c("ref_event_time"))

  weighted_V2 <- subsets_for_avgs[, list(ref_event_time, weighted_estimate_V2)]
  weighted_V2[, ref_onset_time := "Cohort-Weighted V2"]
  weighted_V2[, rn := "att"]
  setnames(weighted_V2, c("weighted_estimate_V2"), c("estimate"))
  setorderv(weighted_V2, c("ref_event_time"))

  figdata <- rbindlist(list(figdata, unweighted, weighted_V1, weighted_V2), use.names = TRUE, fill=TRUE)
  figdata[is.na(cluster_se), cluster_se := 0]

  if(homogeneous_ATT == FALSE){

    # calculate SEs for weighted avg estimates using catt_coefs and catt_vcov
    # the unique "Weighted" ref_event_time values represent the target list of parameters
    # for each such ref_event_time, want to extract the location (number) of the relevant parameters in catt_coefs
    # then will need to grab the relevant weight, and then construct the formula to supply to delta_method()

    event_times <- setdiff(figdata[, sort(unique(ref_event_time))], omitted_event_time)
    onset_times <- as.integer(figdata[rn == "catt", sort(unique(ref_onset_time))])

    min_onset_time <- min(onset_times)
    max_onset_time <- max(onset_times)

    for(et in event_times){

      if(et < 0){
        lookfor <- sprintf("cattlead%s$", abs(et))
        # crucial to have the end-of-line anchor "$" above; otherwise will find, e.g.,  -1 and -19:-10 event times
      } else{
        lookfor <- sprintf("catt%s$", abs(et))
        # crucial to have the end-of-line anchor "$" above; otherwise will find, e.g.,  1 and 10:19 event times
      }
      coef_indices <- grep(lookfor, names(catt_coefs))
      rm(lookfor)
      temp <- as.data.table(do.call(cbind, list(catt_coefs[coef_indices], coef_indices)), keep.rownames = TRUE)
      setnames(temp, c("V1", "V2"), c("estimate", "coef_index"))
      rm(coef_indices)
      temp[, estimate := NULL]
      temp[, rn := gsub("lead", "-", rn)]
      for (c in min_onset_time:max_onset_time) {
        temp[grepl(sprintf("ref\\_onset\\_time%s", c), rn), ref_onset_time := c]
        temp[grepl(sprintf("ref\\_onset\\_time%s", c), rn), rn := gsub(sprintf("ref\\_onset\\_time%s\\_catt", c), "catt", rn)]
      }
      temp[grepl("catt", rn), ref_event_time := as.integer(gsub("catt", "", rn))]
      temp[, rn := NULL]
      temp[, ref_onset_time := as.character(ref_onset_time)]

      # now merge in the weights
      temp <- merge(temp, figdata[rn == "catt"], by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
      temp <- temp[, list(ref_onset_time, ref_event_time, coef_index, cohort_weight_V1, cohort_weight_V2)]
      temp[, equal_weight := 1 / .N]

      temp[, equal_w_formula_entry := sprintf("(%s*x%s)", equal_weight, coef_index)]
      temp[, cohort_w_v1_formula_entry := sprintf("(%s*x%s)", cohort_weight_V1, coef_index)]
      temp[, cohort_w_v2_formula_entry := sprintf("(%s*x%s)", cohort_weight_V2, coef_index)]

      equal_w_g_formula_input = paste0(temp$equal_w_formula_entry, collapse = "+")
      cohort_w_v1_g_formula_input = paste0(temp$cohort_w_v1_formula_entry, collapse = "+")
      cohort_w_v2_g_formula_input = paste0(temp$cohort_w_v2_formula_entry, collapse = "+")

      figdata[rn == "att" & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Equally-Weighted",
              cluster_se := delta_method(g = as.formula(paste("~", equal_w_g_formula_input)),
                                         mean = catt_coefs, cov = catt_vcov, ses = TRUE
              )
              ]

      figdata[rn == "att" & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Cohort-Weighted",
              cluster_se := delta_method(g = as.formula(paste("~", cohort_w_v1_g_formula_input)),
                                         mean = catt_coefs, cov = catt_vcov, ses = TRUE
              )
              ]

      figdata[rn == "att" & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Cohort-Weighted V2",
              cluster_se := delta_method(g = as.formula(paste("~", cohort_w_v2_g_formula_input)),
                                         mean = catt_coefs, cov = catt_vcov, ses = TRUE
              )
              ]

      rm(temp, equal_w_g_formula_input, cohort_w_v1_g_formula_input, cohort_w_v2_g_formula_input)

    }
  }

  rm(subsets_for_avgs)
  rm(weights)
  rm(unweighted)
  rm(weighted_V1)
  rm(weighted_V2)
  gc()

  # merge in the sample sizes
  figdata <- merge(figdata, event_time_total_unique_units, by = "ref_event_time", all.x = TRUE, sort = FALSE)

  # Now we calculate the collapsed estimates, if relevant
  if(calculate_collapse_estimates == TRUE & homogeneous_ATT == FALSE){

    # recall, ran 'collapse_input_dt <- copy(collapse_inputs)' earlier
    # and 'setnames(collapse_input_dt, c("name", "event_times"))

    ew <- list()
    cw <- list()
    cw2 <- list()
    j <- 0

    for(g in unique(na.omit(collapse_input_dt[["name"]]))){

      j <- j + 1

      # extract event_times and results corresponding to grouping
      # as we won't have an estimate for the omitted_event_time, exclude it below
      group_event_times <- setdiff(unique(na.omit(unlist(collapse_input_dt[name == g][[2]]))), omitted_event_time)
      dt <- figdata[(ref_event_time %in% group_event_times) & (rn == "catt")]
      dt[, grouping := g]
      dt[, unweighted_estimate := mean(estimate, na.rm = TRUE), by = list(ref_event_time)]
      dt[, weighted_estimate_V1 := weighted.mean(x = estimate, w = cohort_weight_V1, na.rm = TRUE), by = list(ref_event_time)]
      dt[, weighted_estimate_V2 := weighted.mean(x = estimate, w = cohort_weight_V2, na.rm = TRUE), by = list(ref_event_time)]
      dt[, rowid := seq_len(.N), by = list(ref_event_time)]
      result <- dt[rowid == 1 | is.na(rowid)]
      dt <- dt[, list(ref_event_time, ref_onset_time, cohort_weight_V1, cohort_weight_V2, grouping)]

      ew[[j]] <- result[, list(grouping, unweighted_estimate)]
      ew[[j]][, unweighted_estimate := mean(unweighted_estimate, na.rm = TRUE), by = list(grouping)]
      ew[[j]][, ref_onset_time := "Equally-Weighted + Collapsed"]
      ew[[j]][, rn := "att"]
      ew[[j]][, rowid := seq_len(.N)]
      setnames(ew[[j]], c("unweighted_estimate"), c("estimate"))
      ew[[j]] <- ew[[j]][rowid == 1 | is.na(rowid)]
      ew[[j]][, rowid := NULL]
      ew[[j]][, cluster_se := 0]

      cw[[j]] <- result[, list(grouping, weighted_estimate_V1)]
      cw[[j]][, weighted_estimate_V1 := mean(weighted_estimate_V1, na.rm = TRUE), by = list(grouping)]
      cw[[j]][, ref_onset_time := "Cohort-Weighted + Collapsed"]
      cw[[j]][, rn := "att"]
      cw[[j]][, rowid := seq_len(.N)]
      setnames(cw[[j]], c("weighted_estimate_V1"), c("estimate"))
      cw[[j]] <- cw[[j]][rowid == 1 | is.na(rowid)]
      cw[[j]][, rowid := NULL]
      cw[[j]][, cluster_se := 0]

      cw2[[j]] <- result[, list(grouping, weighted_estimate_V2)]
      cw2[[j]][, weighted_estimate_V2 := mean(weighted_estimate_V2, na.rm = TRUE), by = list(grouping)]
      cw2[[j]][, ref_onset_time := "Cohort-Weighted V2 + Collapsed"]
      cw2[[j]][, rn := "att"]
      cw2[[j]][, rowid := seq_len(.N)]
      setnames(cw2[[j]], c("weighted_estimate_V2"), c("estimate"))
      cw2[[j]] <- cw2[[j]][rowid == 1 | is.na(rowid)]
      cw2[[j]][, rowid := NULL]
      cw2[[j]][, cluster_se := 0]

      rm(result)

      templist = list()
      i = 0
      for(et in group_event_times){

        i = i + 1

        if(et < 0){
          lookfor <- sprintf("cattlead%s$", abs(et))
          # crucial to have the end-of-line anchor "$" above; otherwise will find, e.g.,  -1 and -19:-10 event times
        } else{
          lookfor <- sprintf("catt%s$", abs(et))
          # crucial to have the end-of-line anchor "$" above; otherwise will find, e.g.,  1 and 10:19 event times
        }
        coef_indices <- grep(lookfor, names(catt_coefs))
        rm(lookfor)
        temp <- as.data.table(do.call(cbind, list(catt_coefs[coef_indices], coef_indices)), keep.rownames = TRUE)
        setnames(temp, c("V1", "V2"), c("estimate", "coef_index"))
        rm(coef_indices)
        temp[, estimate := NULL]
        temp[, rn := gsub("lead", "-", rn)]
        for (c in min_onset_time:max_onset_time) {
          temp[grepl(sprintf("ref\\_onset\\_time%s", c), rn), ref_onset_time := c]
          temp[grepl(sprintf("ref\\_onset\\_time%s", c), rn), rn := gsub(sprintf("ref\\_onset\\_time%s\\_catt", c), "catt", rn)]
        }
        temp[grepl("catt", rn), ref_event_time := as.integer(gsub("catt", "", rn))]
        temp[, rn := NULL]
        temp[, ref_onset_time := as.character(ref_onset_time)]

        # now merge in the within-event-time weights
        temp <- merge(temp, dt, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
        temp[, weight_V0 := 1 / .N]

        templist[[i]] <- copy(temp)
        rm(temp)
        gc()

      }
      rm(i)
      rm(dt)

      templist <- rbindlist(templist, use.names = TRUE)

      # Now add the across-event-time weights and calculate full (multiplicative) weights
      templist[, across_weight := (1 / uniqueN(ref_event_time))]
      templist[, full_weight_V0 := weight_V0 * across_weight]
      templist[, full_weight_V1 := cohort_weight_V1 * across_weight]
      templist[, full_weight_V2 := cohort_weight_V2 * across_weight]

      templist[, equal_w_formula_entry := sprintf("(%s*x%s)", full_weight_V0, coef_index)]
      templist[, cohort_w_v1_formula_entry := sprintf("(%s*x%s)", full_weight_V1, coef_index)]
      templist[, cohort_w_v2_formula_entry := sprintf("(%s*x%s)", full_weight_V2, coef_index)]

      formula_input_ew = paste0(templist$equal_w_formula_entry, collapse = "+")
      formula_input_cw = paste0(templist$cohort_w_v1_formula_entry, collapse = "+")
      formula_input_cw2 = paste0(templist$cohort_w_v2_formula_entry, collapse = "+")

      rm(templist)

      ew[[j]][grouping == g,
                   cluster_se := delta_method(g = as.formula(paste("~", formula_input_ew)),
                                              mean = catt_coefs,
                                              cov = catt_vcov,
                                              ses = TRUE
                                              )
                   ]

      cw[[j]][grouping == g,
                   cluster_se := delta_method(g = as.formula(paste("~", formula_input_cw)),
                                              mean = catt_coefs,
                                              cov = catt_vcov,
                                              ses = TRUE
                   )
                   ]

      cw2[[j]][grouping == g,
                   cluster_se := delta_method(g = as.formula(paste("~", formula_input_cw2)),
                                              mean = catt_coefs,
                                              cov = catt_vcov,
                                              ses = TRUE
                   )
                   ]

      rm(formula_input_ew)
      rm(formula_input_cw)
      rm(formula_input_cw2)
    }
    rm(j)
    ew <- rbindlist(ew, use.names = TRUE)
    cw <- rbindlist(cw, use.names = TRUE)
    cw2 <- rbindlist(cw2, use.names = TRUE)

    # merge in sample size for each collapse grouping
    ew <- merge(ew, collapsed_estimate_total_unique_units, by = "grouping", all.x = TRUE, sort = FALSE)
    cw <- merge(cw, collapsed_estimate_total_unique_units, by = "grouping", all.x = TRUE, sort = FALSE)
    cw2 <- merge(cw2, collapsed_estimate_total_unique_units, by = "grouping", all.x = TRUE, sort = FALSE)

    figdata <- rbindlist(list(figdata, ew, cw, cw2), use.names = TRUE, fill = TRUE)
    gc()
  }

  # add overall unique units count
  figdata[, total_unique_units := unique_units]

  # start bootstrap run if relevant
  if(bootstrapES == TRUE){

    # all arguments of bootstrap_ES are passed but for 'iter', which will be supplied by 1:bootstrap_iters below
    boot_results <- rbindlist(parallel::mclapply(X = 1:bootstrap_iters,
                                                 FUN = bootstrap_ES,
                                                 mc.silent = FALSE,
                                                 mc.cores = bootstrap_num_cores,
                                                 mc.set.seed = TRUE,
                                                 long_data = original_sample, outcomevar = outcomevar, unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var, cluster_vars = cluster_vars,
                                                 omitted_event_time= omitted_event_time, anticipation = anticipation, min_control_gap = min_control_gap, max_control_gap = max_control_gap,
                                                 linearize_pretrends = linearize_pretrends, fill_zeros = fill_zeros, residualize_covariates = residualize_covariates,
                                                 control_subset_var = control_subset_var, control_subset_event_time = control_subset_event_time,
                                                 treated_subset_var = treated_subset_var, treated_subset_event_time = treated_subset_event_time,
                                                 control_subset_var2 = control_subset_var2, control_subset_event_time2 = control_subset_event_time2,
                                                 treated_subset_var2 = treated_subset_var2, treated_subset_event_time2 = treated_subset_event_time2,
                                                 discrete_covars = discrete_covars, cont_covars = cont_covars, never_treat_action = never_treat_action,
                                                 homogeneous_ATT = homogeneous_ATT, reg_weights = reg_weights, add_unit_fes = add_unit_fes,
                                                 ipw = ipw, ipw_model = ipw_model, ipw_composition_change = ipw_composition_change, ipw_keep_data = FALSE, ipw_ps_lower_bound = ipw_ps_lower_bound, ipw_ps_upper_bound = ipw_ps_upper_bound,
                                                 event_vs_noevent = event_vs_noevent,
                                                 ref_discrete_covars = ref_discrete_covars, ref_discrete_covar_event_time = ref_discrete_covar_event_time, ref_cont_covars = ref_cont_covars, ref_cont_covar_event_time = ref_cont_covar_event_time,
                                                 calculate_collapse_estimates = calculate_collapse_estimates, collapse_inputs = collapse_inputs,
                                                 ref_reg_weights = ref_reg_weights, ref_reg_weights_event_time = ref_reg_weights_event_time,
                                                 ntile_var = ntile_var, ntile_event_time = ntile_event_time, ntiles = ntiles, ntile_var_value = ntile_var_value, ntile_avg = ntile_avg),
                              use.names = TRUE
                              )

    if(calculate_collapse_estimates == TRUE & homogeneous_ATT == FALSE){
      boot_results_collapse_estimates <- boot_results[!is.na(grouping)]
      boot_results <- boot_results[is.na(grouping)]
      boot_results[, grouping := NULL]
      gc()
    }

    # Construct bootstrap CI for all of the estimates
    boot_results_ses <- boot_results[, list(bootstrap_se = sd(estimate)), by = list(ref_onset_time, ref_event_time)]
    order_to_restore <- na.omit(unique(c(copy(colnames(figdata)),copy(colnames(boot_results_ses)))))
    figdata <- merge(figdata, boot_results_ses, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
    setcolorder(figdata, neworder = order_to_restore)
    rm(order_to_restore)
    figdata[rn == "treatment_means", bootstrap_se := NA]

    if(calculate_collapse_estimates == TRUE & homogeneous_ATT == FALSE){
      boot_results_collapse_estimates_ses <- boot_results_collapse_estimates[, list(bootstrap_se_collapse = sd(estimate)), by = list(ref_onset_time, grouping)]
      order_to_restore <- na.omit(unique(c(copy(colnames(figdata)),copy(colnames(boot_results_collapse_estimates_ses)))))
      figdata <- merge(figdata, boot_results_collapse_estimates_ses, by = c("ref_onset_time", "grouping"), all.x = TRUE, sort = FALSE)
      setcolorder(figdata, neworder = order_to_restore)
      rm(order_to_restore)
      figdata[!is.na(bootstrap_se_collapse), bootstrap_se := bootstrap_se_collapse]
      figdata[, bootstrap_se_collapse := NULL]
    }

    original_sample <- NULL
    gc()

  }

  if(ipw_keep_data == TRUE){figdata <- rbindlist(list(figdata, ipw_dt), use.names = TRUE, fill = TRUE)}

  flog.info('ES is finished.')

  return(figdata)
}
