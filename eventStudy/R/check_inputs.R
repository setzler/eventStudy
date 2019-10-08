#' @export
# function to conduct various preliminary checks of long_data and arguments of ES()
# - the first failed check generates an informative error/warning and breaks the program
# - if all checks are successfully passed, function returns NULL, but invisibly
ES_check_inputs <-
  function(long_data,
           outcomevar,
           unit_var,
           cal_time_var,
           onset_time_var,
           cluster_vars,
           omitted_event_time = -2,
           anticipation = 0,
           min_control_gap = 1,
           max_control_gap = Inf,
           linearize_pretrends = FALSE,
           fill_zeros = FALSE,
           residualize_covariates = FALSE,
           control_subset_var = NA,
           control_subset_event_time = 0,
           treated_subset_var = NA,
           treated_subset_event_time = 0,
           control_subset_var2 = NA,
           control_subset_event_time2 = 0,
           treated_subset_var2 = NA,
           treated_subset_event_time2 = 0,
           discrete_covars = NULL,
           cont_covars = NULL,
           never_treat_action = 'none',
           homogeneous_ATT = FALSE,
           reg_weights = NULL,
           add_unit_fes = FALSE,
           bootstrapES = FALSE,
           bootstrap_iters = 1,
           bootstrap_num_cores = 1,
           ipw = FALSE,
           ipw_model = 'linear',
           ipw_composition_change = FALSE,
           ipw_keep_data = FALSE,
           ipw_ps_lower_bound = 0,
           ipw_ps_upper_bound = 1,
           event_vs_noevent = FALSE,
           ref_discrete_covars = NULL,
           ref_discrete_covar_event_time = 0,
           ref_cont_covars = NULL,
           ref_cont_covar_event_time = 0,
           calculate_collapse_estimates = FALSE,
           collapse_inputs = NULL,
           ref_reg_weights = NULL,
           ref_reg_weights_event_time = 0,
           ntile_var = NULL,
           ntile_event_time = -2,
           ntiles = NA,
           ntile_var_value = NA,
           ntile_avg = FALSE,
           endog_var = NULL,
           linearDiD = FALSE,
           linearDiD_treat_var = NULL,
           cohort_by_cohort = FALSE,
           cohort_by_cohort_num_cores = 1) {

    # type checks
    assertDataTable(long_data)
    assertCharacter(outcomevar, len = 1)
    assertCharacter(unit_var, len = 1)
    assertCharacter(cal_time_var, len = 1)
    assertCharacter(onset_time_var, len = 1)
    if (!is.null(cluster_vars)) {
      assertCharacter(cluster_vars)
    }
    assertIntegerish(omitted_event_time, len = 1, upper = -1)
    assertIntegerish(anticipation, len = 1, lower = 0)
    assertIntegerish(min_control_gap, len = 1, lower = 1)
    if (!any(
      testIntegerish(max_control_gap, len = 1, lower = min_control_gap),
      is.infinite(max_control_gap)
    )) {
      assertIntegerish(max_control_gap, len = 1, lower = min_control_gap)
    }
    if (!any(
      testIntegerish(max_control_gap, len = 1, lower = min_control_gap),
      is.infinite(max_control_gap)
    )) {
      assertIntegerish(max_control_gap, len = 1, lower = min_control_gap)
    }
    if (!is.na(control_subset_var)) {
      assertCharacter(control_subset_var, len = 1)
    }
    if (!is.na(treated_subset_var)) {
      assertCharacter(treated_subset_var, len = 1)
    }
    if (!is.na(control_subset_var2)) {
      assertCharacter(control_subset_var2, len = 1)
    }
    if (!is.na(treated_subset_var2)) {
      assertCharacter(treated_subset_var2, len = 1)
    }
    if (!is.null(reg_weights)) {
      assertCharacter(reg_weights, len = 1)
    }
    if (!is.null(ref_reg_weights)) {
      assertCharacter(ref_reg_weights, len = 1)
    }
    assertIntegerish(control_subset_event_time, len = 1)
    assertIntegerish(treated_subset_event_time, len = 1)
    assertIntegerish(control_subset_event_time2, len = 1)
    assertIntegerish(treated_subset_event_time2, len = 1)
    assertIntegerish(ref_discrete_covar_event_time, len = 1)
    assertIntegerish(ref_cont_covar_event_time, len = 1)
    assertIntegerish(ref_reg_weights_event_time, len = 1)
    assertFlag(linearize_pretrends)
    assertFlag(fill_zeros)
    assertFlag(residualize_covariates)
    if (residualize_covariates) {
      if (!any(testCharacter(discrete_covars),
               testCharacter(cont_covars))) {
        stop(
          "Since residualize_covariates=TRUE, either discrete_covars or cont_covars must be provided as a character vector."
        )
      }
    }
    assertFlag(homogeneous_ATT)
    assertFlag(add_unit_fes)
    assertFlag(bootstrapES)
    assertIntegerish(bootstrap_iters, len = 1, lower = 1)
    assertIntegerish(bootstrap_num_cores, len = 1, lower = 1)
    assertFlag(ipw)
    assertFlag(ipw_composition_change)
    assertNumber(ipw_ps_lower_bound, lower = 0, upper = 1)
    assertNumber(ipw_ps_upper_bound, lower = 0, upper = 1)
    assertFlag(calculate_collapse_estimates)
    if (calculate_collapse_estimates == TRUE) {
      assertDataTable(
        collapse_inputs,
        ncols = 2,
        any.missing = FALSE,
        types = c("character", "list")
      )
    }
    if (!is.null(ntile_var)) {
      assertIntegerish(ntile_event_time, len = 1)
      assertIntegerish(ntiles, len = 1, lower = 2)
      assertIntegerish(ntile_var_value,
                       len = 1,
                       lower = 1,
                       upper = ntiles)
      assertFlag(ntile_avg)
    }
    if(!is.null(endog_var)){
      assertCharacter(endog_var,len=1)
    }
    assertFlag(linearDiD)
    if ((!is.null(linearDiD_treat_var) & linearDiD == FALSE)) {
      stop(
        sprintf(
          "Supplied linearDiD_treat_var='%s' but either didn't set linearDiD or set linearDiD='%s'. \n Let me suggest linearDiD=TRUE.",
          linearDiD_treat_var,
          linearDiD
        )
      )
    }
    assertFlag(cohort_by_cohort)
    assertIntegerish(cohort_by_cohort_num_cores,len=1,lower=1)

    # check that anticipation choice and omitted_event_time choice don't conflict
    if (omitted_event_time + anticipation > -1) {
      stop(
        sprintf(
          "omitted_event_time='%s' and anticipation='%s' implies overlap of pre-treatment and anticipation periods. Let me suggest omitted_event_time<='%s'",
          omitted_event_time,
          anticipation,
          ((-1 * anticipation) - 1)
        )
      )
    }

    # check that all of these variables are actually in the data.table, and provide custom error messages.
    # also check that the variable types within long_data are appropriate for the program
    if (!(outcomevar %in% names(long_data))) {
      stop(
        sprintf(
          "Variable outcomevar='%s' is not in the long_data you provided.",
          outcomevar
        )
      )
    }
    if (!(unit_var %in% names(long_data))) {
      stop(sprintf(
        "Variable unit_var='%s' is not in the long_data you provided.",
        unit_var
      ))
    }
    if (!(cal_time_var %in% names(long_data))) {
      stop(
        sprintf(
          "Variable cal_time_var='%s' is not in the long_data you provided.",
          cal_time_var
        )
      )
    }
    if (!(long_data[, typeof(get(cal_time_var))] == "integer")) {
      # we need this as we will take max/min of these when constructing the stacked data in ES_clean_data()
      stop(
        sprintf(
          "Variable cal_time_var='%s' must be of type 'integer' in long_data.",
          cal_time_var
        )
      )
    }
    if (!(onset_time_var %in% names(long_data))) {
      stop(
        sprintf(
          "Variable onset_time_var='%s' is not in the long_data you provided.",
          onset_time_var
        )
      )
    }
    if (!(long_data[, typeof(get(onset_time_var))] == "integer")) {
      # we need this as we will take max/min of these when constructing the stacked data in ES_clean_data()
      stop(
        sprintf(
          "Variable onset_time_var='%s' must be of type 'integer' in long_data.",
          onset_time_var
        )
      )
    }
    if (!is.null(cluster_vars)) {
      for (vv in cluster_vars) {
        if (!(vv %in% names(long_data))) {
          stop(
            sprintf(
              "Variable cluster_vars='%s' is not in the long_data you provided. Let me suggest cluster_vars='%s'.",
              vv,
              unit_var
            )
          )
        }
      }
    }
    if (!is.na(control_subset_var)) {
      if (!(control_subset_var %in% names(long_data))) {
        stop(
          sprintf(
            "Variable control_subset_var='%s' is not in the long_data you provided.",
            control_subset_var
          )
        )
      }
      if (!(long_data[, typeof(get(control_subset_var))] == "logical")) {
        stop(
          sprintf(
            "Variable control_subset_var='%s' must be of type 'logical' in long_data (i.e., only TRUE or FALSE values).",
            control_subset_var
          )
        )
      }
    }
    if (!is.na(treated_subset_var)) {
      if (!(treated_subset_var %in% names(long_data))) {
        stop(
          sprintf(
            "Variable treated_subset_var='%s' is not in the long_data you provided.",
            treated_subset_var
          )
        )
      }
      if (!(long_data[, typeof(get(treated_subset_var))] == "logical")) {
        stop(
          sprintf(
            "Variable treated_subset_var='%s' must be of type 'logical' in long_data (i.e., only TRUE or FALSE values).",
            treated_subset_var
          )
        )
      }
    }
    if (!is.na(control_subset_var2)) {
      if (!(control_subset_var2 %in% names(long_data))) {
        stop(
          sprintf(
            "Variable control_subset_var2='%s' is not in the long_data you provided.",
            control_subset_var2
          )
        )
      }
      if (!(long_data[, typeof(get(control_subset_var2))] == "logical")) {
        stop(
          sprintf(
            "Variable control_subset_var2='%s' must be of type 'logical' in long_data (i.e., only TRUE or FALSE values).",
            control_subset_var2
          )
        )
      }
    }
    if (!is.na(treated_subset_var2)) {
      if (!(treated_subset_var2 %in% names(long_data))) {
        stop(
          sprintf(
            "Variable treated_subset_var2='%s' is not in the long_data you provided.",
            treated_subset_var2
          )
        )
      }
      if (!(long_data[, typeof(get(treated_subset_var2))] == "logical")) {
        stop(
          sprintf(
            "Variable treated_subset_var2='%s' must be of type 'logical' in long_data (i.e., only TRUE or FALSE values).",
            treated_subset_var2
          )
        )
      }
    }
    if (testCharacter(discrete_covars)) {
      for (vv in discrete_covars) {
        if (!(vv %in% names(long_data))) {
          stop(
            sprintf(
              "Variable discrete_covars='%s' is not in the long_data you provided.",
              vv
            )
          )
        }
      }
    }
    if (testCharacter(cont_covars)) {
      for (vv in cont_covars) {
        if (!(vv %in% names(long_data))) {
          stop(sprintf(
            "Variable cont_covars='%s' is not in the long_data you provided.",
            vv
          ))
        }
      }
    }
    if (testCharacter(ref_discrete_covars)) {
      for (vv in ref_discrete_covars) {
        if (!(vv %in% names(long_data))) {
          stop(
            sprintf(
              "Variable ref_discrete_covars='%s' is not in the long_data you provided.",
              vv
            )
          )
        }
      }
    }
    if (testCharacter(ref_cont_covars)) {
      for (vv in ref_cont_covars) {
        if (!(vv %in% names(long_data))) {
          stop(
            sprintf(
              "Variable ref_cont_covars='%s' is not in the long_data you provided.",
              vv
            )
          )
        }
      }
    }
    if (testCharacter(reg_weights)) {
      if (!(reg_weights %in% names(long_data))) {
        stop(
          sprintf(
            "Variable reg_weights='%s' is not in the long_data you provided.",
            reg_weights
          )
        )
      }
    }
    if (testCharacter(ref_reg_weights)) {
      if (!(ref_reg_weights %in% names(long_data))) {
        stop(
          sprintf(
            "Variable ref_reg_weights='%s' is not in the long_data you provided.",
            ref_reg_weights
          )
        )
      }
    }
    if (testCharacter(ntile_var)) {
      if (!(ntile_var %in% names(long_data))) {
        stop(
          sprintf(
            "Variable ntile_var='%s' is not in the long_data you provided.",
            ntile_var
          )
        )
      }
    }
    if(testCharacter(endog_var)){
      if(!(endog_var %in% names(long_data))){stop(sprintf("Variable endog_var='%s' is not in the long_data you provided.",endog_var))}
    }
    if (linearDiD == TRUE) {
      if (testCharacter(linearDiD_treat_var)) {
        if (!(linearDiD_treat_var %in% names(long_data))) {
          stop(
            sprintf(
              "Variable linearDiD_treat_var='%s' is not in the long_data you provided.",
              linearDiD_treat_var
            )
          )
        }
      }
    }

    # temporary check -- code is currently not set up to accept ntile_event_time != omitted_event_time
    # only bites in the case where !is.null(ntile_var) and omitted_event_time != ntile_event_time
    if ((!is.null(ntile_var)) &
        omitted_event_time != ntile_event_time) {
      stop(
        sprintf(
          "ntile_event_time='%s', but currently code can only accept ntile_event_time == omitted_event_time (e.g., %s).",
          ntile_event_time,
          omitted_event_time
        )
      )
    }

    # check that calculate_collapse_estimates == TRUE only if homogeneous_ATT == FALSE
    if (calculate_collapse_estimates == TRUE &
        homogeneous_ATT == TRUE) {
      stop(
        "Cannot have calculate_collapse_estimates == TRUE & homogeneous_ATT == TRUE. Consider setting homogeneous_ATT = FALSE."
      )
    }

    # check that control variables don't overlap with design variables (e.g., cal_time_var, and onset_time_var or unit_var)
    if (add_unit_fes == TRUE) {
      design_vars <- c(cal_time_var, unit_var)
    } else{
      design_vars <- c(cal_time_var, onset_time_var)
    }
    if (testCharacter(discrete_covars)) {
      for (vv in discrete_covars) {
        if (vv %in% design_vars) {
          stop(
            sprintf(
              "Variable discrete_covars='%s' is among c('%s','%s') which are already controlled in the design.",
              vv,
              design_vars[[1]],
              design_vars[[2]]
            )
          )
        }
      }
    }
    if (testCharacter(cont_covars)) {
      for (vv in cont_covars) {
        if (vv %in% design_vars) {
          stop(
            sprintf(
              "Variable cont_covars='%s' is among c('%s','%s') which are already controlled in the design.",
              vv,
              design_vars[[1]],
              design_vars[[2]]
            )
          )
        }
      }
    }
    if (testCharacter(ref_discrete_covars)) {
      for (vv in ref_discrete_covars) {
        if (vv %in% design_vars) {
          stop(
            sprintf(
              "Variable ref_discrete_covars='%s' is among c('%s','%s') which are already controlled in the design.",
              vv,
              design_vars[[1]],
              design_vars[[2]]
            )
          )
        }
      }
    }
    if (testCharacter(ref_cont_covars)) {
      for (vv in ref_cont_covars) {
        if (vv %in% design_vars) {
          stop(
            sprintf(
              "Variable ref_cont_covars='%s' is among c('%s','%s') which are already controlled in the design.",
              vv,
              design_vars[[1]],
              design_vars[[2]]
            )
          )
        }
      }
    }

    # check that user only supplied one of reg_weights / ref_reg_weights, but not both
    if (!is.null(reg_weights) & !is.null(ref_reg_weights)) {
      stop(
        "Supplied variables for both reg_weights and ref_reg_weights, but ES() only admits using one or the other type of weight (or neither)."
      )
    }

    # check that user correctly input what to do with never treated
    if (!(never_treat_action %in% c('none', 'exclude', 'keep', 'only'))) {
      stop(
        sprintf(
          "never_treat_action='%s' is not among allowed values (c('none', 'exclude', 'keep', 'only')).",
          never_treat_action
        )
      )
    }
    if (never_treat_action == 'none' &
        dim(long_data[is.na(get(onset_time_var))])[1] > 0) {
      stop(
        sprintf(
          "never_treat_action='%s' but some units have %s=NA. Please edit supplied long_data or consider another option for never_treat_action.",
          never_treat_action,
          onset_time_var
        )
      )
    }
    if (never_treat_action != 'none' &
        dim(long_data[is.na(get(onset_time_var))])[1] == 0) {
      stop(
        sprintf(
          "never_treat_action='%s' but no units have %s=NA. Let me suggest never_treat_action='none'.",
          never_treat_action,
          onset_time_var
        )
      )
    }

    # warning if cluster_vars = NULL
    if (is.null(cluster_vars)) {
      warning(
        sprintf(
          "Supplied cluster_vars = NULL; given stacking in ES(), standard errors may be too small. Consider cluster_vars='%s' instead.",
          unit_var
        )
      )
    }

    # warning if supplied bootstrap_num_cores*cohort_by_cohort_num_cores exceeds detectCores() - 1
    # not a perfect test (even as an upper bound) as results of detectCores() may be OS-dependent, may not respect cluster allocation limits, etc.
    # see help for detectCores() and mc.cores() for more information
    if(bootstrap_num_cores > 1 & cohort_by_cohort_num_cores == 1){
      if(bootstrap_num_cores > (parallel::detectCores() - 1)){
        warning(sprintf("Supplied bootstrap_num_cores='%s'; this exceeds typical system limits and may cause issues.", bootstrap_num_cores))
      }
    } else if(bootstrap_num_cores == 1 & cohort_by_cohort_num_cores > 1){
      if(cohort_by_cohort_num_cores > (parallel::detectCores() - 1)){
        warning(sprintf("Supplied cohort_by_cohort_num_cores='%s'; this exceeds typical system limits and may cause issues.", cohort_by_cohort_num_cores))
      }
    } else if(bootstrap_num_cores > 1 & cohort_by_cohort_num_cores > 1){
      if(as.integer(bootstrap_num_cores*cohort_by_cohort_num_cores) > (parallel::detectCores() - 1)){
        warning(sprintf("Supplied bootstrap_num_cores='%s' & cohort_by_cohort_num_cores='%s'; the product (%s) exceeds typical system limits and may cause issues.", bootstrap_num_cores, cohort_by_cohort_num_cores, as.integer(bootstrap_num_cores*cohort_by_cohort_num_cores)))
      }
    }

    # check that ipw model conforms to available options, and that propensity score cutoffs are sensible
    if (ipw == TRUE) {
      if (!(ipw_model %in% c('linear', 'logit', 'probit'))) {
        stop(
          sprintf(
            "ipw_model='%s' is not among allowed values (c('linear', 'logit', 'probit')).",
            ipw_model
          )
        )
      }
    }
    if (ipw_ps_lower_bound > ipw_ps_upper_bound) {
      stop(
        sprintf(
          "ipw_ps_lower_bound='%s' & ipw_ps_upper_bound='%s', which means ipw_ps_lower_bound > ipw_ps_upper_bound. Consider revising these cutoffs.",
          ipw_ps_lower_bound,
          ipw_ps_upper_bound
        )
      )
    }

    return(invisible(NULL))
  }
