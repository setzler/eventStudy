#' @export
Wald_ES <- function(long_data, outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars,
               omitted_event_time= -2, anticipation = 0, min_control_gap=1, max_control_gap=Inf,
               linearize_pretrends=FALSE, fill_zeros=FALSE, residualize_covariates = FALSE,
               control_subset_var=NA, control_subset_event_time=0,
               treated_subset_var=NA, treated_subset_event_time=0,
               control_subset_var2=NA, control_subset_event_time2=0,
               treated_subset_var2=NA, treated_subset_event_time2=0,
               discrete_covars = NULL, cont_covars = NULL, never_treat_action = 'none',
               homogeneous_ATT = FALSE, reg_weights = NULL, add_unit_fes = FALSE,
               bootstrapES = FALSE, bootstrap_iters = 1, bootstrap_num_cores = 1, keep_all_bootstrap_results = FALSE,
               ipw = FALSE, ipw_model = 'linear', ipw_composition_change = FALSE, ipw_keep_data = FALSE, ipw_ps_lower_bound = 0, ipw_ps_upper_bound = 1,
               event_vs_noevent = FALSE,
               ref_discrete_covars = NULL, ref_discrete_covar_event_time=0, ref_cont_covars = NULL, ref_cont_covar_event_time=0,
               calculate_collapse_estimates = FALSE, collapse_inputs = NULL,
               ref_reg_weights = NULL, ref_reg_weights_event_time=0,
               ntile_var = NULL, ntile_event_time = -2, ntiles = NA, ntile_var_value = NA, ntile_avg = FALSE,
               endog_var = NULL, cohort_by_cohort = FALSE, cohort_by_cohort_num_cores = 1){

  flog.info("Beginning Wald_ES.")

  # type checks
  assertDataTable(long_data)
  assertCharacter(outcomevar,len=1)
  assertCharacter(unit_var,len=1)
  assertCharacter(cal_time_var,len=1)
  if(!is.null(cluster_vars)){assertCharacter(cluster_vars)}
  assertIntegerish(omitted_event_time,len=1,upper=-1)
  assertIntegerish(anticipation,len=1,lower=0)
  assertIntegerish(min_control_gap,len=1,lower=1)
  if(!any(testIntegerish(max_control_gap,len=1,lower=min_control_gap),is.infinite(max_control_gap))){
    assertIntegerish(max_control_gap,len=1,lower=min_control_gap)
  }
  if(!any(testIntegerish(max_control_gap,len=1,lower=min_control_gap),is.infinite(max_control_gap))){
    assertIntegerish(max_control_gap,len=1,lower=min_control_gap)
  }
  if(!is.na(control_subset_var)){
    assertCharacter(control_subset_var,len=1)
  }
  if(!is.na(treated_subset_var)){
    assertCharacter(treated_subset_var,len=1)
  }
  if(!is.na(control_subset_var2)){
    assertCharacter(control_subset_var2,len=1)
  }
  if(!is.na(treated_subset_var2)){
    assertCharacter(treated_subset_var2,len=1)
  }
  if(!is.null(reg_weights)){
    assertCharacter(reg_weights,len=1)
  }
  if(!is.null(ref_reg_weights)){
    assertCharacter(ref_reg_weights,len=1)
  }
  assertIntegerish(control_subset_event_time,len=1)
  assertIntegerish(treated_subset_event_time,len=1)
  assertIntegerish(control_subset_event_time2,len=1)
  assertIntegerish(treated_subset_event_time2,len=1)
  assertIntegerish(ref_discrete_covar_event_time,len=1)
  assertIntegerish(ref_cont_covar_event_time,len=1)
  assertIntegerish(ref_reg_weights_event_time,len=1)
  assertFlag(linearize_pretrends)
  assertFlag(fill_zeros)
  assertFlag(residualize_covariates)
  if(residualize_covariates){
    if(!any(testCharacter(discrete_covars), testCharacter(cont_covars))){
      stop("Since residualize_covariates=TRUE, either discrete_covars or cont_covars must be provided as a character vector.")
    }
  }
  assertFlag(homogeneous_ATT)
  assertFlag(add_unit_fes)
  assertFlag(bootstrapES)
  assertIntegerish(bootstrap_iters,len=1,lower=1)
  assertIntegerish(bootstrap_num_cores,len=1,lower=1)
  assertFlag(keep_all_bootstrap_results)
  assertFlag(ipw)
  assertFlag(ipw_composition_change)
  assertNumber(ipw_ps_lower_bound,lower=0,upper=1)
  assertNumber(ipw_ps_upper_bound,lower=0,upper=1)
  assertFlag(calculate_collapse_estimates)
  if(calculate_collapse_estimates == TRUE){
    assertDataTable(collapse_inputs, ncols = 2, any.missing = FALSE, types = c("character", "list"))
  }
  if(!is.null(ntile_var)){
    assertIntegerish(ntile_event_time,len=1)
    assertIntegerish(ntiles,len=1,lower=2)
    assertIntegerish(ntile_var_value,len=1,lower=1,upper=ntiles)
    assertFlag(ntile_avg)
  }
  if(!is.null(endog_var)){
    assertCharacter(endog_var,len=1)
  }
  assertFlag(cohort_by_cohort)
  assertIntegerish(cohort_by_cohort_num_cores,len=1,lower=1)

  # check that anticipation choice and omitted_event_time choice don't conflict
  if(omitted_event_time + anticipation > -1){
    stop(sprintf("omitted_event_time='%s' and anticipation='%s' implies overlap of pre-treatment and anticipation periods. Let me suggest omitted_event_time<='%s'",
                 omitted_event_time, anticipation, ((-1 * anticipation) - 1)))
  }

  # check that all of these variables are actually in the data.table, and provide custom error messages.
  if(!(outcomevar %in% names(long_data))){stop(sprintf("Variable outcomevar='%s' is not in the long_data you provided.",outcomevar))}
  if(!(unit_var %in% names(long_data))){stop(sprintf("Variable unit_var='%s' is not in the long_data you provided.",unit_var))}
  if(!(cal_time_var %in% names(long_data))){stop(sprintf("Variable cal_time_var='%s' is not in the long_data you provided.",cal_time_var))}
  if(!(onset_time_var %in% names(long_data))){stop(sprintf("Variable onset_time_var='%s' is not in the long_data you provided.",onset_time_var))}
  if(!is.null(cluster_vars)){
    for(vv in cluster_vars){
      if(!(vv %in% names(long_data))){stop(sprintf("Variable cluster_vars='%s' is not in the long_data you provided. Let me suggest cluster_vars='%s'.",vv,unit_var))}
    }
  }
  if(!is.na(control_subset_var)){
    if(!(control_subset_var %in% names(long_data))){stop(sprintf("Variable control_subset_var='%s' is not in the long_data you provided.",control_subset_var))}
    if(!(long_data[,typeof(get(control_subset_var))]=="logical")){stop(sprintf("Variable control_subset_var='%s' must be of type logical (i.e., only TRUE or FALSE values).",control_subset_var))}
  }
  if(!is.na(treated_subset_var)){
    if(!(treated_subset_var %in% names(long_data))){stop(sprintf("Variable treated_subset_var='%s' is not in the long_data you provided.",treated_subset_var))}
    if(!(long_data[,typeof(get(treated_subset_var))]=="logical")){stop(sprintf("Variable treated_subset_var='%s' must be of type logical (i.e., only TRUE or FALSE values).",treated_subset_var))}
  }
  if(!is.na(control_subset_var2)){
    if(!(control_subset_var2 %in% names(long_data))){stop(sprintf("Variable control_subset_var2='%s' is not in the long_data you provided.",control_subset_var2))}
    if(!(long_data[,typeof(get(control_subset_var2))]=="logical")){stop(sprintf("Variable control_subset_var2='%s' must be of type logical (i.e., only TRUE or FALSE values).",control_subset_var2))}
  }
  if(!is.na(treated_subset_var2)){
    if(!(treated_subset_var2 %in% names(long_data))){stop(sprintf("Variable treated_subset_var2='%s' is not in the long_data you provided.",treated_subset_var2))}
    if(!(long_data[,typeof(get(treated_subset_var2))]=="logical")){stop(sprintf("Variable treated_subset_var2='%s' must be of type logical (i.e., only TRUE or FALSE values).",treated_subset_var2))}
  }
  if(testCharacter(discrete_covars)){
    for(vv in discrete_covars){
      if(!(vv %in% names(long_data))){stop(sprintf("Variable discrete_covars='%s' is not in the long_data you provided.",vv))}
    }
  }
  if(testCharacter(cont_covars)){
    for(vv in cont_covars){
      if(!(vv %in% names(long_data))){stop(sprintf("Variable cont_covars='%s' is not in the long_data you provided.",vv))}
    }
  }
  if(testCharacter(ref_discrete_covars)){
    for(vv in ref_discrete_covars){
      if(!(vv %in% names(long_data))){stop(sprintf("Variable ref_discrete_covars='%s' is not in the long_data you provided.",vv))}
    }
  }
  if(testCharacter(ref_cont_covars)){
    for(vv in ref_cont_covars){
      if(!(vv %in% names(long_data))){stop(sprintf("Variable ref_cont_covars='%s' is not in the long_data you provided.",vv))}
    }
  }
  if(testCharacter(reg_weights)){
    if(!(reg_weights %in% names(long_data))){stop(sprintf("Variable reg_weights='%s' is not in the long_data you provided.",reg_weights))}
  }
  if(testCharacter(ref_reg_weights)){
    if(!(ref_reg_weights %in% names(long_data))){stop(sprintf("Variable ref_reg_weights='%s' is not in the long_data you provided.",ref_reg_weights))}
  }
  if(testCharacter(ntile_var)){
    if(!(ntile_var %in% names(long_data))){stop(sprintf("Variable ntile_var='%s' is not in the long_data you provided.",ntile_var))}
  }
  if(testCharacter(endog_var)){
    if(!(endog_var %in% names(long_data))){stop(sprintf("Variable endog_var='%s' is not in the long_data you provided.",endog_var))}
  }

  # temporary check -- code is currently not set up to accept ntile_event_time != omitted_event_time
  # only bites in the case where !is.null(ntile_var) and omitted_event_time != ntile_event_time
  if((!is.null(ntile_var)) & omitted_event_time != ntile_event_time){
    stop(sprintf("ntile_event_time='%s', but currently code can only accept ntile_event_time == omitted_event_time (e.g., %s).", ntile_event_time, omitted_event_time))
  }

  # check that calculate_collapse_estimates == TRUE only if homogeneous_ATT == FALSE
  if(calculate_collapse_estimates == TRUE & homogeneous_ATT == TRUE){
    stop("Cannot have calculate_collapse_estimates == TRUE & homogeneous_ATT == TRUE. Consider setting homogeneous_ATT = FALSE.")
  }

  # check that control variables don't overlap with design variables (e.g., cal_time_var, and onset_time_var or unit_var)
  if(add_unit_fes == TRUE){
    design_vars <- c(cal_time_var, unit_var)
  } else{
    design_vars <- c(cal_time_var, onset_time_var)
  }
  if(testCharacter(discrete_covars)){
    for(vv in discrete_covars){
      if(vv %in% design_vars){stop(sprintf("Variable discrete_covars='%s' is among c('%s','%s') which are already controlled in the design.",vv, design_vars[[1]], design_vars[[2]]))}
    }
  }
  if(testCharacter(cont_covars)){
    for(vv in cont_covars){
      if(vv %in% design_vars){stop(sprintf("Variable cont_covars='%s' is among c('%s','%s') which are already controlled in the design.",vv, design_vars[[1]], design_vars[[2]]))}
    }
  }
  if(testCharacter(ref_discrete_covars)){
    for(vv in ref_discrete_covars){
      if(vv %in% design_vars){stop(sprintf("Variable ref_discrete_covars='%s' is among c('%s','%s') which are already controlled in the design.",vv, design_vars[[1]], design_vars[[2]]))}
    }
  }
  if(testCharacter(ref_cont_covars)){
    for(vv in ref_cont_covars){
      if(vv %in% design_vars){stop(sprintf("Variable ref_cont_covars='%s' is among c('%s','%s') which are already controlled in the design.",vv, design_vars[[1]], design_vars[[2]]))}
    }
  }

  # check that user only supplied one of reg_weights / ref_reg_weights, but not both
  if(!is.null(reg_weights) & !is.null(ref_reg_weights)){
    stop("Supplied variables for both reg_weights and ref_reg_weights, but Wald_ES() only admits using one or the other type of weight (or neither).")
  }

  # check that user correctly input what to do with never treated
  if(!(never_treat_action %in% c('none', 'exclude', 'keep', 'only'))){
    stop(sprintf("never_treat_action='%s' is not among allowed values (c('none', 'exclude', 'keep', 'only')).", never_treat_action))
  }
  if(never_treat_action=='none' & dim(long_data[is.na(get(onset_time_var))])[1] > 0){
    stop(sprintf("never_treat_action='%s' but some units have %s=NA. Please edit supplied long_data or consider another option for never_treat_action.", never_treat_action, onset_time_var))
  }
  if(never_treat_action!='none' & dim(long_data[is.na(get(onset_time_var))])[1] == 0){
    stop(sprintf("never_treat_action='%s' but no units have %s=NA. Let me suggest never_treat_action='none'.", never_treat_action, onset_time_var))
  }

  # warning if cluster_vars = NULL
  if(is.null(cluster_vars)){warning(sprintf("Supplied cluster_vars = NULL; given stacking in Wald_ES(), standard errors may be too small. Consider cluster_vars='%s' instead.", unit_var))}

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
  if(ipw == TRUE){
    if(!(ipw_model %in% c('linear', 'logit', 'probit'))){
      stop(sprintf("ipw_model='%s' is not among allowed values (c('linear', 'logit', 'probit')).", ipw_model))
    }
  }
  if(ipw_ps_lower_bound > ipw_ps_upper_bound){
    stop(sprintf("ipw_ps_lower_bound='%s' & ipw_ps_upper_bound='%s', which means ipw_ps_lower_bound > ipw_ps_upper_bound. Consider revising these cutoffs.", ipw_ps_lower_bound, ipw_ps_upper_bound))
  }

  # if user will be producing bootstrap SEs, want to preserve a copy of the data as it is now;
  # otherwise, potential for changes to underlying sample in subsequent steps
  if(bootstrapES == TRUE){
    original_sample <- copy(long_data)
  }

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

    flog.info((sprintf("\n Warning: Droppped %s from long_data due to missing or extreme reg_weights. \n Of those dropped, the breakdown is: \n 1) %s%% had NA reg_weights \n 2) %s%% had Inf/-Inf reg_weights \n 3) %s%% had reg_weights <= 0",
                       format(count_dropped, scientific = FALSE, big.mark = ","),
                       round(((count_reg_weights_isna / count_dropped) * 100), digits = 4),
                       round(((count_reg_weights_isinf / count_dropped) * 100), digits = 4),
                       round(((count_reg_weights_nonpositive / count_dropped) * 100), digits = 4)
    )
    )
    )
    rm(count_reg_weights_isna, count_reg_weights_isinf, count_reg_weights_nonpositive, count_dropped)
    gc()

  }

  # linearize pre-trends; never-treated will be treated as a single cohort
  # NOTE: this will alter supplied long_data
  if(linearize_pretrends){
    flog.info("NOT CURRENTLY Linearizing pre-trends.")
  }

  # extrapolate contribution of covariates using pre-treated period
  # NOTE: this will alter supplied long_data
  if(residualize_covariates){
    flog.info("NOT CURRENTLY Residualizing on covariates.")
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
                           endog_var = endog_var)

  # construct discrete covariates specific to a ref_event_time
  # will be time invariant for a given ref_onset_time == ref_discrete_covar_event_time, but time-varying across ref_onset_times

  if(!is.null(ref_discrete_covars)){

    start_cols <- copy(colnames(ES_data))

    for(var in ref_discrete_covars){

      # If, for some reason, there is overlap (in variable name) between discrete_covars and ref_discrete_covars, want to make sure not to overwrite
      if(var %in% discrete_covars){

        varname <- sprintf("%s_dyn", var)

        # For a variable that is a character, will need to generate a numeric
        # Note: if there are missings, this will assign them all to a single level
        if(class(ES_data[[var]]) == "character"){
          ES_data[, (varname) := .GRP, by = var]
        } else{
          ES_data[, (varname) := get(var)]
        }

        # Define the time-invariant outcome for all units within a ref_onset_time
        # use sum() instead of max() in case 'varname' takes on negative values (as max() would then return 0)
        ES_data[, (varname) := sum(as.integer(get(varname)*(ref_event_time==ref_discrete_covar_event_time))), by=c(unit_var, "ref_onset_time")]

      } else{

        varname <- var

        # For a variable that is a character, will need to generate a numeric
        if(class(ES_data[[var]]) == "character"){
          ES_data[, varnum := .GRP, by = var]
          ES_data[, (varname) := NULL]
          setnames(ES_data, "varnum", varname)
          # ^need to do it this way because (varname) is character and redefining and assigning in one step is a pain
        }

        # Define the time-invariant outcome for all units within a ref_onset_time
        # use sum() instead of max() in case 'varname' takes on negative values (as max() would then return 0)
        ES_data[, (varname) := sum(as.integer(get(varname)*(ref_event_time==ref_discrete_covar_event_time))), by=c(unit_var, "ref_onset_time")]

      }

    }

    # If any columns were generated above, we can add them to ref_discrete_covars
    # When these are used later in estimation, the loop is over unique values in the intersection of discrete_covars and ref_discrete_covars
    ref_discrete_covars <- unique(na.omit(c(ref_discrete_covars, setdiff(colnames(ES_data), start_cols))))
    rm(start_cols)
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

        # Define the time-invariant outcome for all units within a ref_onset_time
        if(class(ES_data[[var]]) == "integer"){
          # needed to include type checks as sum() often converts from integer to numeric
          # use sum() instead of max() in case 'var' takes on negative values (as max() would then return 0)
          ES_data[, (varname) := sum(as.integer(get(var)*(ref_event_time==ref_cont_covar_event_time))), by=c(unit_var, "ref_onset_time")]
        } else{
          ES_data[, (varname) := sum(get(var)*(ref_event_time==ref_cont_covar_event_time)), by=c(unit_var, "ref_onset_time")]
        }
      } else{

        # Define the time-invariant outcome for all units within a ref_onset_time
        if(class(ES_data[[var]]) == "integer"){
          # needed to include type checks as sum() often converts from integer to numeric
          # use sum() instead of max() in case 'var' takes on negative values (as max() would then return 0)
          ES_data[, (var) := sum(as.integer(get(var)*(ref_event_time==ref_cont_covar_event_time))), by=c(unit_var, "ref_onset_time")]
        } else{
          ES_data[, (var) := sum(get(var)*(ref_event_time==ref_cont_covar_event_time)), by=c(unit_var, "ref_onset_time")]
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
      # needed to include type checks as sum() often converts from integer to numeric
      # use sum() instead of max() just to parallel earlier code; wouldn't anticipate negative values in a weight variable
      ES_data[, (ref_reg_weights) := sum(as.integer(get(ref_reg_weights)*(ref_event_time==ref_reg_weights_event_time))), by=c(unit_var, "ref_onset_time")]
    } else{
      ES_data[, (ref_reg_weights) := sum(get(ref_reg_weights)*(ref_event_time==ref_reg_weights_event_time)), by=c(unit_var, "ref_onset_time")]
    }

    # remove any cases with ref_reg_weights <= 0, -Inf, Inf, or NA
    count_ES_initial <- dim(ES_data)[1]
    count_ref_reg_weights_isna <- dim(ES_data[is.na(get(ref_reg_weights))])[1]
    count_ref_reg_weights_isinf <- dim(ES_data[is.infinite(get(ref_reg_weights))])[1]
    count_ref_reg_weights_nonpositive <- dim(ES_data[get(ref_reg_weights) <= 0])[1]

    ES_data <- ES_data[(!is.na(get(ref_reg_weights))) & (!is.infinite(get(ref_reg_weights))) & (get(ref_reg_weights) > 0)]
    gc()

    count_dropped <- (count_ES_initial-dim(ES_data)[1])
    rm(count_ES_initial)

    flog.info((sprintf("\n Warning: Droppped %s from stacked data due to missing or extreme ref_reg_weights. \n Of those dropped, the breakdown is: \n 1) %s%% had NA ref_reg_weights \n 2) %s%% had Inf/-Inf ref_reg_weights \n 3) %s%% had ref_reg_weights <= 0",
                       format(count_dropped, scientific = FALSE, big.mark = ","),
                       round(((count_ref_reg_weights_isna / count_dropped) * 100), digits = 4),
                       round(((count_ref_reg_weights_isinf / count_dropped) * 100), digits = 4),
                       round(((count_ref_reg_weights_nonpositive / count_dropped) * 100), digits = 4)
    )
    )
    )
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

    # reduced form; homogeneous_ATT = FALSE
    # NOTE: need to copy(ES_data) as by_cohort_ES_estimate_ATT changes the underlying data
    ES_results_hetero_rf <- by_cohort_ES_estimate_ATT(ES_data = copy(ES_data),
                                                      outcomevar = outcomevar,
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
                                                      cohort_by_cohort = cohort_by_cohort,
                                                      cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)
    gc()

    catt_coefs_rf <- ES_results_hetero_rf[[2]]
    catt_vcov_rf <- ES_results_hetero_rf[[3]]

    ES_results_hetero_rf <- ES_results_hetero_rf[[1]]
    gc()

    setnames(ES_results_hetero_rf,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

    if(!is.null(endog_var)){
      # first_stage; homogeneous_ATT = FALSE
      # NOTE: need to copy(ES_data) as by_cohort_ES_estimate_ATT changes the underlying data
      ES_results_hetero_fs <- by_cohort_ES_estimate_ATT(ES_data = copy(ES_data),
                                                        outcomevar = endog_var,
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
                                                        cohort_by_cohort = cohort_by_cohort,
                                                        cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)
      gc()

      catt_coefs_fs <- ES_results_hetero_fs[[2]]
      catt_vcov_fs <- ES_results_hetero_fs[[3]]

      ES_results_hetero_fs <- ES_results_hetero_fs[[1]]
      gc()

      setnames(ES_results_hetero_fs,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

    } else{
      ES_results_hetero_fs <- NULL
    }
  } else{
    ES_results_hetero_rf <- NULL
    ES_results_hetero_fs <- NULL
  }

  # reduced form; homogeneous_ATT = TRUE
  # NOTE: need to copy(ES_data) as by_cohort_ES_estimate_ATT changes the underlying data
  ES_results_homo_rf <- by_cohort_ES_estimate_ATT(ES_data = copy(ES_data),
                                                  outcomevar = outcomevar,
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
                                                  cohort_by_cohort = cohort_by_cohort,
                                                  cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)[[1]]
  gc()
  setnames(ES_results_homo_rf,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

  if(!is.null(endog_var)){
    # first_stage; homogeneous_ATT = TRUE
    # NOTE: need to copy(ES_data) as by_cohort_ES_estimate_ATT changes the underlying data
    ES_results_homo_fs <- by_cohort_ES_estimate_ATT(ES_data = copy(ES_data),
                                                    outcomevar = endog_var,
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
                                                    cohort_by_cohort = cohort_by_cohort,
                                                    cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)[[1]]
    gc()
    setnames(ES_results_homo_fs,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))
  } else{
    ES_results_homo_fs <- NULL
  }

  # collect levels by treatment/control
  if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
    # STILL NEED TO FIX THE SEs to be WEIGHTED; for now use unweighted
    ES_treatcontrol_means_rf <- ES_data[,list(rn="treatment_means", estimate = weighted.mean(get(outcomevar), get(na.omit(unique(c(reg_weights, ref_reg_weights))))), cluster_se = sd(get(outcomevar))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]
    if(!is.null(endog_var)){
      ES_treatcontrol_means_fs <- ES_data[,list(rn="treatment_means", estimate = weighted.mean(get(endog_var), get(na.omit(unique(c(reg_weights, ref_reg_weights))))), cluster_se = sd(get(endog_var))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]
    } else{
      ES_treatcontrol_means_fs <- NULL
    }
  } else{
    ES_treatcontrol_means_rf <- ES_data[,list(rn="treatment_means", estimate = mean(get(outcomevar)), cluster_se = sd(get(outcomevar))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]
    if(!is.null(endog_var)){
      ES_treatcontrol_means_fs <- ES_data[,list(rn="treatment_means", estimate = mean(get(endog_var)), cluster_se = sd(get(endog_var))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]
    } else{
      ES_treatcontrol_means_fs <- NULL
    }
  }

  # Collect CATT-specific-sample means of the outcome in the omitted_event_time (for scaling results)
  means_rf <- ES_data[ref_event_time == omitted_event_time, list(mean_outcome = mean(get(outcomevar))), by = list(ref_onset_time,catt_specific_sample)][order(ref_onset_time,catt_specific_sample)]
  means_rf[, pooled_mean := ES_data[ref_event_time == omitted_event_time, mean(get(outcomevar))]]
  treat_means_rf <- ES_data[ref_event_time == omitted_event_time & treated == 1, list(treat_mean_outcome = mean(get(outcomevar))), by = list(ref_onset_time,catt_specific_sample)][order(ref_onset_time,catt_specific_sample)]
  treat_means_rf[, treat_pooled_mean := ES_data[ref_event_time == omitted_event_time & treated == 1, mean(get(outcomevar))]]

  mapping <- ES_data[, .N, by = list(ref_onset_time,catt_specific_sample, ref_event_time)][order(ref_onset_time,catt_specific_sample,ref_event_time)]
  mapping <- mapping[ref_event_time != omitted_event_time]

  means_rf <- merge(means_rf, mapping, by = c("ref_onset_time", "catt_specific_sample"), all.x = TRUE)
  treat_means_rf <- merge(treat_means_rf, mapping, by = c("ref_onset_time", "catt_specific_sample"), all.x = TRUE)

  if(!is.null(endog_var)){
    means_fs <- ES_data[ref_event_time == omitted_event_time, list(mean_outcome = mean(get(endog_var))), by = list(ref_onset_time,catt_specific_sample)][order(ref_onset_time,catt_specific_sample)]
    means_fs[, pooled_mean := ES_data[ref_event_time == omitted_event_time, mean(get(endog_var))]]
    treat_means_fs <- ES_data[ref_event_time == omitted_event_time & treated == 1, list(treat_mean_outcome = mean(get(endog_var))), by = list(ref_onset_time,catt_specific_sample)][order(ref_onset_time,catt_specific_sample)]
    treat_means_fs[, treat_pooled_mean := ES_data[ref_event_time == omitted_event_time & treated == 1, mean(get(endog_var))]]
    means_fs <- merge(means_fs, mapping, by = c("ref_onset_time", "catt_specific_sample"), all.x = TRUE)
    treat_means_fs <- merge(treat_means_fs, mapping, by = c("ref_onset_time", "catt_specific_sample"), all.x = TRUE)
  } else{
    means_fs <- NULL
    treat_means_fs <- NULL
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

  figdata_rf <- rbindlist(list(ES_results_hetero_rf, ES_results_homo_rf, ES_treatcontrol_means_rf), use.names = TRUE, fill=TRUE)
  rm(ES_results_hetero_rf, ES_results_homo_rf, ES_treatcontrol_means_rf)
  figdata_rf[, model := "rf"]

  if(!is.null(endog_var)){
    figdata_fs <- rbindlist(list(ES_results_hetero_fs, ES_results_homo_fs, ES_treatcontrol_means_fs), use.names = TRUE, fill=TRUE)
    rm(ES_results_hetero_fs, ES_results_homo_fs, ES_treatcontrol_means_fs)
    figdata_fs[, model := "fs"]
  } else{
    figdata_fs <- NULL
  }

  figdata <- rbindlist(list(figdata_rf, figdata_fs), use.names = TRUE, fill = TRUE)
  rm(figdata_rf)
  rm(figdata_fs)

  # merge various *_unique_units into figdata
  # exclude catt_treated_unique_units/catt_total_unique_units as these will be added as part of (un)weighted means below

  # calculate unweighted and V1/V2 weighted means
  catt_treated_unique_units[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  catt_total_unique_units[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  figdata <- merge(figdata, catt_treated_unique_units, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
  figdata <- merge(figdata, catt_total_unique_units, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)

  subsets_for_avgs <- figdata[rn %in% c("catt")]
  subsets_for_avgs[, unweighted_estimate := mean(estimate, na.rm = TRUE), by = list(ref_event_time, model)]
  subsets_for_avgs[, cohort_weight_V1 := catt_treated_unique_units / sum(catt_treated_unique_units, na.rm = TRUE), by = list(ref_event_time, model)]
  subsets_for_avgs[, cohort_weight_V2 := catt_total_unique_units / sum(catt_total_unique_units, na.rm = TRUE), by = list(ref_event_time, model)]
  subsets_for_avgs[, weighted_estimate_V1 := weighted.mean(x = estimate, w = cohort_weight_V1, na.rm = TRUE), by = list(ref_event_time, model)]
  subsets_for_avgs[, weighted_estimate_V2 := weighted.mean(x = estimate, w = cohort_weight_V2, na.rm = TRUE), by = list(ref_event_time, model)]

  # merge weights into figdata
  weights <- subsets_for_avgs[, list(ref_event_time, ref_onset_time, model, cohort_weight_V1, cohort_weight_V2)]
  figdata <- merge(figdata, weights, by = c("ref_onset_time", "ref_event_time", "model"), all.x = TRUE, sort = FALSE)

  # rbind the (un)weighted avg estimates to figdata much like the "Pooled" data
  subsets_for_avgs[, rowid := seq_len(.N), by = list(ref_event_time, model)]
  subsets_for_avgs <- subsets_for_avgs[rowid == 1 | is.na(rowid)]

  unweighted <- subsets_for_avgs[, list(ref_event_time, unweighted_estimate, model)]
  unweighted[, ref_onset_time := "Equally-Weighted"]
  unweighted[, rn := "att"]
  setnames(unweighted, c("unweighted_estimate"), c("estimate"))
  setorderv(unweighted, c("ref_event_time", "model"))

  weighted_V1 <- subsets_for_avgs[, list(ref_event_time, weighted_estimate_V1, model)]
  weighted_V1[, ref_onset_time := "Cohort-Weighted"]
  weighted_V1[, rn := "att"]
  setnames(weighted_V1, c("weighted_estimate_V1"), c("estimate"))
  setorderv(weighted_V1, c("ref_event_time", "model"))

  weighted_V2 <- subsets_for_avgs[, list(ref_event_time, weighted_estimate_V2, model)]
  weighted_V2[, ref_onset_time := "Cohort-Weighted V2"]
  weighted_V2[, rn := "att"]
  setnames(weighted_V2, c("weighted_estimate_V2"), c("estimate"))
  setorderv(weighted_V2, c("ref_event_time", "model"))

  figdata <- rbindlist(list(figdata, unweighted, weighted_V1, weighted_V2), use.names = TRUE, fill=TRUE)
  figdata[is.na(cluster_se), cluster_se := 0]

  rm(subsets_for_avgs)
  rm(weights)
  rm(unweighted)
  rm(weighted_V1)
  rm(weighted_V2)
  rm(catt_treated_unique_units)
  rm(catt_total_unique_units)
  gc()

  if(homogeneous_ATT == FALSE){

    # calculate SEs for weighted avg estimates using catt_coefs and catt_vcov
    # the unique "Weighted" ref_event_time values represent the target list of parameters
    # for each such ref_event_time, want to extract the location (number) of the relevant parameters in catt_coefs
    # then will need to grab the relevant weight, and then construct the formula to supply to delta_method()

    event_times <- setdiff(figdata[, sort(unique(ref_event_time))], omitted_event_time)
    onset_times <- as.integer(figdata[rn == "catt", sort(unique(ref_onset_time))])

    min_onset_time <- min(onset_times)
    max_onset_time <- max(onset_times)

    if(!is.null(endog_var)){
      model_vals <- c("rf", "fs")
    } else{
      model_vals <- "rf"
    }

    for(mm in model_vals){

      for(et in event_times){

        if(et < 0){
          lookfor <- sprintf("cattlead%s$", abs(et))
          # crucial to have the end-of-line anchor "$" above; otherwise will find, e.g.,  -1 and -19:-10 event times
        } else{
          lookfor <- sprintf("catt%s$", abs(et))
          # crucial to have the end-of-line anchor "$" above; otherwise will find, e.g.,  1 and 10:19 event times
        }
        coef_indices <- grep(lookfor, names(get(sprintf("catt_coefs_%s", mm))))
        rm(lookfor)
        temp <- as.data.table(do.call(cbind, list(get(sprintf("catt_coefs_%s", mm))[coef_indices], coef_indices)), keep.rownames = TRUE)
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
        temp <- merge(temp, figdata[rn == "catt" & model == mm], by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
        temp <- temp[, list(ref_onset_time, ref_event_time, coef_index, cohort_weight_V1, cohort_weight_V2)]
        temp[, equal_weight := 1 / .N]

        temp[, equal_w_formula_entry := sprintf("(%s*x%s)", equal_weight, coef_index)]
        temp[, cohort_w_v1_formula_entry := sprintf("(%s*x%s)", cohort_weight_V1, coef_index)]
        temp[, cohort_w_v2_formula_entry := sprintf("(%s*x%s)", cohort_weight_V2, coef_index)]

        equal_w_g_formula_input = paste0(temp$equal_w_formula_entry, collapse = "+")
        cohort_w_v1_g_formula_input = paste0(temp$cohort_w_v1_formula_entry, collapse = "+")
        cohort_w_v2_g_formula_input = paste0(temp$cohort_w_v2_formula_entry, collapse = "+")

        figdata[rn == "att" & model == mm & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Equally-Weighted",
                cluster_se := delta_method(g = as.formula(paste("~", equal_w_g_formula_input)),
                                           mean = get(sprintf("catt_coefs_%s", mm)),
                                           cov = get(sprintf("catt_vcov_%s", mm)),
                                           ses = TRUE
                )
                ]

        figdata[rn == "att" & model == mm & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Cohort-Weighted",
                cluster_se := delta_method(g = as.formula(paste("~", cohort_w_v1_g_formula_input)),
                                           mean = get(sprintf("catt_coefs_%s", mm)),
                                           cov = get(sprintf("catt_vcov_%s", mm)),
                                           ses = TRUE
                )
                ]

        figdata[rn == "att" & model == mm & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Cohort-Weighted V2",
                cluster_se := delta_method(g = as.formula(paste("~", cohort_w_v2_g_formula_input)),
                                           mean = get(sprintf("catt_coefs_%s", mm)),
                                           cov = get(sprintf("catt_vcov_%s", mm)),
                                           ses = TRUE
                )
                ]

        rm(temp, equal_w_g_formula_input, cohort_w_v1_g_formula_input, cohort_w_v2_g_formula_input)

      }
    }
    rm(et, mm)
    gc()
  }

  # merge in the sample sizes
  figdata <- merge(figdata, event_time_total_unique_units, by = "ref_event_time", all.x = TRUE, sort = FALSE)
  rm(event_time_total_unique_units)

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
      dt[, unweighted_estimate := mean(estimate, na.rm = TRUE), by = list(ref_event_time, model)]
      dt[, weighted_estimate_V1 := weighted.mean(x = estimate, w = cohort_weight_V1, na.rm = TRUE), by = list(ref_event_time, model)]
      dt[, weighted_estimate_V2 := weighted.mean(x = estimate, w = cohort_weight_V2, na.rm = TRUE), by = list(ref_event_time, model)]
      dt[, rowid := seq_len(.N), by = list(ref_event_time, model)]
      result <- dt[rowid == 1 | is.na(rowid)]
      dt <- dt[, list(ref_event_time, ref_onset_time, model, cohort_weight_V1, cohort_weight_V2, grouping)]

      ew[[j]] <- result[, list(grouping, model, unweighted_estimate)]
      ew[[j]][, unweighted_estimate := mean(unweighted_estimate, na.rm = TRUE), by = list(grouping, model)]
      ew[[j]][, ref_onset_time := "Equally-Weighted + Collapsed"]
      ew[[j]][, rn := "att"]
      ew[[j]][, rowid := seq_len(.N), by = list(model)]
      setnames(ew[[j]], c("unweighted_estimate"), c("estimate"))
      ew[[j]] <- ew[[j]][rowid == 1 | is.na(rowid)]
      ew[[j]][, rowid := NULL]
      ew[[j]][, cluster_se := 0]

      cw[[j]] <- result[, list(grouping, model, weighted_estimate_V1)]
      cw[[j]][, weighted_estimate_V1 := mean(weighted_estimate_V1, na.rm = TRUE), by = list(grouping, model)]
      cw[[j]][, ref_onset_time := "Cohort-Weighted + Collapsed"]
      cw[[j]][, rn := "att"]
      cw[[j]][, rowid := seq_len(.N), by = list(model)]
      setnames(cw[[j]], c("weighted_estimate_V1"), c("estimate"))
      cw[[j]] <- cw[[j]][rowid == 1 | is.na(rowid)]
      cw[[j]][, rowid := NULL]
      cw[[j]][, cluster_se := 0]

      cw2[[j]] <- result[, list(grouping, model, weighted_estimate_V2)]
      cw2[[j]][, weighted_estimate_V2 := mean(weighted_estimate_V2, na.rm = TRUE), by = list(grouping, model)]
      cw2[[j]][, ref_onset_time := "Cohort-Weighted V2 + Collapsed"]
      cw2[[j]][, rn := "att"]
      cw2[[j]][, rowid := seq_len(.N), by = list(model)]
      setnames(cw2[[j]], c("weighted_estimate_V2"), c("estimate"))
      cw2[[j]] <- cw2[[j]][rowid == 1 | is.na(rowid)]
      cw2[[j]][, rowid := NULL]
      cw2[[j]][, cluster_se := 0]

      rm(result)

      if(!is.null(endog_var)){
        model_vals <- c("rf", "fs")
      } else{
        model_vals <- "rf"
      }

      for(mm in model_vals){

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
          coef_indices <- grep(lookfor, names(get(sprintf("catt_coefs_%s", mm))))
          rm(lookfor)
          temp <- as.data.table(do.call(cbind, list(get(sprintf("catt_coefs_%s", mm))[coef_indices], coef_indices)), keep.rownames = TRUE)
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
          temp <- merge(temp, dt[model == mm], by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
          temp[, weight_V0 := 1 / .N]

          templist[[i]] <- copy(temp)
          rm(temp)
          gc()

        }
        rm(i)

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

        ew[[j]][grouping == g & model == mm,
                cluster_se := delta_method(g = as.formula(paste("~", formula_input_ew)),
                                           mean = get(sprintf("catt_coefs_%s", mm)),
                                           cov = get(sprintf("catt_vcov_%s", mm)),
                                           ses = TRUE
                )
                ]

        cw[[j]][grouping == g & model == mm,
                cluster_se := delta_method(g = as.formula(paste("~", formula_input_cw)),
                                           mean = get(sprintf("catt_coefs_%s", mm)),
                                           cov = get(sprintf("catt_vcov_%s", mm)),
                                           ses = TRUE
                )
                ]

        cw2[[j]][grouping == g & model == mm,
                 cluster_se := delta_method(g = as.formula(paste("~", formula_input_cw2)),
                                            mean = get(sprintf("catt_coefs_%s", mm)),
                                            cov = get(sprintf("catt_vcov_%s", mm)),
                                            ses = TRUE
                 )
                 ]

        rm(formula_input_ew)
        rm(formula_input_cw)
        rm(formula_input_cw2)

      }
      rm(dt)

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
    rm(ew, cw, cw2)
    rm(collapsed_estimate_total_unique_units)
    gc()
  }

  # add overall unique units count
  figdata[, total_unique_units := unique_units]
  rm(unique_units)

  # 'grouping' only relevant if calculate_collapse_estimates == TRUE, but easier to just generate a dummy column if not
  if(calculate_collapse_estimates == FALSE){
    figdata[, grouping := NA]
    figdata[, collapsed_estimate_total_unique_units := NA]
  }

  # form ratios across models (rf/fs) by (ref_onset_time, ref_event_time, grouping)
  if(!is.null(endog_var)){

    ratio_data_rf <- figdata[rn %in% c("att", "catt") & model == "rf", list(ref_onset_time, ref_event_time, grouping, estimate, rn)]
    setnames(ratio_data_rf, "estimate", "estimate_rf")

    ratio_data_fs <- figdata[rn %in% c("att", "catt") & model == "fs", list(ref_onset_time, ref_event_time, grouping, estimate, rn)]
    setnames(ratio_data_fs, "estimate", "estimate_fs")

    ratio_data <- merge(ratio_data_rf, ratio_data_fs, by = c("ref_onset_time", "ref_event_time", "grouping", "rn"), sort = FALSE)
    rm(ratio_data_rf, ratio_data_fs)
    ratio_data[, estimate := (estimate_rf / estimate_fs)]
    ratio_data[, c("estimate_rf", "estimate_fs") := NULL]

    # merge in various counts to ratio data, then rbindlist back with figdata
    # counts and weights will be the same across 'model' by (ref_onset_time, ref_event_time, grouping), so just need one
    counts_and_weights <- figdata[rn %in% c("att", "catt") & model == "rf", list(ref_onset_time, ref_event_time, grouping, rn, reg_sample_size, catt_treated_unique_units, catt_total_unique_units, cohort_weight_V1, cohort_weight_V2, event_time_total_unique_units, collapsed_estimate_total_unique_units, total_unique_units)]
    ratio_data <- merge(ratio_data, counts_and_weights, by = c("ref_onset_time", "ref_event_time", "grouping", "rn"), sort = FALSE)
    ratio_data[, model := "ratio"]
    ratio_data[, cluster_se := 0]

    figdata <- rbindlist(list(figdata, ratio_data), use.names = TRUE, fill = TRUE)
    rm(ratio_data)
  }

  # addin model name / renaming to make the models crystal clear; model == "ratio" already pretty clear
  figdata[model == "rf", model := "reduced_form"]
  figdata[model == "fs", model := "first_stage"]

  means_rf[, model := "reduced_form"]
  treat_means_rf[, model := "reduced_form"]

  if(!is.null(endog_var)){
    means_fs[, model := "first_stage"]
    treat_means_fs[, model := "first_stage"]
  }

  ###
  # TO-DO: two-sample IV / delta-method SES
  ###

  # start bootstrap run if relevant
  if(bootstrapES == TRUE){

    # all arguments of Wald_bootstrap_ES are passed but for 'iter', which will be supplied by 1:bootstrap_iters below
    boot_results <- rbindlist(parallel::mclapply(X = 1:bootstrap_iters,
                                                 FUN = Wald_bootstrap_ES,
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
                                                 ntile_var = ntile_var, ntile_event_time = ntile_event_time, ntiles = ntiles, ntile_var_value = ntile_var_value, ntile_avg = ntile_avg,
                                                 endog_var = endog_var, cohort_by_cohort = cohort_by_cohort, cohort_by_cohort_num_cores = cohort_by_cohort_num_cores),
                              use.names = TRUE
    )

    if(keep_all_bootstrap_results == TRUE){
      bootstrap_output <- copy(boot_results)
    }

    if(calculate_collapse_estimates == TRUE & homogeneous_ATT == FALSE){
      boot_results_collapse_estimates <- boot_results[!is.na(grouping)]
      boot_results <- boot_results[is.na(grouping)]
      boot_results[, grouping := NULL]
      gc()
    }

    # Construct bootstrap CI for all of the estimates
    boot_results_ses <- boot_results[, list(bootstrap_se = sd(estimate)), by = list(ref_onset_time, ref_event_time, model)]
    order_to_restore <- na.omit(unique(c(copy(colnames(figdata)),copy(colnames(boot_results_ses)))))
    figdata <- merge(figdata, boot_results_ses, by = c("ref_onset_time", "ref_event_time", "model"), all.x = TRUE, sort = FALSE)
    setcolorder(figdata, neworder = order_to_restore)
    rm(order_to_restore)
    figdata[rn == "treatment_means", bootstrap_se := NA]

    if(calculate_collapse_estimates == TRUE & homogeneous_ATT == FALSE){
      boot_results_collapse_estimates_ses <- boot_results_collapse_estimates[, list(bootstrap_se_collapse = sd(estimate)), by = list(ref_onset_time, grouping, model)]
      order_to_restore <- na.omit(unique(c(copy(colnames(figdata)),copy(colnames(boot_results_collapse_estimates_ses)))))
      figdata <- merge(figdata, boot_results_collapse_estimates_ses, by = c("ref_onset_time", "grouping", "model"), all.x = TRUE, sort = FALSE)
      setcolorder(figdata, neworder = order_to_restore)
      rm(order_to_restore)
      figdata[!is.na(bootstrap_se_collapse), bootstrap_se := bootstrap_se_collapse]
      figdata[, bootstrap_se_collapse := NULL]
    }

    original_sample <- NULL
    gc()

  }

  if(ipw_keep_data == TRUE){
    figdata <- rbindlist(list(figdata, ipw_dt), use.names = TRUE, fill = TRUE)
  }

  # Wrapping up all results into a list

  return_list = list()
  return_list[[1]] <- figdata
  return_list[[2]] <- list(means_rf, means_fs)
  return_list[[3]] <- list(treat_means_rf, treat_means_fs)

  if(homogeneous_ATT == FALSE){
    return_list[[4]] <- list(catt_coefs_rf, catt_vcov_rf)
    if(!is.null(endog_var)){
      return_list[[5]] <- list(catt_coefs_fs, catt_vcov_fs)
    }
  }

  if(ipw_keep_data == TRUE){
    return_list[[6]] <- ipw_dt
  }

  if(keep_all_bootstrap_results == TRUE){
    return_list[[7]] <- bootstrap_output
  }

  flog.info('Wald_ES is finished.')

  return(return_list)

}

#' @export
# Need 'iter' to come first given how R reads the lapply/mclapply implementing Wald_bootstrap_ES()
Wald_bootstrap_ES <- function(iter, long_data, outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars,
                              omitted_event_time= -2, anticipation = 0, min_control_gap=1, max_control_gap=Inf,
                              linearize_pretrends=FALSE, fill_zeros=FALSE, residualize_covariates = FALSE,
                              control_subset_var=NA, control_subset_event_time=0,
                              treated_subset_var=NA, treated_subset_event_time=0,
                              control_subset_var2=NA, control_subset_event_time2=0,
                              treated_subset_var2=NA, treated_subset_event_time2=0,
                              discrete_covars = NULL, cont_covars = NULL, never_treat_action = 'none',
                              homogeneous_ATT = FALSE, reg_weights = NULL, add_unit_fes = FALSE,
                              ipw = FALSE, ipw_model = 'linear', ipw_composition_change = FALSE, ipw_keep_data = FALSE, ipw_ps_lower_bound = 0, ipw_ps_upper_bound = 1,
                              event_vs_noevent = FALSE,
                              ref_discrete_covars = NULL, ref_discrete_covar_event_time=0, ref_cont_covars = NULL, ref_cont_covar_event_time=0,
                              calculate_collapse_estimates = FALSE, collapse_inputs = NULL,
                              ref_reg_weights = NULL, ref_reg_weights_event_time=0,
                              ntile_var = NULL, ntile_event_time = -2, ntiles = NA, ntile_var_value = NA, ntile_avg = FALSE,
                              endog_var = NULL, cohort_by_cohort = FALSE, cohort_by_cohort_num_cores = 1){

  # Function mirrors Wald_ES(), but with less printing/warning and skipping sections responsible for analytical standard errors and the like which aren't relevant
  # We keep various checks to make sure no issues arise due to a potential strange bootstrap sample draw

  bs_long_data <- block_sample(long_data = copy(long_data), unit_var = unit_var, cal_time_var = cal_time_var)
  gc()

  flog.info(sprintf("Beginning Wald_bootstrap_ES iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))

  # type checks
  assertDataTable(bs_long_data)
  assertCharacter(outcomevar,len=1)
  assertCharacter(unit_var,len=1)
  assertCharacter(cal_time_var,len=1)
  if(!is.null(cluster_vars)){assertCharacter(cluster_vars)}
  assertIntegerish(omitted_event_time,len=1,upper=-1)
  assertIntegerish(anticipation,len=1,lower=0)
  assertIntegerish(min_control_gap,len=1,lower=1)
  if(!any(testIntegerish(max_control_gap,len=1,lower=min_control_gap),is.infinite(max_control_gap))){
    assertIntegerish(max_control_gap,len=1,lower=min_control_gap)
  }
  if(!any(testIntegerish(max_control_gap,len=1,lower=min_control_gap),is.infinite(max_control_gap))){
    assertIntegerish(max_control_gap,len=1,lower=min_control_gap)
  }
  if(!is.na(control_subset_var)){
    assertCharacter(control_subset_var,len=1)
  }
  if(!is.na(treated_subset_var)){
    assertCharacter(treated_subset_var,len=1)
  }
  if(!is.na(control_subset_var2)){
    assertCharacter(control_subset_var2,len=1)
  }
  if(!is.na(treated_subset_var2)){
    assertCharacter(treated_subset_var2,len=1)
  }
  if(!is.null(reg_weights)){
    assertCharacter(reg_weights,len=1)
  }
  if(!is.null(ref_reg_weights)){
    assertCharacter(ref_reg_weights,len=1)
  }
  assertIntegerish(control_subset_event_time,len=1)
  assertIntegerish(treated_subset_event_time,len=1)
  assertIntegerish(control_subset_event_time2,len=1)
  assertIntegerish(treated_subset_event_time2,len=1)
  assertIntegerish(ref_discrete_covar_event_time,len=1)
  assertIntegerish(ref_cont_covar_event_time,len=1)
  assertIntegerish(ref_reg_weights_event_time,len=1)
  assertFlag(linearize_pretrends)
  assertFlag(fill_zeros)
  assertFlag(residualize_covariates)
  if(residualize_covariates){
    if(!any(testCharacter(discrete_covars), testCharacter(cont_covars))){
      stop("Since residualize_covariates=TRUE, either discrete_covars or cont_covars must be provided as a character vector.")
    }
  }
  assertFlag(homogeneous_ATT)
  assertFlag(add_unit_fes)
  assertFlag(ipw)
  assertFlag(ipw_composition_change)
  assertNumber(ipw_ps_lower_bound,lower=0,upper=1)
  assertNumber(ipw_ps_upper_bound,lower=0,upper=1)
  assertFlag(calculate_collapse_estimates)
  if(calculate_collapse_estimates == TRUE){
    assertDataTable(collapse_inputs, ncols = 2, any.missing = FALSE, types = c("character", "list"))
  }
  if(!is.null(ntile_var)){
    assertIntegerish(ntile_event_time,len=1)
    assertIntegerish(ntiles,len=1,lower=2)
    assertIntegerish(ntile_var_value,len=1,lower=1,upper=ntiles)
    assertFlag(ntile_avg)
  }
  if(!is.null(endog_var)){
    assertCharacter(endog_var,len=1)
  }
  assertFlag(cohort_by_cohort)
  assertIntegerish(cohort_by_cohort_num_cores,len=1,lower=1)

  # check that anticipation choice and omitted_event_time choice don't conflict
  if(omitted_event_time + anticipation > -1){
    stop(sprintf("omitted_event_time='%s' and anticipation='%s' implies overlap of pre-treatment and anticipation periods. Let me suggest omitted_event_time<='%s'",
                 omitted_event_time, anticipation, ((-1 * anticipation) - 1)))
  }

  # check that all of these variables are actually in the data.table, and provide custom error messages.
  if(!(outcomevar %in% names(bs_long_data))){stop(sprintf("Variable outcomevar='%s' is not in the bs_long_data you provided.",outcomevar))}
  if(!(unit_var %in% names(bs_long_data))){stop(sprintf("Variable unit_var='%s' is not in the bs_long_data you provided.",unit_var))}
  if(!(cal_time_var %in% names(bs_long_data))){stop(sprintf("Variable cal_time_var='%s' is not in the bs_long_data you provided.",cal_time_var))}
  if(!(onset_time_var %in% names(bs_long_data))){stop(sprintf("Variable onset_time_var='%s' is not in the bs_long_data you provided.",onset_time_var))}
  if(!is.null(cluster_vars)){
    for(vv in cluster_vars){
      if(!(vv %in% names(bs_long_data))){stop(sprintf("Variable cluster_vars='%s' is not in the bs_long_data you provided. Let me suggest cluster_vars='%s'.",vv,unit_var))}
    }
  }
  if(!is.na(control_subset_var)){
    if(!(control_subset_var %in% names(bs_long_data))){stop(sprintf("Variable control_subset_var='%s' is not in the bs_long_data you provided.",control_subset_var))}
    if(!(bs_long_data[,typeof(get(control_subset_var))]=="logical")){stop(sprintf("Variable control_subset_var='%s' must be of type logical (i.e., only TRUE or FALSE values).",control_subset_var))}
  }
  if(!is.na(treated_subset_var)){
    if(!(treated_subset_var %in% names(bs_long_data))){stop(sprintf("Variable treated_subset_var='%s' is not in the bs_long_data you provided.",treated_subset_var))}
    if(!(bs_long_data[,typeof(get(treated_subset_var))]=="logical")){stop(sprintf("Variable treated_subset_var='%s' must be of type logical (i.e., only TRUE or FALSE values).",treated_subset_var))}
  }
  if(!is.na(control_subset_var2)){
    if(!(control_subset_var2 %in% names(bs_long_data))){stop(sprintf("Variable control_subset_var2='%s' is not in the bs_long_data you provided.",control_subset_var2))}
    if(!(bs_long_data[,typeof(get(control_subset_var2))]=="logical")){stop(sprintf("Variable control_subset_var2='%s' must be of type logical (i.e., only TRUE or FALSE values).",control_subset_var2))}
  }
  if(!is.na(treated_subset_var2)){
    if(!(treated_subset_var2 %in% names(bs_long_data))){stop(sprintf("Variable treated_subset_var2='%s' is not in the bs_long_data you provided.",treated_subset_var2))}
    if(!(bs_long_data[,typeof(get(treated_subset_var2))]=="logical")){stop(sprintf("Variable treated_subset_var2='%s' must be of type logical (i.e., only TRUE or FALSE values).",treated_subset_var2))}
  }
  if(testCharacter(discrete_covars)){
    for(vv in discrete_covars){
      if(!(vv %in% names(bs_long_data))){stop(sprintf("Variable discrete_covars='%s' is not in the bs_long_data you provided.",vv))}
    }
  }
  if(testCharacter(cont_covars)){
    for(vv in cont_covars){
      if(!(vv %in% names(bs_long_data))){stop(sprintf("Variable cont_covars='%s' is not in the bs_long_data you provided.",vv))}
    }
  }
  if(testCharacter(ref_discrete_covars)){
    for(vv in ref_discrete_covars){
      if(!(vv %in% names(bs_long_data))){stop(sprintf("Variable ref_discrete_covars='%s' is not in the bs_long_data you provided.",vv))}
    }
  }
  if(testCharacter(ref_cont_covars)){
    for(vv in ref_cont_covars){
      if(!(vv %in% names(bs_long_data))){stop(sprintf("Variable ref_cont_covars='%s' is not in the bs_long_data you provided.",vv))}
    }
  }
  if(testCharacter(reg_weights)){
    if(!(reg_weights %in% names(bs_long_data))){stop(sprintf("Variable reg_weights='%s' is not in the bs_long_data you provided.",reg_weights))}
  }
  if(testCharacter(ref_reg_weights)){
    if(!(ref_reg_weights %in% names(bs_long_data))){stop(sprintf("Variable ref_reg_weights='%s' is not in the bs_long_data you provided.",ref_reg_weights))}
  }
  if(testCharacter(ntile_var)){
    if(!(ntile_var %in% names(bs_long_data))){stop(sprintf("Variable ntile_var='%s' is not in the bs_long_data you provided.",ntile_var))}
  }
  if(testCharacter(endog_var)){
    if(!(endog_var %in% names(bs_long_data))){stop(sprintf("Variable endog_var='%s' is not in the bs_long_data you provided.",endog_var))}
  }

  # temporary check -- code is currently not set up to accept ntile_event_time != omitted_event_time
  # only bites in the case where !is.null(ntile_var) and omitted_event_time != ntile_event_time
  if((!is.null(ntile_var)) & omitted_event_time != ntile_event_time){
    stop(sprintf("ntile_event_time='%s', but currently code can only accept ntile_event_time == omitted_event_time (e.g., %s).", ntile_event_time, omitted_event_time))
  }

  # check that calculate_collapse_estimates == TRUE only if homogeneous_ATT == FALSE
  if(calculate_collapse_estimates == TRUE & homogeneous_ATT == TRUE){
    stop("Cannot have calculate_collapse_estimates == TRUE & homogeneous_ATT == TRUE. Consider setting homogeneous_ATT = FALSE.")
  }

  # check that control variables don't overlap with design variables (e.g., cal_time_var, and onset_time_var or unit_var)
  if(add_unit_fes == TRUE){
    design_vars <- c(cal_time_var, unit_var)
  } else{
    design_vars <- c(cal_time_var, onset_time_var)
  }
  if(testCharacter(discrete_covars)){
    for(vv in discrete_covars){
      if(vv %in% design_vars){stop(sprintf("Variable discrete_covars='%s' is among c('%s','%s') which are already controlled in the design.",vv, design_vars[[1]], design_vars[[2]]))}
    }
  }
  if(testCharacter(cont_covars)){
    for(vv in cont_covars){
      if(vv %in% design_vars){stop(sprintf("Variable cont_covars='%s' is among c('%s','%s') which are already controlled in the design.",vv, design_vars[[1]], design_vars[[2]]))}
    }
  }
  if(testCharacter(ref_discrete_covars)){
    for(vv in ref_discrete_covars){
      if(vv %in% design_vars){stop(sprintf("Variable ref_discrete_covars='%s' is among c('%s','%s') which are already controlled in the design.",vv, design_vars[[1]], design_vars[[2]]))}
    }
  }
  if(testCharacter(ref_cont_covars)){
    for(vv in ref_cont_covars){
      if(vv %in% design_vars){stop(sprintf("Variable ref_cont_covars='%s' is among c('%s','%s') which are already controlled in the design.",vv, design_vars[[1]], design_vars[[2]]))}
    }
  }

  # check that user only supplied one of reg_weights / ref_reg_weights, but not both
  if(!is.null(reg_weights) & !is.null(ref_reg_weights)){
    stop("Supplied variables for both reg_weights and ref_reg_weights, but Wald_bootstrap_ES() only admits using one or the other type of weight (or neither).")
  }

  # check that user correctly input what to do with never treated
  if(!(never_treat_action %in% c('none', 'exclude', 'keep', 'only'))){
    stop(sprintf("never_treat_action='%s' is not among allowed values (c('none', 'exclude', 'keep', 'only')).", never_treat_action))
  }
  if(never_treat_action=='none' & dim(bs_long_data[is.na(get(onset_time_var))])[1] > 0){
    stop(sprintf("never_treat_action='%s' but some units have %s=NA. Please edit supplied bs_long_data or consider another option for never_treat_action.", never_treat_action, onset_time_var))
  }
  if(never_treat_action!='none' & dim(bs_long_data[is.na(get(onset_time_var))])[1] == 0){
    stop(sprintf("never_treat_action='%s' but no units have %s=NA. Let me suggest never_treat_action='none'.", never_treat_action, onset_time_var))
  }

  # check that ipw model conforms to available options, and that propensity score cutoffs are sensible
  if(ipw == TRUE){
    if(!(ipw_model %in% c('linear', 'logit', 'probit'))){
      stop(sprintf("ipw_model='%s' is not among allowed values (c('linear', 'logit', 'probit')).", ipw_model))
    }
  }
  if(ipw_ps_lower_bound > ipw_ps_upper_bound){
    stop(sprintf("ipw_ps_lower_bound='%s' & ipw_ps_upper_bound='%s', which means ipw_ps_lower_bound > ipw_ps_upper_bound. Consider revising these cutoffs.", ipw_ps_lower_bound, ipw_ps_upper_bound))
  }

  # edit bs_long_data in line with supplied never_treat_action option
  if(never_treat_action == 'exclude'){
    never_treat_val = NA
    bs_long_data <- bs_long_data[!is.na(get(onset_time_var))]
    gc()
  } else if(never_treat_action %in% c('keep', 'only')){
    # assign the never-treated a unique onset time that ensures they are always part of the control group
    never_treat_val <- max(max(bs_long_data[[onset_time_var]], na.rm = TRUE),
                           max(bs_long_data[[cal_time_var]], na.rm = TRUE)
    ) + min_control_gap + anticipation + 1
    bs_long_data[is.na(get(onset_time_var)), (onset_time_var) := never_treat_val]
  }

  # fill with zeros
  if(fill_zeros){
    flog.info(sprintf("Filling in zeros in Wald_bootstrap_ES iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
    bs_long_data <- ES_expand_to_balance(long_data = bs_long_data,
                                         vars_to_fill = na.omit(unique(c(outcomevar, endog_var))),
                                         unit_var = unit_var,
                                         cal_time_var = cal_time_var,
                                         onset_time_var = onset_time_var)
  }

  # Check that there exist cohorts with observations at omitted_event_time
  if(is.infinite(suppressWarnings(bs_long_data[get(cal_time_var) - get(onset_time_var) == omitted_event_time, min(get(onset_time_var))]))){
    stop(sprintf("Variable onset_time_var='%s' has no treated groups with observations at pre-treatment event time %s.",onset_time_var, omitted_event_time))
  }

  # remove any cases with reg_weights <= 0, -Inf, Inf, or NA (if relevant)
  if(!(is.null(reg_weights))){
    bs_long_data <- bs_long_data[(!is.na(get(reg_weights))) & (!is.infinite(get(reg_weights))) & (get(reg_weights) > 0)]
    gc()
  }

  # linearize pre-trends; never-treated will be treated as a single cohort
  if(linearize_pretrends){
    # Not straightforward what to do in the case where !is.null(endog_var)
    flog.info("NOT CURRENTLY Linearizing pre-trends.")
  }

  # extrapolate contribution of covariates using pre-treated period
  if(residualize_covariates){
    # Not straightforward what to do in the case where !is.null(endog_var)
    flog.info("NOT CURRENTLY Residualizing on covariates.")
  }

  # process data
  flog.info(sprintf("Beginning data stacking in Wald_bootstrap_ES iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
  bs_ES_data <- ES_clean_data(long_data = bs_long_data, outcomevar = outcomevar,
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
                              endog_var = endog_var)

  # construct discrete covariates specific to a ref_event_time
  # will be time invariant for a given ref_onset_time == ref_discrete_covar_event_time, but time-varying across ref_onset_times

  if(!is.null(ref_discrete_covars)){

    start_cols <- copy(colnames(bs_ES_data))

    for(var in ref_discrete_covars){

      # If, for some reason, there is overlap (in variable name) between discrete_covars and ref_discrete_covars, want to make sure not to overwrite
      if(var %in% discrete_covars){

        varname <- sprintf("%s_dyn", var)

        # For a variable that is a character, will need to generate a numeric
        # Note: if there are missings, this will assign them all to a single level
        if(class(bs_ES_data[[var]]) == "character"){
          bs_ES_data[, (varname) := .GRP, by = var]
        } else{
          bs_ES_data[, (varname) := get(var)]
        }

        # Define the time-invariant outcome for all units within a ref_onset_time
        # use sum() instead of max() in case 'varname' takes on negative values (as max() would then return 0)
        bs_ES_data[, (varname) := sum(as.integer(get(varname)*(ref_event_time==ref_discrete_covar_event_time))), by=c(unit_var, "ref_onset_time")]

      } else{

        varname <- var

        # For a variable that is a character, will need to generate a numeric
        if(class(bs_ES_data[[var]]) == "character"){
          bs_ES_data[, varnum := .GRP, by = var]
          bs_ES_data[, (varname) := NULL]
          setnames(bs_ES_data, "varnum", varname)
          # ^need to do it this way because (varname) is character and redefining and assigning in one step is a pain
        }

        # Define the time-invariant outcome for all units within a ref_onset_time
        # use sum() instead of max() in case 'varname' takes on negative values (as max() would then return 0)
        bs_ES_data[, (varname) := sum(as.integer(get(varname)*(ref_event_time==ref_discrete_covar_event_time))), by=c(unit_var, "ref_onset_time")]

      }

    }

    # If any columns were generated above, we can add them to ref_discrete_covars
    # When these are used later in estimation, the loop is over unique values in the intersection of discrete_covars and ref_discrete_covars
    ref_discrete_covars <- unique(na.omit(c(ref_discrete_covars, setdiff(colnames(bs_ES_data), start_cols))))
    rm(start_cols)
  }

  # construct continuous covariates specific to a ref_event_time
  # will be time invariant for a given ref_onset_time, but time-varying across ref_onset_times

  if(!is.null(ref_cont_covars)){

    # Will construct this control variable which will be time-invariant, but defined as of the ref_cont_covar_event_time relative to each ref_onset_time
    # Shouldn't have any missings, but in the code, will just treat this as another category

    start_cols <- copy(colnames(bs_ES_data))

    for(var in ref_cont_covars){

      # If, for some reason, there is overlap (in variable name) between cont_covars and ref_cont_covars, want to make sure not to overwrite
      if(var %in% cont_covars){

        varname <- sprintf("%s_dyn", var)

        # Define the time-invariant outcome for all units within a ref_onset_time
        if(class(bs_ES_data[[var]]) == "integer"){
          # needed to include type checks as sum() often converts from integer to numeric
          # use sum() instead of max() in case 'var' takes on negative values (as max() would then return 0)
          bs_ES_data[, (varname) := sum(as.integer(get(var)*(ref_event_time==ref_cont_covar_event_time))), by=c(unit_var, "ref_onset_time")]
        } else{
          bs_ES_data[, (varname) := sum(get(var)*(ref_event_time==ref_cont_covar_event_time)), by=c(unit_var, "ref_onset_time")]
        }
      } else{

        # Define the time-invariant outcome for all units within a ref_onset_time
        if(class(bs_ES_data[[var]]) == "integer"){
          # needed to include type checks as sum() often converts from integer to numeric
          # use sum() instead of max() in case 'var' takes on negative values (as max() would then return 0)
          bs_ES_data[, (var) := sum(as.integer(get(var)*(ref_event_time==ref_cont_covar_event_time))), by=c(unit_var, "ref_onset_time")]
        } else{
          bs_ES_data[, (var) := sum(get(var)*(ref_event_time==ref_cont_covar_event_time)), by=c(unit_var, "ref_onset_time")]
        }

      }

    }

    # If any columns were generated above, we can add them to ref_cont_covars
    # When these are used later in estimation, the loop is over unique values in the intersection of cont_covars and ref_cont_covars
    ref_cont_covars <- unique(na.omit(c(ref_cont_covars, setdiff(colnames(bs_ES_data), start_cols))))
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
    if(class(bs_ES_data[[ref_reg_weights]]) == "integer"){
      # needed to include type checks as sum() often converts from integer to numeric
      # use sum() instead of max() just to parallel earlier code; wouldn't anticipate negative values in a weight variable
      bs_ES_data[, (ref_reg_weights) := sum(as.integer(get(ref_reg_weights)*(ref_event_time==ref_reg_weights_event_time))), by=c(unit_var, "ref_onset_time")]
    } else{
      bs_ES_data[, (ref_reg_weights) := sum(get(ref_reg_weights)*(ref_event_time==ref_reg_weights_event_time)), by=c(unit_var, "ref_onset_time")]
    }

    # remove any cases with ref_reg_weights <= 0, -Inf, Inf, or NA
    bs_ES_data <- bs_ES_data[(!is.na(get(ref_reg_weights))) & (!is.infinite(get(ref_reg_weights))) & (get(ref_reg_weights) > 0)]
    gc()
  }

  # estimate inverse probability weights, if relevant
  if(ipw == TRUE){

    # Within each DiD sample, estimate weights to balance provided covariates
    bs_ES_data[, did_id := .GRP, by = list(ref_onset_time, catt_specific_sample)]
    did_ids <- bs_ES_data[, sort(unique(did_id))]

    if(ipw_composition_change == FALSE){

      bs_ES_data[, pr_temp := as.numeric(NA)]

      for(i in did_ids){

        ipw_dt <- ES_make_ipw_dt(did_dt = copy(bs_ES_data[did_id == i]),
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

        bs_ES_data <- merge(bs_ES_data, ipw_dt, by = c(unit_var, cal_time_var, "did_id"), all.x = TRUE, sort = FALSE)
        bs_ES_data[is.na(pr_temp) & !is.na(pr), pr_temp := pr]
        bs_ES_data[, pr := NULL]
        ipw_dt <- NULL
        gc()
      }

      setnames(bs_ES_data, "pr_temp", "pr")

      # Restrict to propensity scores strictly in [a,b], a>0 and b<1 (default: a=0, b=1)
      bs_ES_data <- bs_ES_data[!is.na(pr) & !is.infinite(pr) & (between(pr, ipw_ps_lower_bound, ipw_ps_upper_bound, incbounds = FALSE))]
      gc()

      bs_ES_data[treated == 1, pw := 1]
      bs_ES_data[treated == 0, pw := (pr / (1 - pr))]

      if(ipw_keep_data == TRUE){
        ipw_dt <- bs_ES_data[, list(get(unit_var), ref_onset_time, ref_event_time, catt_specific_sample, treated, pr)]
        setnames(ipw_dt, "V1", unit_var)
      }
    }
  }

  # collect ATT estimates
  if(homogeneous_ATT == FALSE){

    # reduced form; homogeneous_ATT = FALSE
    # NOTE: need to copy(ES_data) as by_cohort_ES_estimate_ATT changes the underlying data
    ES_results_hetero_rf <- by_cohort_ES_estimate_ATT(ES_data = copy(bs_ES_data),
                                                      outcomevar = outcomevar,
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
                                                      cohort_by_cohort = cohort_by_cohort,
                                                      cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)[[1]]
    gc()
    setnames(ES_results_hetero_rf,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

    if(!is.null(endog_var)){
      # first_stage; homogeneous_ATT = FALSE
      # NOTE: need to copy(ES_data) as by_cohort_ES_estimate_ATT changes the underlying data
      ES_results_hetero_fs <- by_cohort_ES_estimate_ATT(ES_data = copy(bs_ES_data),
                                                        outcomevar = endog_var,
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
                                                        cohort_by_cohort = cohort_by_cohort,
                                                        cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)[[1]]
      gc()
      setnames(ES_results_hetero_fs,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))
    } else{
      ES_results_hetero_fs <- NULL
    }

  } else{
    ES_results_hetero_rf <- NULL
    ES_results_hetero_fs <- NULL
  }

  # reduced form; homogeneous_ATT = TRUE
  # NOTE: need to copy(ES_data) as by_cohort_ES_estimate_ATT changes the underlying data
  ES_results_homo_rf <- by_cohort_ES_estimate_ATT(ES_data = copy(bs_ES_data),
                                                  outcomevar = outcomevar,
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
                                                  cohort_by_cohort = cohort_by_cohort,
                                                  cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)[[1]]
  gc()
  setnames(ES_results_homo_rf,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

  if(!is.null(endog_var)){
    # first_stage; homogeneous_ATT = TRUE
    # NOTE: need to copy(ES_data) as by_cohort_ES_estimate_ATT changes the underlying data
    ES_results_homo_fs <- by_cohort_ES_estimate_ATT(ES_data = copy(bs_ES_data),
                                                    outcomevar = endog_var,
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
                                                    cohort_by_cohort = cohort_by_cohort,
                                                    cohort_by_cohort_num_cores = cohort_by_cohort_num_cores)[[1]]
    gc()
    setnames(ES_results_homo_fs,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))
  } else{
    ES_results_homo_fs <- NULL
  }

  # collect count of treated units by each (ref_onset_time, ref_event_time) for V1 of population-weighted ATTs
  # as we won't have an estimate for the omitted_event_time, exclude it below
  if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
    catt_treated_unique_units <- bs_ES_data[treated == 1 & (ref_event_time != omitted_event_time),list(catt_treated_unique_units = sum(get(na.omit(unique(c(reg_weights, ref_reg_weights)))))), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  } else{
    catt_treated_unique_units <- bs_ES_data[treated == 1 & (ref_event_time != omitted_event_time),list(catt_treated_unique_units = .N), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  }

  # collect count of treated+control units by each (ref_onset_time, ref_event_time) for V2 of population-weighted ATTs
  # as we won't have an estimate for the omitted_event_time, exclude it below
  if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
    catt_total_unique_units <- bs_ES_data[ref_event_time != omitted_event_time,list(catt_total_unique_units = sum(get(na.omit(unique(c(reg_weights, ref_reg_weights)))))), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  } else{
    catt_total_unique_units <- bs_ES_data[ref_event_time != omitted_event_time,list(catt_total_unique_units = .N), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  }

  # collect count of unique units that will be (implicitly) used to estimate cohort- and equally-weighted avgs for each event time
  # as we won't have an estimate for the omitted_event_time, exclude it below
  # will merge these onto figdata at the very end
  event_time_total_unique_units <- bs_ES_data[ref_event_time != omitted_event_time, list(event_time_total_unique_units = uniqueN(get(unit_var), na.rm = TRUE)), by = list(ref_event_time)][order(ref_event_time)]

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
      count <- bs_ES_data[(ref_event_time != omitted_event_time) & (ref_event_time %in% unique(na.omit(unlist(collapse_input_dt[name == g][[2]])))),
                          uniqueN(get(unit_var), na.rm = TRUE)
                          ]
      collapsed_estimate_total_unique_units[[i]] <- data.table(grouping = g, collapsed_estimate_total_unique_units = count)
    }
    rm(i)

    collapsed_estimate_total_unique_units <- rbindlist(collapsed_estimate_total_unique_units, use.names = TRUE)

  }

  # finally, outermost count of unique units in full bs_ES_data
  unique_units <- bs_ES_data[, uniqueN(get(unit_var), na.rm = TRUE)]
  gc()

  bs_ES_data <- NULL
  gc()

  figdata_rf <- rbindlist(list(ES_results_hetero_rf, ES_results_homo_rf), use.names = TRUE, fill=TRUE)
  rm(ES_results_hetero_rf, ES_results_homo_rf)
  figdata_rf[, model := "rf"]

  if(!is.null(endog_var)){
    figdata_fs <- rbindlist(list(ES_results_hetero_fs, ES_results_homo_fs), use.names = TRUE, fill=TRUE)
    rm(ES_results_hetero_fs, ES_results_homo_fs)
    figdata_fs[, model := "fs"]
  } else{
    figdata_fs <- NULL
  }

  figdata <- rbindlist(list(figdata_rf, figdata_fs), use.names = TRUE, fill = TRUE)
  rm(figdata_rf)
  rm(figdata_fs)

  # merge various *_unique_units into figdata
  # exclude catt_treated_unique_units/catt_total_unique_units as these will be added as part of (un)weighted means below

  # calculate unweighted and V1/V2 weighted means
  catt_treated_unique_units[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  catt_total_unique_units[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  figdata <- merge(figdata, catt_treated_unique_units, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
  figdata <- merge(figdata, catt_total_unique_units, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)

  subsets_for_avgs <- figdata[rn %in% c("catt")]
  subsets_for_avgs[, unweighted_estimate := mean(estimate, na.rm = TRUE), by = list(ref_event_time, model)]
  subsets_for_avgs[, cohort_weight_V1 := catt_treated_unique_units / sum(catt_treated_unique_units, na.rm = TRUE), by = list(ref_event_time, model)]
  subsets_for_avgs[, cohort_weight_V2 := catt_total_unique_units / sum(catt_total_unique_units, na.rm = TRUE), by = list(ref_event_time, model)]
  subsets_for_avgs[, weighted_estimate_V1 := weighted.mean(x = estimate, w = cohort_weight_V1, na.rm = TRUE), by = list(ref_event_time, model)]
  subsets_for_avgs[, weighted_estimate_V2 := weighted.mean(x = estimate, w = cohort_weight_V2, na.rm = TRUE), by = list(ref_event_time, model)]

  # merge weights into figdata
  weights <- subsets_for_avgs[, list(ref_event_time, ref_onset_time, model, cohort_weight_V1, cohort_weight_V2)]
  figdata <- merge(figdata, weights, by = c("ref_onset_time", "ref_event_time", "model"), all.x = TRUE, sort = FALSE)

  # rbind the (un)weighted avg estimates to figdata much like the "Pooled" data
  subsets_for_avgs[, rowid := seq_len(.N), by = list(ref_event_time, model)]
  subsets_for_avgs <- subsets_for_avgs[rowid == 1 | is.na(rowid)]

  unweighted <- subsets_for_avgs[, list(ref_event_time, unweighted_estimate, model)]
  unweighted[, ref_onset_time := "Equally-Weighted"]
  unweighted[, rn := "att"]
  setnames(unweighted, c("unweighted_estimate"), c("estimate"))
  setorderv(unweighted, c("ref_event_time", "model"))

  weighted_V1 <- subsets_for_avgs[, list(ref_event_time, weighted_estimate_V1, model)]
  weighted_V1[, ref_onset_time := "Cohort-Weighted"]
  weighted_V1[, rn := "att"]
  setnames(weighted_V1, c("weighted_estimate_V1"), c("estimate"))
  setorderv(weighted_V1, c("ref_event_time", "model"))

  weighted_V2 <- subsets_for_avgs[, list(ref_event_time, weighted_estimate_V2, model)]
  weighted_V2[, ref_onset_time := "Cohort-Weighted V2"]
  weighted_V2[, rn := "att"]
  setnames(weighted_V2, c("weighted_estimate_V2"), c("estimate"))
  setorderv(weighted_V2, c("ref_event_time", "model"))

  figdata <- rbindlist(list(figdata, unweighted, weighted_V1, weighted_V2), use.names = TRUE, fill=TRUE)
  figdata[is.na(cluster_se), cluster_se := 0]

  rm(subsets_for_avgs)
  rm(weights)
  rm(unweighted)
  rm(weighted_V1)
  rm(weighted_V2)
  rm(catt_treated_unique_units)
  rm(catt_total_unique_units)
  gc()

  # merge in the sample sizes
  figdata <- merge(figdata, event_time_total_unique_units, by = "ref_event_time", all.x = TRUE, sort = FALSE)
  rm(event_time_total_unique_units)

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
      dt[, unweighted_estimate := mean(estimate, na.rm = TRUE), by = list(ref_event_time, model)]
      dt[, weighted_estimate_V1 := weighted.mean(x = estimate, w = cohort_weight_V1, na.rm = TRUE), by = list(ref_event_time, model)]
      dt[, weighted_estimate_V2 := weighted.mean(x = estimate, w = cohort_weight_V2, na.rm = TRUE), by = list(ref_event_time, model)]
      dt[, rowid := seq_len(.N), by = list(ref_event_time, model)]
      result <- dt[rowid == 1 | is.na(rowid)]
      dt <- dt[, list(ref_event_time, ref_onset_time, model, cohort_weight_V1, cohort_weight_V2, grouping)]

      ew[[j]] <- result[, list(grouping, model, unweighted_estimate)]
      ew[[j]][, unweighted_estimate := mean(unweighted_estimate, na.rm = TRUE), by = list(grouping, model)]
      ew[[j]][, ref_onset_time := "Equally-Weighted + Collapsed"]
      ew[[j]][, rn := "att"]
      ew[[j]][, rowid := seq_len(.N), by = list(model)]
      setnames(ew[[j]], c("unweighted_estimate"), c("estimate"))
      ew[[j]] <- ew[[j]][rowid == 1 | is.na(rowid)]
      ew[[j]][, rowid := NULL]
      ew[[j]][, cluster_se := 0]

      cw[[j]] <- result[, list(grouping, model, weighted_estimate_V1)]
      cw[[j]][, weighted_estimate_V1 := mean(weighted_estimate_V1, na.rm = TRUE), by = list(grouping, model)]
      cw[[j]][, ref_onset_time := "Cohort-Weighted + Collapsed"]
      cw[[j]][, rn := "att"]
      cw[[j]][, rowid := seq_len(.N), by = list(model)]
      setnames(cw[[j]], c("weighted_estimate_V1"), c("estimate"))
      cw[[j]] <- cw[[j]][rowid == 1 | is.na(rowid)]
      cw[[j]][, rowid := NULL]
      cw[[j]][, cluster_se := 0]

      cw2[[j]] <- result[, list(grouping, model, weighted_estimate_V2)]
      cw2[[j]][, weighted_estimate_V2 := mean(weighted_estimate_V2, na.rm = TRUE), by = list(grouping, model)]
      cw2[[j]][, ref_onset_time := "Cohort-Weighted V2 + Collapsed"]
      cw2[[j]][, rn := "att"]
      cw2[[j]][, rowid := seq_len(.N), by = list(model)]
      setnames(cw2[[j]], c("weighted_estimate_V2"), c("estimate"))
      cw2[[j]] <- cw2[[j]][rowid == 1 | is.na(rowid)]
      cw2[[j]][, rowid := NULL]
      cw2[[j]][, cluster_se := 0]

      rm(dt, result)
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
    rm(ew, cw, cw2)
    rm(collapsed_estimate_total_unique_units)
    gc()
  }

  # add overall unique units count
  figdata[, total_unique_units := unique_units]
  rm(unique_units)

  # 'grouping' and'collapsed_estimate_total_unique_units' only relevant if calculate_collapse_estimates == TRUE, but easier to just generate a dummy columns if not
  if(calculate_collapse_estimates == FALSE){
    figdata[, grouping := NA]
    figdata[, collapsed_estimate_total_unique_units := NA]
  }

  # form ratios across models (rf/fs) by (ref_onset_time, ref_event_time, grouping)
  if(!is.null(endog_var)){

    ratio_data_rf <- figdata[rn %in% c("att", "catt") & model == "rf", list(ref_onset_time, ref_event_time, grouping, estimate, rn)]
    setnames(ratio_data_rf, "estimate", "estimate_rf")

    ratio_data_fs <- figdata[rn %in% c("att", "catt") & model == "fs", list(ref_onset_time, ref_event_time, grouping, estimate, rn)]
    setnames(ratio_data_fs, "estimate", "estimate_fs")

    ratio_data <- merge(ratio_data_rf, ratio_data_fs, by = c("ref_onset_time", "ref_event_time", "grouping", "rn"), sort = FALSE)
    rm(ratio_data_rf, ratio_data_fs)
    ratio_data[, estimate := (estimate_rf / estimate_fs)]
    ratio_data[, c("estimate_rf", "estimate_fs") := NULL]

    # merge in various counts to ratio data, then rbindlist back with figdata
    # counts and weights will be the same across 'model' by (ref_onset_time, ref_event_time, grouping), so just need one
    counts_and_weights <- figdata[rn %in% c("att", "catt") & model == "rf", list(ref_onset_time, ref_event_time, grouping, rn, reg_sample_size, catt_treated_unique_units, catt_total_unique_units, cohort_weight_V1, cohort_weight_V2, event_time_total_unique_units, collapsed_estimate_total_unique_units, total_unique_units)]
    ratio_data <- merge(ratio_data, counts_and_weights, by = c("ref_onset_time", "ref_event_time", "grouping", "rn"), sort = FALSE)
    ratio_data[, model := "ratio"]
    ratio_data[, cluster_se := 0]

    figdata <- rbindlist(list(figdata, ratio_data), use.names = TRUE, fill = TRUE)
    rm(ratio_data)
  }

  figdata[, bootstrap_sample := iter]

  # renaming to make the models crystal clear; model == "ratio" already pretty clear
  figdata[model == "rf", model := "reduced_form"]
  figdata[model == "fs", model := "first_stage"]

  flog.info(sprintf("Wald_bootstrap_ES iteration %s is finished.", format(iter, scientific = FALSE, big.mark = ",")))

  return(figdata)
}
