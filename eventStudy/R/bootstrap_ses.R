#' @export
# Need 'iter' to come first given how R reads the lapply/mclapply implementing bootstrap_ES()
bootstrap_ES <- function(iter, long_data, outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars,
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
                         ntile_var = NULL, ntile_event_time = -2, ntiles = NA, ntile_var_value = NA, ntile_avg = FALSE){

  # Function mirrors ES(), but with less printing/warning and skipping sections responsible for analytical standard errors and the like which aren't relevant
  # We keep various checks to make sure no issues arise after the bootstrap sample is drawn

  bs_long_data <- block_sample(long_data = copy(long_data), unit_var = unit_var, cal_time_var = cal_time_var)
  gc()

  flog.info(sprintf("Beginning ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))

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
    stop("Supplied variables for both reg_weights and ref_reg_weights, but ES() only admits using one or the other type of weight (or neither).")
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
    flog.info(sprintf("Filling in zeros in ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
    bs_long_data <- ES_expand_to_balance(long_data = bs_long_data,
                                         vars_to_fill = outcomevar,
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
    flog.info(sprintf("Linearizing pre-trends in ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
    bs_long_data <- ES_parallelize_trends(long_data = bs_long_data, outcomevar = outcomevar,
                                          unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                                          anticipation = anticipation, reg_weights = reg_weights)
  }

  # extrapolate contribution of covariates using pre-treated period
  if(residualize_covariates){
    flog.info(sprintf("Residualizing on covariates in ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
    bs_long_data <- ES_residualize_covariates(long_data = bs_long_data, outcomevar = outcomevar,
                                              unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                                              anticipation = anticipation, discrete_covars = discrete_covars, cont_covars = cont_covars,
                                              reg_weights = reg_weights)
  }

  # process data
  flog.info(sprintf("Beginning data stacking in ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
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
                              ntile_var = ntile_var, ntile_event_time = ntile_event_time, ntiles = ntiles, ntile_var_value = ntile_var_value, ntile_avg = ntile_avg)

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
        bs_ES_data[, (varname) := sum(as.integer(get(varname)*(ref_event_time==ref_discrete_covar_event_time)), na.rm = TRUE), by=c(unit_var, "ref_onset_time")]

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
        bs_ES_data[, (varname) := sum(as.integer(get(varname)*(ref_event_time==ref_discrete_covar_event_time)), na.rm = TRUE), by=c(unit_var, "ref_onset_time")]

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
          bs_ES_data[, (varname) := sum(as.integer(get(var)*(ref_event_time==ref_cont_covar_event_time)), na.rm = TRUE), by=c(unit_var, "ref_onset_time")]
        } else{
          bs_ES_data[, (varname) := sum(get(var)*(ref_event_time==ref_cont_covar_event_time), na.rm = TRUE), by=c(unit_var, "ref_onset_time")]
        }
      } else{

        # Define the time-invariant outcome for all units within a ref_onset_time
        if(class(bs_ES_data[[var]]) == "integer"){
          # needed to include type checks as sum() often converts from integer to numeric
          # use sum() instead of max() in case 'var' takes on negative values (as max() would then return 0)
          bs_ES_data[, (var) := sum(as.integer(get(var)*(ref_event_time==ref_cont_covar_event_time)), na.rm = TRUE), by=c(unit_var, "ref_onset_time")]
        } else{
          bs_ES_data[, (var) := sum(get(var)*(ref_event_time==ref_cont_covar_event_time), na.rm = TRUE), by=c(unit_var, "ref_onset_time")]
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
      bs_ES_data[, (ref_reg_weights) := sum(as.integer(get(ref_reg_weights)*(ref_event_time==ref_reg_weights_event_time)), na.rm = TRUE), by=c(unit_var, "ref_onset_time")]
    } else{
      bs_ES_data[, (ref_reg_weights) := sum(get(ref_reg_weights)*(ref_event_time==ref_reg_weights_event_time), na.rm = TRUE), by=c(unit_var, "ref_onset_time")]
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
    ES_results_hetero <- ES_estimate_ATT(ES_data = copy(bs_ES_data),
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
                                         add_unit_fes = add_unit_fes)[[1]]
    gc()
    setnames(ES_results_hetero,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))
  } else{
    ES_results_hetero <- NULL
  }
  ES_results_homo <- ES_estimate_ATT(ES_data = copy(bs_ES_data),
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
                                     add_unit_fes = add_unit_fes)[[1]]
  gc()
  setnames(ES_results_homo,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

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

  figdata <- rbindlist(list(ES_results_hetero, ES_results_homo), use.names = TRUE, fill=TRUE)

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
    }
    rm(j)
    ew <- rbindlist(ew, use.names = TRUE)
    cw <- rbindlist(cw, use.names = TRUE)
    cw2 <- rbindlist(cw2, use.names = TRUE)

    figdata <- rbindlist(list(figdata, ew, cw, cw2), use.names = TRUE, fill = TRUE)
    gc()
  }

  # add overall unique units count
  figdata[, total_unique_units := unique_units]

  figdata[, bootstrap_sample := iter]

  flog.info(sprintf("ES_bootstrap iteration %s is finished.", format(iter, scientific = FALSE, big.mark = ",")))

  return(figdata)
}
