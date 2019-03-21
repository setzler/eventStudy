#' @export
ES <- function(long_data, outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars,
               omitted_event_time= -2, anticipation = 0, min_control_gap=1, max_control_gap=Inf, linearize_pretrends=FALSE,
               control_subset_var=NA, control_subset_event_time=0, fill_zeros=FALSE,
               residualize_covariates = FALSE, discrete_covars = NULL, cont_covars = NULL, never_treat_action = 'none',
               homogeneous_ATT = FALSE, num_cores = 1, reg_weights = NULL, require_balanced_control_diff = TRUE,
               bootstrapES = FALSE, bootstrap_iters = 1,
               ipw = FALSE, ipw_model = 'linear', ipw_composition_change = FALSE){

  flog.info("Beginning ES.")

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
  assertIntegerish(control_subset_event_time,len=1)
  assertFlag(linearize_pretrends)
  assertFlag(fill_zeros)
  assertFlag(residualize_covariates)
  if(residualize_covariates){
    if(!any(testCharacter(discrete_covars), testCharacter(cont_covars))){
      stop("Since residualize_covariates=TRUE, either discrete_covars or cont_covars must be provided as a character vector.")
    }
  }
  assertFlag(homogeneous_ATT)
  assertIntegerish(num_cores,len=1,lower=1)
  assertFlag(bootstrapES)
  assertIntegerish(bootstrap_iters,len=1,lower=1)
  assertFlag(ipw)
  assertFlag(ipw_composition_change)

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
  if(testCharacter(reg_weights)){
    if(!(reg_weights %in% names(long_data))){stop(sprintf("Variable reg_weights='%s' is not in the long_data you provided.",reg_weights))}
  }

  # check that control variables don't overlap with design variables (e.g., cal_time_var, and onset_time_var)
  design_vars <- c(cal_time_var, onset_time_var)
  if(testCharacter(discrete_covars)){
    for(vv in discrete_covars){
      if(vv %in% design_vars){stop(sprintf("Variable discrete_covars='%s' is among %s, which are already controlled in the design.",vv))}
    }
  }
  if(testCharacter(cont_covars)){
    for(vv in cont_covars){
      if(vv %in% design_vars){stop(sprintf("Variable cont_covars='%s' is among %s, which are already controlled in the design.",vv))}
    }
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
  if(is.null(cluster_vars)){warning(sprintf("Supplied cluster_vars = NULL; given stacking in ES(), standard errors may be too small. Consider cluster_vars='%s' instead.", unit_var))}

  # warning if supplied num_cores exceeds detectCores() - 1
  # not a perfect test (even as an upper bound) as results of detectCores() may be OS-dependent, do not respect cluster allocation limits, etc.
  # see help for detectCores() and mc.cores() for more information
  if(num_cores > (parallel::detectCores() - 1)){warning(sprintf("Supplied num_cores='%s'; this exceeds typical system limits and may cause issues.", num_cores))}

  # check that ipw model conforms to available options
  if(ipw == TRUE){
    if(!(ipw_model %in% c('linear', 'logit', 'probit'))){
      stop(sprintf("ipw_model='%s' is not among allowed values (c('linear', 'logit', 'probit')).", ipw_model))
    }
  }

  # if user will be producing bootstrap SEs, want to preserve a copy of the data as it is now;
  # otherwise, potential for changes to underlying sample in subsequent steps
  if(bootstrapES == TRUE){
    orig_sample <- copy(long_data)
  }

  # edit long_data in line with supplied never_treat_action option
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
  if(fill_zeros){
    flog.info("Filling in zeros.")
    long_data <- ES_expand_to_balance(long_data = long_data,
                                      vars_to_fill = outcomevar,
                                      unit_var = unit_var,
                                      cal_time_var = cal_time_var,
                                      onset_time_var = onset_time_var)
  }

  # Check that there exist cohorts with observations at omitted_event_time
  if(is.infinite(suppressWarnings(long_data[get(cal_time_var) - get(onset_time_var) == omitted_event_time, min(get(onset_time_var))]))){
    stop(sprintf("Variable onset_time_var='%s' has no treated groups with observations at pre-treatment event time %s.",onset_time_var, omitted_event_time))
  }

  # linearize pre-trends; never-treated will be treated as a single cohort
  if(linearize_pretrends){
    flog.info("Linearizing pre-trends.")
    long_data <- ES_parallelize_trends(long_data = long_data, outcomevar = outcomevar,
                                       unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                                       anticipation = anticipation, reg_weights = reg_weights)
  }

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
                           never_treat_action = never_treat_action, never_treat_val = never_treat_val,
                           cluster_vars = cluster_vars, discrete_covars = discrete_covars, cont_covars = cont_covars, reg_weights = reg_weights)

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
                                 omitted_event_time = omitted_event_time,
                                 ipw_model = ipw_model,
                                 reg_weights = reg_weights,
                                 ipw_composition_change = ipw_composition_change
        )
        ipw_dt[, did_id := i]

        ES_data <- merge(ES_data, ipw_dt, by = c(unit_var, cal_time_var, "did_id"), all.x = TRUE, sort = FALSE)
        ES_data[is.na(pr_temp) & !is.na(pr), pr_temp := pr]
        ES_data[, pr := NULL]
        gc()
      }

      setnames(ES_data, "pr_temp", "pr")

    }
  }

  # collect ATT estimates
  if(homogeneous_ATT == FALSE){
    ES_results_hetero <- ES_estimate_ATT(ES_data = ES_data,
                                         outcomevar=outcomevar,
                                         onset_time_var = onset_time_var,
                                         cluster_vars = cluster_vars,
                                         homogeneous_ATT = homogeneous_ATT,
                                         omitted_event_time = omitted_event_time,
                                         discrete_covars = discrete_covars,
                                         cont_covars = cont_covars,
                                         residualize_covariates = residualize_covariates,
                                         reg_weights = reg_weights,
                                         ipw = ipw,
                                         ipw_composition_change = ipw_composition_change)

    catt_coefs <- ES_results_hetero[[2]]
    catt_vcov <- ES_results_hetero[[3]]

    ES_results_hetero <- ES_results_hetero[[1]]

    setnames(ES_results_hetero,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))
  } else{
    ES_results_hetero = NULL
  }
  ES_results_homo <- ES_estimate_ATT(ES_data = ES_data,
                                     outcomevar=outcomevar,
                                     onset_time_var = onset_time_var,
                                     cluster_vars = cluster_vars,
                                     homogeneous_ATT = TRUE,
                                     omitted_event_time = omitted_event_time,
                                     discrete_covars = discrete_covars,
                                     cont_covars = cont_covars,
                                     residualize_covariates = residualize_covariates,
                                     reg_weights = reg_weights,
                                     ipw = ipw,
                                     ipw_composition_change = ipw_composition_change)[[1]]
  setnames(ES_results_homo,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

  # collect levels by treatment/control
  if(!(is.null(reg_weights))){
    # STILL NEED TO FIX THE SEs to be WEIGHTEDs; for now use unweighted
    ES_treatcontrol_means <- ES_data[,list(rn="treatment_means", estimate = weighted.mean(get(outcomevar), get(reg_weights)), cluster_se = sd(get(outcomevar))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]
  } else{
    ES_treatcontrol_means <- ES_data[,list(rn="treatment_means", estimate = mean(get(outcomevar)), cluster_se = sd(get(outcomevar))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]
  }

  # collect count of treated units by each (ref_onset_time, ref_event_time) for V1 of population-weighted ATTs
  if(!(is.null(reg_weights))){
    # to revisit
    ES_treat_count_V1 <- ES_data[treated == 1,list(treat_count_V1 = sum(get(reg_weights))), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  } else{
    ES_treat_count_V1 <- ES_data[treated == 1,list(treat_count_V1 = .N), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  }
  # omitted_event_time is excluded, so let's ensure it is excluded from the above as well
  ES_treat_count_V1 <- ES_treat_count_V1[ref_event_time != omitted_event_time]

  # collect count of treated+control units by each (ref_onset_time, ref_event_time) for V2 of population-weighted ATTs
  if(!(is.null(reg_weights))){
    # to revisit
    ES_treat_count_V2 <- ES_data[,list(treat_count_V2 = sum(get(reg_weights))), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  } else{
    ES_treat_count_V2 <- ES_data[,list(treat_count_V2 = .N), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  }
  # omitted_event_time is excluded, so let's ensure it is excluded from the above as well
  ES_treat_count_V2 <- ES_treat_count_V2[ref_event_time != omitted_event_time]

  ES_data <- NULL
  gc()

  figdata <- rbindlist(list(ES_results_hetero, ES_results_homo, ES_treatcontrol_means), use.names = TRUE, fill=TRUE)

  # calculate unweighted and V1/V2 weighted means
  ES_treat_count_V1[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  ES_treat_count_V2[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  figdata <- merge(figdata, ES_treat_count_V1, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
  figdata <- merge(figdata, ES_treat_count_V2, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)

  subsets_for_avgs <- figdata[rn %in% c("catt")]
  subsets_for_avgs[, unweighted_estimate := mean(estimate, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weight_V1 := treat_count_V1 / sum(treat_count_V1, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weight_V2 := treat_count_V2 / sum(treat_count_V2, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weighted_estimate_V1 := weighted.mean(x = estimate, w = weight_V1, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weighted_estimate_V2 := weighted.mean(x = estimate, w = weight_V2, na.rm = TRUE), by = list(ref_event_time)]

  # merge weights into figdata
  weights <- subsets_for_avgs[, list(ref_event_time, ref_onset_time, weight_V1, weight_V2)]
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

  # calculate SEs for weighted avg estimates using catt_coefs and catt_vcov
  # the unique "Weighted" ref_event_time values represent the target list of parameters
  # for each such ref_event_time, want to extract the location (number) of the relevant parameters in catt_coefs
  # then will need to grab the relevant weight, and then construct the formula to supply to deltamethod()

  event_times <- setdiff(figdata[, sort(unique(ref_event_time))], omitted_event_time)
  onset_times <- as.integer(figdata[rn == "catt", sort(unique(ref_onset_time))])

  min_onset_time <- min(onset_times)
  max_onset_time <- max(onset_times)

  for(et in event_times){

    if(et < 0){
      lookfor <- sprintf("cattlead%s", abs(et))
    } else{
      lookfor <- sprintf("catt%s", abs(et))
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
    temp <- temp[, list(ref_onset_time, ref_event_time, coef_index, weight_V1, weight_V2)]
    temp[, weight_V0 := 1 / .N]

    temp[, equal_w_formula_entry := sprintf("(%s*x%s)", weight_V0, coef_index)]
    temp[, cohort_w_v1_formula_entry := sprintf("(%s*x%s)", weight_V1, coef_index)]
    temp[, cohort_w_v2_formula_entry := sprintf("(%s*x%s)", weight_V2, coef_index)]

    equal_w_g_formula_input = paste0(temp$equal_w_formula_entry, collapse = "+")
    cohort_w_v1_g_formula_input = paste0(temp$cohort_w_v1_formula_entry, collapse = "+")
    cohort_w_v2_g_formula_input = paste0(temp$cohort_w_v2_formula_entry, collapse = "+")

    figdata[rn == "att" & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Equally-Weighted",
            cluster_se := deltamethod(g = as.formula(paste("~", equal_w_g_formula_input)),
                                      mean = catt_coefs, cov = catt_vcov, ses = TRUE
            )
            ]

    figdata[rn == "att" & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Cohort-Weighted",
            cluster_se := deltamethod(g = as.formula(paste("~", cohort_w_v1_g_formula_input)),
                                      mean = catt_coefs, cov = catt_vcov, ses = TRUE
            )
            ]

    figdata[rn == "att" & cluster_se == 0 & ref_event_time == et & ref_onset_time == "Cohort-Weighted V2",
            cluster_se := deltamethod(g = as.formula(paste("~", cohort_w_v2_g_formula_input)),
                                      mean = catt_coefs, cov = catt_vcov, ses = TRUE
            )
            ]

    rm(temp, equal_w_g_formula_input, cohort_w_v1_g_formula_input, cohort_w_v2_g_formula_input)

  }

  subsets_for_avgs <- NULL
  weights <- NULL
  unweighted <- NULL
  weighted_V1 <- NULL
  weighted_V2 <- NULL
  gc()

  # start bootstrap run if relevant
  if(bootstrapES == TRUE){

    boot_results <- rbindlist(parallel::mclapply(X = 1:bootstrap_iters,
                                                 FUN = bootstrap_ES,
                                                 mc.silent = FALSE,
                                                 mc.cores = num_cores,
                                                 mc.set.seed = TRUE,
                                                 long_data = orig_sample, outcomevar = outcomevar, unit_var = unit_var, cal_time_var = cal_time_var,
                                                 onset_time_var = onset_time_var, cluster_vars = cluster_vars, omitted_event_time= omitted_event_time,
                                                 anticipation = anticipation, min_control_gap = min_control_gap, max_control_gap = max_control_gap,
                                                 linearize_pretrends = linearize_pretrends, control_subset_var = control_subset_var, control_subset_event_time = control_subset_event_time,
                                                 fill_zeros = fill_zeros, residualize_covariates = residualize_covariates, discrete_covars = discrete_covars, cont_covars = cont_covars,
                                                 never_treat_action = never_treat_action, homogeneous_ATT = homogeneous_ATT, reg_weights = reg_weights,
                                                 require_balanced_control_diff = require_balanced_control_diff,
                                                 ipw = ipw, ipw_model = ipw_model, ipw_composition_change = ipw_composition_change),
                              use.names = TRUE
                              )

    # Construct bootstrap CI for all of the estimates
    boot_results_ses <- boot_results[, list(bootstrap_se = sd(estimate)), by = list(ref_onset_time, ref_event_time)]
    figdata <- merge(figdata, boot_results_ses, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)

    orig_sample <- NULL
    gc()

  }

  flog.info('ES is finished.')

  return(figdata)
}



#' @export
ES_plot_ATTs <- function(figdata, lower_event = -3, upper_event = 5, ci_factor = 1.96, homogeneous_ATT = FALSE, omitted_event_time = -2, bootstrap_ses = FALSE){

  # As this function changes the underlying results data and figdata tends to be a small file, make a copy now
  figdata_temp <- copy(figdata)

  figdata_temp <- figdata_temp[rn %in% c("att","catt")]
  figdata_temp[, ref_event_time := as.numeric(ref_event_time)]
  figdata_temp <- figdata_temp[ ref_event_time >= lower_event & ref_event_time <= upper_event]
  figdata_temp[, jitter := .GRP, by = ref_onset_time]
  jitter_center <- figdata_temp[,median(unique(jitter))]
  jitter_scale <- figdata_temp[,length(unique(jitter))]
  figdata_temp[, jitter_event_time := ref_event_time + (jitter - jitter_center) / jitter_scale]

  if(bootstrap_ses == TRUE){
    if(!("bootstrap_se" %in% colnames(figdata_temp))){
      stop(print("Variable 'bootstrap_se' is not found in results data you provided. Let me suggest bootstrap_ses = FALSE."))
    } else{
      figdata_temp[, cluster_se := NULL]
      setnames(figdata_temp, "bootstrap_se", "cluster_se")
    }
  }

  if(homogeneous_ATT){
    fig <- ggplot(aes( x = ref_event_time, y = estimate), data = figdata_temp[ref_onset_time=="Pooled"]) +
      geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - ci_factor * cluster_se, ymax = estimate + ci_factor * cluster_se)) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time") +
      annotate(geom = "text", x = omitted_event_time, y = 0, label = "(Omitted)", size = 4)
  } else {
    fig <- ggplot(aes( x = jitter_event_time, y = estimate, colour = factor(ref_onset_time)), data = figdata_temp) +
      geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - ci_factor * cluster_se, ymax = estimate + ci_factor * cluster_se)) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time", color = "Cohort") +
      annotate(geom = "text", x = omitted_event_time, y = 0, label = "(Omitted)", size = 4)
  }

  return(fig)

}

#' @export
ES_plot_levels <- function(figdata, cohort_subset = NA, lower_event = -3, upper_event = 5, omitted_event_time = -2){

  # As this function changes the underlying results data and figdata tends to be a small file, make a copy now
  figdata_temp <- copy(figdata)

  figdata_temp <- figdata_temp[rn %in% c("treatment_means")]
  figdata_temp[, estimate := estimate - estimate[ref_event_time == omitted_event_time], list(ref_onset_time,treated)]
  figdata_temp[, ref_onset_time :=  as.integer(as.character(ref_onset_time))]
  figdata_temp[, ref_event_time := as.integer(as.character(ref_event_time))]
  figdata_temp[treated==1, treatment := "Treatment"]
  figdata_temp[treated==0, treatment := "Control"]
  figdata_temp[, treatment := factor(treatment,levels=c("Treatment","Control"))]
  figdata_temp[, year := ref_onset_time + ref_event_time]
  figdata_temp <- figdata_temp[ref_event_time >= lower_event & ref_event_time <= upper_event]

  if(!is.na(cohort_subset)){
    figdata_temp <- figdata_temp[ref_onset_time %in% cohort_subset]
  }

  gg <- ggplot(aes(x=year,y=estimate,colour=factor(ref_onset_time),linetype=treatment),data=figdata_temp) +
    geom_line() + geom_point() + geom_vline(aes(xintercept=ref_onset_time,colour=factor(ref_onset_time)),linetype=3) +
    theme_bw(base_size = 16) + labs(colour="Cohort",linetype="",x="Year") +
    scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

  return(gg)

}

#' @export
block_sample <- function(long_data, unit_var, cal_time_var){
  # Thanks (without implicating) to Michael Graber for sharing this function
  # vector of unique cluster IDs and total number of clusters
  setkeyv(long_data, c(unit_var, cal_time_var))
  ids <- unique(long_data[, get(unit_var)])
  N   <- length(ids)
  # sample ID's with replacement:
  # we sample ID's with replacement, then collapse the resulting draws to
  # a vector of unique IDs with a column representing the number of times the cluster
  # has been sampled. e.g. if unit_var = c("A","B") and count = c(1,3) then
  # cluster "A" has been sampled once, while cluster "B" has been sampled three times.
  clusters <- data.table(sample(ids, size = N, replace = TRUE))
  names(clusters) <- unit_var
  setkeyv(clusters, unit_var)
  clusters[,count := .N, by = get(unit_var)]
  clusters <- unique(clusters)
  # repeated merges to create block bootstrap sample (bs)
  max_iter <- max(clusters$count)
  bs       <- data.table()
  for (i in 1 : max_iter) {
    tmp <- merge(clusters[count >= i,], long_data, by = unit_var, all.x = TRUE,
                 all.y = FALSE, sort = FALSE)
    tmp[, eval(unit_var) := paste0(get(unit_var), "-", as.character(i))]  # create new unique ID: original ID "-" index of draw
    bs <- rbind(bs, tmp)
  }
  bs[,count := NULL]
  setkeyv(bs, c(unit_var, cal_time_var))
  gc()

  return(bs)
}

#' @export
bootstrap_ES <- function(long_data, outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars,
                         omitted_event_time= -2, anticipation = 0, min_control_gap=1, max_control_gap=Inf, linearize_pretrends=FALSE,
                         control_subset_var=NA, control_subset_event_time=0, fill_zeros=FALSE,
                         residualize_covariates = FALSE, discrete_covars = NULL, cont_covars = NULL, never_treat_action = 'none',
                         homogeneous_ATT = FALSE, reg_weights = NULL, require_balanced_control_diff = TRUE,
                         ipw = FALSE, ipw_model = 'linear', ipw_composition_change = FALSE,
                         iter){


  bs_long_data <- block_sample(long_data = copy(long_data), unit_var = unit_var, cal_time_var = cal_time_var)

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
  assertIntegerish(control_subset_event_time,len=1)
  assertFlag(linearize_pretrends)
  assertFlag(fill_zeros)
  assertFlag(residualize_covariates)
  if(residualize_covariates){
    if(!any(testCharacter(discrete_covars), testCharacter(cont_covars))){
      stop("Since residualize_covariates=TRUE, either discrete_covars or cont_covars must be provided as a character vector.")
    }
  }
  assertFlag(homogeneous_ATT)
  assertFlag(ipw)
  assertFlag(ipw_composition_change)

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
      if(!(vv %in% names(long_data))){stop(sprintf("Variable cluster_vars='%s' is not in the long_data you provided. Let me suggest cluster_vars='%s'.",vv,unit_var))}
    }
  }
  if(!is.na(control_subset_var)){
    if(!(control_subset_var %in% names(bs_long_data))){stop(sprintf("Variable control_subset_var='%s' is not in the bs_long_data you provided.",control_subset_var))}
    if(!(bs_long_data[,typeof(get(control_subset_var))]=="logical")){stop(sprintf("Variable control_subset_var='%s' must be of type logical (i.e., only TRUE or FALSE values).",control_subset_var))}
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
  if(testCharacter(reg_weights)){
    if(!(reg_weights %in% names(bs_long_data))){stop(sprintf("Variable reg_weights='%s' is not in the bs_long_data you provided.",reg_weights))}
  }

  # check that control variables don't overlap with design variables (e.g., cal_time_var, and onset_time_var)
  design_vars <- c(cal_time_var, onset_time_var)
  if(testCharacter(discrete_covars)){
    for(vv in discrete_covars){
      if(vv %in% design_vars){stop(sprintf("Variable discrete_covars='%s' is among %s, which are already controlled in the design.",vv))}
    }
  }
  if(testCharacter(cont_covars)){
    for(vv in cont_covars){
      if(vv %in% design_vars){stop(sprintf("Variable cont_covars='%s' is among %s, which are already controlled in the design.",vv))}
    }
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

  # check that ipw model conforms to available options
  if(ipw == TRUE){
    if(!(ipw_model %in% c('linear', 'logit', 'probit'))){
      stop(sprintf("ipw_model='%s' is not among allowed values (c('linear', 'logit', 'probit')).", ipw_model))
    }
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
    flog.info(sprintf("Filling in zeros ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
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

  # linearize pre-trends; never-treated will be treated as a single cohort
  if(linearize_pretrends){
    flog.info(sprintf("Linearizing pre-trends in ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
    bs_long_data <- ES_parallelize_trends(long_data = bs_long_data, outcomevar = outcomevar,
                                          unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                                          anticipation = anticipation, reg_weights = reg_weights)
  }

  if(residualize_covariates){
    flog.info(sprintf("Residualizing on covariates in ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
    bs_long_data <- ES_residualize_covariates(long_data = bs_long_data, outcomevar = outcomevar,
                                              unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                                              anticipation = anticipation, discrete_covars = discrete_covars, cont_covars = cont_covars,
                                              reg_weights = reg_weights)
  }

  # process data
  flog.info(sprintf("Beginning data stacking in ES_bootstrap iteration %s.", format(iter, scientific = FALSE, big.mark = ",")))
  ES_data <- ES_clean_data(long_data = bs_long_data, outcomevar = outcomevar,
                           unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                           anticipation = anticipation, min_control_gap = min_control_gap, max_control_gap = max_control_gap, omitted_event_time = omitted_event_time,
                           control_subset_var = control_subset_var, control_subset_event_time = control_subset_event_time,
                           never_treat_action = never_treat_action, never_treat_val = never_treat_val,
                           cluster_vars = NULL, discrete_covars = discrete_covars, cont_covars = cont_covars, reg_weights = reg_weights)

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
                                 omitted_event_time = omitted_event_time,
                                 ipw_model = ipw_model,
                                 reg_weights = reg_weights,
                                 ipw_composition_change = ipw_composition_change
        )
        ipw_dt[, did_id := i]

        ES_data <- merge(ES_data, ipw_dt, by = c(unit_var, cal_time_var, "did_id"), all.x = TRUE, sort = FALSE)
        ES_data[is.na(pr_temp) & !is.na(pr), pr_temp := pr]
        ES_data[, pr := NULL]
        gc()
      }

      setnames(ES_data, "pr_temp", "pr")

    }
  }

  # collect ATT estimates
  if(homogeneous_ATT == FALSE){
    ES_results_hetero <- ES_estimate_ATT(ES_data = ES_data,
                                         outcomevar=outcomevar,
                                         onset_time_var = onset_time_var,
                                         cluster_vars = NULL,
                                         homogeneous_ATT = homogeneous_ATT,
                                         omitted_event_time = omitted_event_time,
                                         discrete_covars = discrete_covars,
                                         cont_covars = cont_covars,
                                         residualize_covariates = residualize_covariates,
                                         reg_weights = reg_weights,
                                         ipw = ipw,
                                         ipw_composition_change = ipw_composition_change)[[1]]
    setnames(ES_results_hetero,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))
  } else{
    ES_results_hetero = NULL
  }
  ES_results_homo <- ES_estimate_ATT(ES_data = ES_data,
                                     outcomevar=outcomevar,
                                     onset_time_var = onset_time_var,
                                     cluster_vars = NULL,
                                     homogeneous_ATT = TRUE,
                                     omitted_event_time = omitted_event_time,
                                     discrete_covars = discrete_covars,
                                     cont_covars = cont_covars,
                                     residualize_covariates = residualize_covariates,
                                     reg_weights = reg_weights,
                                     ipw = ipw,
                                     ipw_composition_change = ipw_composition_change)[[1]]
  setnames(ES_results_homo,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

  # collect count of treated units by each (ref_onset_time, ref_event_time) for V1 of population-weighted ATTs
  if(!(is.null(reg_weights))){
    # to revisit
    ES_treat_count_V1 <- ES_data[treated == 1,list(treat_count_V1 = sum(get(reg_weights))), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  } else{
    ES_treat_count_V1 <- ES_data[treated == 1,list(treat_count_V1 = .N), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  }
  # omitted_event_time is excluded, so let's ensure it is excluded from the above as well
  ES_treat_count_V1 <- ES_treat_count_V1[ref_event_time != omitted_event_time]

  # collect count of treated+control units by each (ref_onset_time, ref_event_time) for V2 of population-weighted ATTs
  if(!(is.null(reg_weights))){
    # to revisit
    ES_treat_count_V2 <- ES_data[,list(treat_count_V2 = sum(get(reg_weights))), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  } else{
    ES_treat_count_V2 <- ES_data[,list(treat_count_V2 = .N), by = list(ref_onset_time,ref_event_time)][order(ref_onset_time,ref_event_time)]
  }
  # omitted_event_time is excluded, so let's ensure it is excluded from the above as well
  ES_treat_count_V2 <- ES_treat_count_V2[ref_event_time != omitted_event_time]

  ES_data <- NULL
  gc()

  figdata <- rbindlist(list(ES_results_hetero, ES_results_homo), use.names = TRUE, fill=TRUE)

  # calculate unweighted and V1/V2 weighted means
  ES_treat_count_V1[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  ES_treat_count_V2[, ref_onset_time := as.character(ref_onset_time)] # to match figdata
  figdata <- merge(figdata, ES_treat_count_V1, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)
  figdata <- merge(figdata, ES_treat_count_V2, by = c("ref_onset_time", "ref_event_time"), all.x = TRUE, sort = FALSE)

  subsets_for_avgs <- figdata[rn %in% c("catt")]
  subsets_for_avgs[, unweighted_estimate := mean(estimate, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weight_V1 := treat_count_V1 / sum(treat_count_V1, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weight_V2 := treat_count_V2 / sum(treat_count_V2, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weighted_estimate_V1 := weighted.mean(x = estimate, w = weight_V1, na.rm = TRUE), by = list(ref_event_time)]
  subsets_for_avgs[, weighted_estimate_V2 := weighted.mean(x = estimate, w = weight_V2, na.rm = TRUE), by = list(ref_event_time)]

  # merge weights into figdata
  weights <- subsets_for_avgs[, list(ref_event_time, ref_onset_time, weight_V1, weight_V2)]
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

  subsets_for_avgs <- NULL
  weights <- NULL
  unweighted <- NULL
  weighted_V1 <- NULL
  weighted_V2 <- NULL
  gc()

  figdata[, bootstrap_sample := iter]

  flog.info(sprintf("ES_bootstrap iteration %s is finished.", format(iter, scientific = FALSE, big.mark = ",")))

  return(figdata)
}



