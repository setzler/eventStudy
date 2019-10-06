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
