

#' @export
ES_estimate_ATT <- function(ES_data,
                            outcomevar,
                            onset_time_var,
                            cluster_vars,
                            residualize_covariates = FALSE,
                            discrete_covars = NULL,
                            cont_covars = NULL,
                            homogeneous_ATT = TRUE,
                            omitted_event_time = -2,
                            reg_weights = NULL,
                            ipw = FALSE,
                            ipw_composition_change = FALSE
                            ) {

  onset_times <- sort(unique(ES_data[, .N, by = eval(onset_time_var)][[onset_time_var]]))
  min_onset_time <- min(onset_times)
  max_onset_time <- max(onset_times)

  ES_data[, ref_onset_ref_event_time := .GRP, by = list(ref_onset_time, ref_event_time)]
  ES_data[, unit_sample := .GRP, by = list(get(onset_time_var), catt_specific_sample, ref_onset_time)]
  ES_data[, cluster_on_this := .GRP, by = cluster_vars]

  if(residualize_covariates == FALSE & ipw == FALSE){

    if(!is.null(discrete_covars)){

      # For felm(), want to tell which factor combination will be sent to the intercept
      # Let's set it to the min

      i <- 0
      reference_lookup <- list()
      discrete_covar_formula_input <- c()
      for(var in discrete_covars){
        i <- i + 1
        min_pre_value <- ES_data[, min(get(var))]
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

    if(!is.null(cont_covars)){
      cont_covar_formula_input = paste0(na.omit(cont_covars), collapse = "+")
    } else{
      cont_covar_formula_input = NA
    }

  }

  start_cols <- copy(colnames(ES_data))

  if (homogeneous_ATT == FALSE) {

    for (h in min(ES_data$ref_onset_time):max(ES_data$ref_onset_time)) {
      for (r in setdiff(min(ES_data[ref_onset_time == h]$ref_event_time):max(ES_data[ref_onset_time == h]$ref_event_time), omitted_event_time)) {
        var <- sprintf("ref_onset_time%s_catt%s", h, r)
        ES_data[, (var) := as.integer(ref_onset_time == h & ref_event_time == r & get(onset_time_var) == h)]
      }
    }

    new_cols <- setdiff(colnames(ES_data), start_cols)
    new_cols_used <- gsub("\\-", "lead", new_cols) # "-" wreaks havoc otherwise, syntactically
    setnames(ES_data, new_cols, new_cols_used)

    # Should now be able to combine all of the above in a natural way
    felm_formula_input <- paste(new_cols_used, collapse = "+")

    if(ipw == TRUE){

      if(ipw_composition_change == FALSE){
        # Construct ATT weights and estimate WLS
        # For ipw_composition_change == TRUE case, make them as part of ES_clean_data

        ES_data[treated == 1, weight := 1]
        ES_data[treated == 0, weight := (pr / (1 - pr))]

        # This would be the place to add additional assumptions
        # E.g., restrict to propensity scores strictly between 0 and 1

        ES_data[, pr_subset := as.integer(between(pr, 0, 1, incbounds = FALSE))]
        gc()
      } else{
        ES_data[, pr_subset := 1L]
      }

      # Pre-processing some issues of numerical precision
      ES_data[weight < 0, weight := 0]

      est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | cluster_on_this")),
                  data = ES_data[pr_subset == 1], weights = ES_data[pr_subset == 1]$weight,
                  nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
      )

    } else if(residualize_covariates == TRUE){

      if(!(is.null(reg_weights))){
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | cluster_on_this")),
                    data = ES_data, weights = ES_data[[reg_weights]],
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      } else{
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | cluster_on_this")),
                    data = ES_data,
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      }

    } else{

      felm_formula_input = paste(na.omit(c(felm_formula_input, cont_covar_formula_input)), collapse = " + ")
      fe_formula = paste(na.omit(c("unit_sample + ref_onset_ref_event_time", discrete_covar_formula_input)), collapse = " + ")

      if(!(is.null(reg_weights))){
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | ", fe_formula," | 0 | cluster_on_this")),
                    data = ES_data, weights = ES_data[[reg_weights]],
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      } else{
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | ", fe_formula," | 0 | cluster_on_this")),
                    data = ES_data,
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      }
    }

    ES_data <- NULL
    gc()

    results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
    est <- NULL
    gc()
    results[!grepl("ref\\_onset\\_time", rn), e := min_onset_time]
    results[, rn := gsub("lead", "-", rn)]
    for (c in min_onset_time:(max_onset_time - 1)) {
      results[grepl(sprintf("ref\\_onset\\_time%s", c), rn), e := c]
      results[grepl(sprintf("ref\\_onset\\_time%s", c), rn), rn := gsub(sprintf("ref\\_onset\\_time%s\\_et", c), "et", rn)]
      results[grepl(sprintf("ref\\_onset\\_time%s", c), rn), rn := gsub(sprintf("ref\\_onset\\_time%s\\_catt", c), "catt", rn)]
    }
    results[grepl("et", rn), event_time := as.integer(gsub("et", "", rn))]
    results[grepl("catt", rn), event_time := as.integer(gsub("catt", "", rn))]
    results[grepl("et", rn), rn := "event_time"]
    results[grepl("catt", rn), rn := "catt"]
    results <- results[rn == "catt"]
  } else {

    for (r in setdiff(min(ES_data$ref_event_time):max(ES_data$ref_event_time), omitted_event_time)) {
      var <- sprintf("att%s", r)
      ES_data[, (var) := as.integer(ref_event_time == r & ref_onset_time == get(onset_time_var))]
    }

    new_cols <- setdiff(colnames(ES_data), start_cols)
    new_cols_used <- gsub("\\-", "lead", new_cols) # "-" wreaks havoc otherwise, syntactically
    setnames(ES_data, new_cols, new_cols_used)

    # Should now be able to combine all of the above in a natural way
    felm_formula_input <- paste(c(new_cols_used), collapse = "+")

    if(ipw == TRUE){

      if(ipw_composition_change == FALSE){
        # Construct ATT weights and estimate WLS
        # For ipw_composition_change == TRUE case, make them as part of ES_clean_data

        ES_data[treated == 1, weight := 1]
        ES_data[treated == 0, weight := (pr / (1 - pr))]

        # This would be the place to add additional assumptions
        # E.g., restrict to propensity scores strictly between 0 and 1

        ES_data[, pr_subset := as.integer(between(pr, 0, 1, incbounds = FALSE))]
        gc()
      } else{
        ES_data[, pr_subset := 1L]
      }

      # Pre-processing some issues of numerical precision
      ES_data[weight < 0, weight := 0]

      est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | cluster_on_this")),
                  data = ES_data[pr_subset == 1], weights = ES_data[pr_subset == 1]$weight,
                  nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
      )

    } else if(residualize_covariates == TRUE){

      if(!(is.null(reg_weights))){
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | cluster_on_this")),
                    data = ES_data, weights = ES_data[[reg_weights]],
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      } else{
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | cluster_on_this")),
                    data = ES_data,
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      }

    } else{

      felm_formula_input = paste(na.omit(c(felm_formula_input, cont_covar_formula_input)), collapse = " + ")
      fe_formula = paste(na.omit(c("unit_sample + ref_onset_ref_event_time", discrete_covar_formula_input)), collapse = " + ")

      if(!(is.null(reg_weights))){
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | ", fe_formula," | 0 | cluster_on_this")),
                    data = ES_data, weights = ES_data[[reg_weights]],
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      } else{
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | ", fe_formula," | 0 | cluster_on_this")),
                    data = ES_data,
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      }

    }

    ES_data <- NULL
    gc()

    results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
    est <- NULL
    gc()
    results <- results[grep("att", rn)]
    results[grep("lead", rn), rn := gsub("lead", "-", rn)]
    results[, e := "Pooled"]
    results[, event_time := as.integer(gsub("att", "", rn))]
    results[, rn := "att"]
  }
  setnames(results, c("e"), onset_time_var)
  setnames(results, c("Estimate", "Cluster s.e."), c("estimate", "cluster_se"))
  setnames(results,"Pr(>|t|)","pval")
  results[, pval := round(pval,8)]

  if(homogeneous_ATT == FALSE & ipw == TRUE){
    flog.info("Estimated heterogeneous case with WLS using IPW.")
  } else if(homogeneous_ATT == FALSE & ipw == FALSE & is.null(reg_weights)){
    flog.info("Estimated heterogeneous case with OLS.")
  } else if(homogeneous_ATT == FALSE & ipw == FALSE & (!(is.null(reg_weights)))){
    flog.info("Estimated heterogeneous case with WLS.")
  } else if(homogeneous_ATT == TRUE & ipw == TRUE){
    flog.info("Estimated homogeneous case with WLS using IPW.")
  } else if(homogeneous_ATT == TRUE & ipw == FALSE & is.null(reg_weights)){
    flog.info("Estimated homogeneous case with OLS.")
  } else{
    flog.info("Estimated homogeneous case with WLS.")
  }

  return(results)
}


ES_estimate_std_did <- function(long_data,
                                outcomevar,
                                unit_var,
                                cal_time_var,
                                onset_time_var,
                                correct_pre_trends = FALSE,
                                omitted_event_time = -2,
                                cluster_vars = NULL,
                                std_subset_var = NULL,
                                std_subset_event_time = -1,
                                time_vary_confounds_low_dim = FALSE,
                                time_vary_confounds_high_dim = FALSE,
                                time_vary_confounds_cont = FALSE,
                                correct_time_vary_confounds = FALSE,
                                time_vary_covar) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  model_dt <- copy(long_data)

  if(!is.null(std_subset_var)){
    model_dt[, subset := max(as.integer(get(std_subset_var)*(event_time==std_subset_event_time))), by=c(unit_var)]
    model_dt <- model_dt[subset == 1]
    model_dt[, subset := NULL]
    gc()
  }

  start_cols <- copy(colnames(model_dt))

  # make event-time dummies by hand to supply them to main felm formula and extract coefficients by name
  # will omit the omitted_event_time (normalization) and the min event_time (otherwise collinear with tax_yr and tin fes)
  for (e in setdiff((min(model_dt$event_time) + 1):max(model_dt$event_time), omitted_event_time)) {
    var <- sprintf("event_time%s", e)
    model_dt[, (var) := as.integer(event_time == e)]
  }

  new_cols <- setdiff(colnames(model_dt), start_cols)
  event_time_cols <- gsub("\\-", "lead", new_cols) # "-" wreaks havoc otherwise, syntactically
  setnames(model_dt, new_cols, event_time_cols)

  # Should now be able to combine all of the above in a natural way
  event_time_form <- paste(event_time_cols, collapse = "+")
  fe_form <- paste(c(unit_var, cal_time_var), collapse = "+")

  if(correct_time_vary_confounds == TRUE & time_vary_confounds_low_dim == TRUE){
    fe_form <- paste(c(unit_var, cal_time_var, time_vary_covar), collapse = "+")
  }

  if( !is.null(cluster_vars) ){
    # cluster-robust SEs
    if( length(cluster_vars) > 1 ){
      cluster_on_this <- "cluster_on_this"
      model_dt[, (cluster_on_this) := .GRP, by = cluster_vars]
    } else if( length(cluster_vars) == 1){
      cluster_on_this <- cluster_vars
    }

    results_se_col_name <- "Cluster s.e."

  } else{
    # White / robust SEs
    cluster_on_this <- 0
    results_se_col_name <- "Robust s.e"
    # ^^ Note: not a typo -- when not clustering, summary(felm_object, robust = TRUE) returns with name "Robust s.e" missing period after "e"
  }

  if (correct_pre_trends == TRUE) {

    # Further omit min event_time + 1
    addl_event_time_to_omit <- abs(min(model_dt$event_time) + 1)
    event_time_form <- gsub(sprintf("event_timelead%s\\+", addl_event_time_to_omit), "", event_time_form)

    est <- felm(as.formula(sprintf("%s ~ %s + factor(%s)*%s | %s | 0 | %s", outcomevar, event_time_form, onset_time_var, cal_time_var, fe_form, cluster_on_this)),
                data = model_dt,
      nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
    )
    model_dt <- NULL
    gc()
  } else {

    est <- felm(as.formula(sprintf("%s ~ %s | %s | 0 | %s", outcomevar, event_time_form, fe_form, cluster_on_this)),
                data = model_dt,
      nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
    )
    model_dt <- NULL
    gc()
  }

  results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
  est <- NULL
  gc()

  results[grep("lead", rn), rn := gsub("lead", "-", rn)]
  results[, e := "Standard DiD"]
  results[, event_time := as.integer(gsub("event\\_time", "", rn))]
  results[, rn := "att"]

  setnames(results, c("Estimate", results_se_col_name), c("estimate", "cluster_se"))
  setnames(results, c("e"), onset_time_var)
  setnames(results,"Pr(>|t|)","pval")
  results[, pval := round(pval,8)]

  return(results)
}
