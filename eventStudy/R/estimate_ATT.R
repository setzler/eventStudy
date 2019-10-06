#' @export
ES_estimate_ATT <- function(ES_data,
                            outcomevar,
                            unit_var,
                            onset_time_var,
                            cluster_vars,
                            residualize_covariates = FALSE,
                            discrete_covars = NULL,
                            cont_covars = NULL,
                            ref_discrete_covars = NULL,
                            ref_cont_covars = NULL,
                            homogeneous_ATT = TRUE,
                            omitted_event_time = -2,
                            reg_weights = NULL,
                            ref_reg_weights = NULL,
                            ipw = FALSE,
                            ipw_composition_change = FALSE,
                            add_unit_fes = FALSE,
                            linearDiD = FALSE,
                            linearDiD_treat_var = NULL) {

  onset_times <- sort(unique(ES_data[, .N, by = eval(onset_time_var)][[onset_time_var]]))
  min_onset_time <- min(onset_times)
  max_onset_time <- max(onset_times)

  ES_data[, sortorder := .I]
  ES_data[, ref_onset_ref_event_time := .GRP, keyby = list(ref_onset_time, ref_event_time)]
  setorderv(ES_data, "sortorder")

  if(add_unit_fes == TRUE){
    ES_data[, unit_sample := .GRP, keyby = c(unit_var, "catt_specific_sample", "ref_onset_time")]
  } else{
    ES_data[, unit_sample := .GRP, keyby = c(onset_time_var, "catt_specific_sample", "ref_onset_time")]
  }
  setorderv(ES_data, "sortorder")

  # felm() will necessarily drop singletons on unit_sample, so we just exclude those now
  # in practice, this step should only have bite if add_unit_fes == TRUE
  ES_data[, unit_var_block_count := .N, by = list(unit_sample)]
  count_dropcases <- dim(ES_data[unit_var_block_count <= 1])[1]
  if(count_dropcases > 0){
    ES_data <- ES_data[unit_var_block_count > 1]
    rm(count_dropcases)
  }
  ES_data[, unit_var_block_count := NULL]
  gc()

  # define the unique reference onset times now, after all known pre-drops of observations have been done
  ref_onset_times <- sort(unique(ES_data[, .N, by = list(ref_onset_time)][["ref_onset_time"]]))

  if(!(is.null(cluster_vars))){
    ES_data[, cluster_on_this := .GRP, keyby = cluster_vars]
    cluster_on_this <- "cluster_on_this"
  } else{
    # White / robust SEs
    cluster_on_this <- 0
  }
  setorderv(ES_data, "sortorder")
  ES_data[, sortorder := NULL]

  if(residualize_covariates == FALSE & ipw == FALSE){

    if((!is.null(discrete_covars)) | (!is.null(ref_discrete_covars)) ){

      # For felm(), want to tell which factor combination will be sent to the intercept
      # Let's set it to the (approximate) median

      i <- 0
      reference_lookup <- list()
      discrete_covar_formula_input <- c()
      for(var in unique(na.omit(c(discrete_covars, ref_discrete_covars)))){
        i <- i + 1
        levels <- sort(ES_data[, unique(get(var))])
        median_pre_value <- levels[round(length(levels)/2)] # approximate median
        reference_lookup[[i]] <- data.table(varname = var, reference_level = as.character(median_pre_value))

        discrete_covar_formula_input <- c(discrete_covar_formula_input,
                                          sprintf("relevel(as.factor(%s), ref = '%s')", var, median_pre_value))

        var <- NULL
        median_pre_value <- NULL
      }
      reference_lookup <- rbindlist(reference_lookup, use.names = TRUE)

    } else{
      discrete_covar_formula_input = NA
    }

    if((!is.null(cont_covars)) | (!is.null(ref_cont_covars))){
      cont_covar_formula_input = paste0(unique(na.omit(c(cont_covars, ref_cont_covars))), collapse = "+")
    } else{
      cont_covar_formula_input = NA
    }

  }

  start_cols <- copy(colnames(ES_data))

  if (homogeneous_ATT == FALSE) {

    for (h in intersect((min(ES_data$ref_onset_time):max(ES_data$ref_onset_time)), ref_onset_times)) {
      for (r in setdiff(min(ES_data[ref_onset_time == h]$ref_event_time):max(ES_data[ref_onset_time == h]$ref_event_time), omitted_event_time)) {
        var <- sprintf("ref_onset_time%s_catt%s", h, r)
        ES_data[, (var) := as.integer(ref_onset_time == h & ref_event_time == r & get(onset_time_var) == h)]
        if(linearDiD == TRUE){
          ES_data[, (var) := as.numeric(get(var) * get(linearDiD_treat_var))]
        }
      }
    }

    new_cols <- setdiff(colnames(ES_data), start_cols)
    new_cols_used <- gsub("\\-", "lead", new_cols) # "-" wreaks havoc otherwise, syntactically
    setnames(ES_data, new_cols, new_cols_used)

    # Should now be able to combine all of the above in a natural way
    felm_formula_input <- paste(new_cols_used, collapse = "+")

    if(ipw == TRUE){

      if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
        ES_data[, pw := pw * get(na.omit(unique(c(reg_weights, ref_reg_weights))))]
      }

      est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | ",  cluster_on_this)),
                  data = ES_data, weights = ES_data$pw,
                  nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
      )

    } else if(residualize_covariates == TRUE){

      if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | ",  cluster_on_this)),
                    data = ES_data, weights = ES_data[[na.omit(unique(c(reg_weights, ref_reg_weights)))]],
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      } else{
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | ",  cluster_on_this)),
                    data = ES_data,
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      }

    } else{

      felm_formula_input = paste(na.omit(c(felm_formula_input, cont_covar_formula_input)), collapse = " + ")
      fe_formula = paste(na.omit(c("unit_sample + ref_onset_ref_event_time", discrete_covar_formula_input)), collapse = " + ")

      if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | ", fe_formula," | 0 | ",  cluster_on_this)),
                    data = ES_data, weights = ES_data[[na.omit(unique(c(reg_weights, ref_reg_weights)))]],
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      } else{
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | ", fe_formula," | 0 | ",  cluster_on_this)),
                    data = ES_data,
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      }
    }

    # Grabbing the estimates and vcov to produce weighted avg SEs
    # Note: vcov already reflects weights used in estimation
    catt_coefs <- coef(est)

    if(!(is.null(cluster_vars))){
      catt_vcov <- est$clustervcv
    } else{
      catt_vcov <- est$robustvcv
    }

    if(!(is.null(cluster_vars))){
      results <- data.table(rn = rownames(est$beta),
                            estimate = as.vector(est$beta),
                            cluster_se = est$cse,
                            t = est$ctval,
                            pval = est$cpval,
                            reg_sample_size = est$N)
    } else{
      results <- data.table(rn = rownames(est$beta),
                            estimate = as.vector(est$beta),
                            cluster_se = est$rse,
                            t = est$rtval,
                            pval = est$rpval,
                            reg_sample_size = est$N)
    }
    est <- NULL
    gc()

    results[!grepl("ref\\_onset\\_time", rn), e := min_onset_time]
    results[, rn := gsub("lead", "-", rn)]
    for (c in min_onset_time:max_onset_time) {
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
      if(linearDiD == TRUE){
        ES_data[, (var) := as.numeric(get(var) * get(linearDiD_treat_var))]
      }
    }

    new_cols <- setdiff(colnames(ES_data), start_cols)
    new_cols_used <- gsub("\\-", "lead", new_cols) # "-" wreaks havoc otherwise, syntactically
    setnames(ES_data, new_cols, new_cols_used)

    # Should now be able to combine all of the above in a natural way
    felm_formula_input <- paste(c(new_cols_used), collapse = "+")

    if(ipw == TRUE){

      if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
        ES_data[, pw := pw * get(na.omit(unique(c(reg_weights, ref_reg_weights))))]
      }

      est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | ",  cluster_on_this)),
                  data = ES_data, weights = ES_data$pw,
                  nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
      )

    } else if(residualize_covariates == TRUE){

      if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | ",  cluster_on_this)),
                    data = ES_data, weights = ES_data[[na.omit(unique(c(reg_weights, ref_reg_weights)))]],
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      } else{
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | ",  cluster_on_this)),
                    data = ES_data,
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      }

    } else{

      felm_formula_input = paste(na.omit(c(felm_formula_input, cont_covar_formula_input)), collapse = " + ")
      fe_formula = paste(na.omit(c("unit_sample + ref_onset_ref_event_time", discrete_covar_formula_input)), collapse = " + ")

      if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | ", fe_formula," | 0 | ",  cluster_on_this)),
                    data = ES_data, weights = ES_data[[na.omit(unique(c(reg_weights, ref_reg_weights)))]],
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      } else{
        est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | ", fe_formula," | 0 | ",  cluster_on_this)),
                    data = ES_data,
                    nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
        )
      }

    }

    # Don't need to grab the estimates and vcov for this homogeneous case, so just do NAs for symmetry
    catt_coefs <- NA
    catt_vcov <- NA

    if(!(is.null(cluster_vars))){
      results <- data.table(rn = rownames(est$beta),
                            estimate = as.vector(est$beta),
                            cluster_se = est$cse,
                            t = est$ctval,
                            pval = est$cpval,
                            reg_sample_size = est$N)
    } else{
      results <- data.table(rn = rownames(est$beta),
                            estimate = as.vector(est$beta),
                            cluster_se = est$rse,
                            t = est$rtval,
                            pval = est$rpval,
                            reg_sample_size = est$N)
    }
    est <- NULL
    gc()

    results <- results[grep("att", rn)]
    results[grep("lead", rn), rn := gsub("lead", "-", rn)]
    results[, e := "Pooled"]
    results[, event_time := as.integer(gsub("att", "", rn))]
    results[, rn := "att"]
  }
  setnames(results, c("e"), onset_time_var)
  results[, pval := round(pval,8)]

  if(homogeneous_ATT == FALSE & ipw == TRUE){
    flog.info("Estimated heterogeneous case with WLS using IPW.")
  } else if(homogeneous_ATT == FALSE & ipw == FALSE & is.null(reg_weights) & is.null(ref_reg_weights)){
    flog.info("Estimated heterogeneous case with OLS.")
  } else if(homogeneous_ATT == FALSE & ipw == FALSE & (!(is.null(reg_weights)) | !(is.null(ref_reg_weights)))){
    flog.info("Estimated heterogeneous case with WLS.")
  } else if(homogeneous_ATT == TRUE & ipw == TRUE){
    flog.info("Estimated homogeneous case with WLS using IPW.")
  } else if(homogeneous_ATT == TRUE & ipw == FALSE & is.null(reg_weights) & is.null(ref_reg_weights)){
    flog.info("Estimated homogeneous case with OLS.")
  } else{
    flog.info("Estimated homogeneous case with WLS.")
  }

  output = list()
  output[[1]] <- results
  output[[2]] <- catt_coefs
  output[[3]] <- catt_vcov

  return(output)
}
