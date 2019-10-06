#' @export
# Calculate propensity score weights within each DiD sample
# ==> two time periods for reference cohort, and two time periods for all later-treated units
ES_make_ipw_dt <- function(did_dt,
                           unit_var,
                           cal_time_var,
                           discrete_covars = NULL,
                           cont_covars = NULL,
                           ref_discrete_covars = NULL,
                           ref_cont_covars = NULL,
                           omitted_event_time = -2,
                           reg_weights = NULL,
                           ref_reg_weights = NULL,
                           ipw_model = "linear",
                           ipw_composition_change = FALSE){

  if(ipw_composition_change == FALSE){

    # Construct 'pre' version of each covariate, which will be time-invariant in did_dt

    if((!is.null(discrete_covars)) | (!is.null(ref_discrete_covars)) ){

      pre_discrete_covars = c()

      for(d in unique(na.omit(c(discrete_covars, ref_discrete_covars)))){
        d_name <- sprintf("%s_pre", d)
        did_dt[, has_d_pre := max(as.integer((!is.na(get(d))) * (ref_event_time == omitted_event_time))), by = list(get(unit_var))]
        did_dt <- did_dt[has_d_pre == 1] # these will be the only obs used anyway when estimating the propensity score
        # use sum() instead of max() below in case 'd' takes on negative values (as max() would then return 0)
        did_dt[, (d_name) := sum(get(d) * (ref_event_time == omitted_event_time)), by = list(get(unit_var))]
        did_dt[, has_d_pre := NULL]
        gc()
        pre_discrete_covars = c(pre_discrete_covars, d_name)
      }

      pre_discrete_covar_formula_input = paste0("factor(", na.omit(pre_discrete_covars), ")", collapse = "+")
      rm(d_name)
    } else{
      pre_discrete_covars = NA
      pre_discrete_covar_formula_input = NA
    }

    if((!is.null(cont_covars)) | (!is.null(ref_cont_covars))){

      pre_cont_covars = c()

      for(c in unique(na.omit(c(cont_covars, ref_cont_covars)))){
        c_name <- sprintf("%s_pre", c)
        did_dt[, has_c_pre := max(as.integer((!is.na(get(c))) * (ref_event_time == omitted_event_time))), by = list(get(unit_var))]
        did_dt <- did_dt[has_c_pre == 1]  # these will be the only obs used anyway when estimating the propensity score
        # use sum() instead of max() below in case 'c' takes on negative values (as max() would then return 0)
        did_dt[, (c_name) := sum(get(c) * (ref_event_time == omitted_event_time)), by = list(get(unit_var))]
        did_dt[, has_c_pre := NULL]
        gc()
        pre_cont_covars = c(pre_cont_covars, c_name)
      }

      pre_cont_covar_formula_input = paste0(na.omit(pre_cont_covars), collapse = "+")
      rm(c_name)
    } else{
      pre_cont_covars = NA
      pre_cont_covar_formula_input = NA
    }

    if(is.na(pre_discrete_covar_formula_input)){
      formula_input <- pre_cont_covar_formula_input
    } else if(is.na(pre_cont_covar_formula_input)){
      formula_input <- pre_discrete_covar_formula_input
    } else{
      formula_input <- paste(pre_discrete_covar_formula_input,
                             pre_cont_covar_formula_input,
                             sep = "+"
      )
    }

    # Now impose common support
    # Can revisit -- for now take a simplified approach
    # Form strata using the discrete covariates (if present)
    # Rely solely on parametric model for continuous covariates

    if((!is.null(discrete_covars)) | (!is.null(ref_discrete_covars)) ){
      did_dt[, strata := .GRP, keyby = c(pre_discrete_covars)]
      treat_strata <- did_dt[treated == 1, sort(unique(strata))]
      control_strata <- did_dt[treated == 0, sort(unique(strata))]
      common_strata <- intersect(treat_strata, control_strata)
      did_dt <- did_dt[strata %in% common_strata]
      gc()
    }

    if(ipw_model == "linear"){
      if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
        did_dt[, pr := predict(lm(as.formula(paste0("treated ~ ", formula_input)),data = did_dt, weights = did_dt[[na.omit(unique(c(reg_weights, ref_reg_weights)))]]), type = "response")]
      } else{
        did_dt[, pr := predict(lm(as.formula(paste0("treated ~ ", formula_input)),data = did_dt), type = "response")]
      }
    } else if(ipw_model %in% c("logit", "probit")){
      if(!(is.null(reg_weights)) | !(is.null(ref_reg_weights))){
        did_dt[, pr := predict(glm(as.formula(paste0("treated ~ ", formula_input)),family = quasibinomial(link = ipw_model), data = did_dt, weights = did_dt[[na.omit(unique(c(reg_weights, ref_reg_weights)))]]), type = "response")]
      } else{
        did_dt[, pr := predict(glm(as.formula(paste0("treated ~ ", formula_input)),family = quasibinomial(link = ipw_model), data = did_dt), type = "response")]
      }
    }

    rm(pre_discrete_covar_formula_input, pre_cont_covar_formula_input, formula_input)
    new_cols <- setdiff(colnames(did_dt), c(unit_var, cal_time_var, "pr"))
    did_dt[, (new_cols) := NULL]
    gc()
  }

  return(did_dt)
}
