

#' @export
ES_estimate_ATT <- function(ES_data, outcomevar, onset_time_var, cluster_vars, homogeneous_ATT = TRUE, omitted_event_time = -2) {

  # Just in case, we immediately make a copy of the input ES_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  model_dt = copy(ES_data)

  onset_times <- sort(unique(model_dt[, .N, by = eval(onset_time_var)][[onset_time_var]]))
  min_onset_time <- min(onset_times)
  max_onset_time <- max(onset_times)

  model_dt[, ref_onset_ref_event_time := .GRP, by = list(ref_onset_time, ref_event_time)]
  model_dt[, unit_sample := .GRP, by = list(get(onset_time_var), catt_specific_sample, ref_onset_time)]
  model_dt[, cluster_on_this := .GRP, by = cluster_vars]

  start_cols <- copy(colnames(model_dt))

  if (homogeneous_ATT == FALSE) {
    # for (h in min(model_dt$ref_onset_time):max(model_dt$ref_onset_time)) {
    #   var <- sprintf("ref_onset_time%s", h)
    #   model_dt[, (var) := as.integer(ref_onset_time == h)]
    #   var <- sprintf("ref_onset_time%s_treated", h)
    #   model_dt[, (var) := as.integer(ref_onset_time == h & get(onset_time_var) == h)]
    # }

    for (h in min(model_dt$ref_onset_time):max(model_dt$ref_onset_time)) {
      for (r in setdiff(min(model_dt[ref_onset_time == h]$ref_event_time):max(model_dt[ref_onset_time == h]$ref_event_time), omitted_event_time)) {
        # var <- sprintf("ref_onset_time%s_et%s", h, r)
        # model_dt[, (var) := as.integer(ref_onset_time == h & ref_event_time == r)]
        var <- sprintf("ref_onset_time%s_catt%s", h, r)
        model_dt[, (var) := as.integer(ref_onset_time == h & ref_event_time == r & get(onset_time_var) == h)]
      }
    }

    new_cols <- setdiff(colnames(model_dt), start_cols)

    # new_cols_used <- setdiff(
    #   new_cols,
    #   c(sprintf("ref_onset_time%s", min_onset_time))
    # )

    new_cols_used <- gsub("\\-", "lead", new_cols) # "-" wreaks havoc otherwise, syntactically
    setnames(model_dt, new_cols, new_cols_used)

    # Should now be able to combine all of the above in a natural way
    felm_formula_input <- paste(new_cols_used, collapse = "+")

    est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | cluster_on_this")),
      data = model_dt,
      nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
    )
    gc()

    all_added_cols <- setdiff(colnames(model_dt), start_cols)

    model_dt[, (all_added_cols) := NULL]
    gc()

    results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
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

    setnames(results, c("e"), onset_time_var)
    setnames(results, c("Estimate", "Cluster s.e."), c("estimate", "cluster_se"))

    results <- results[rn == "catt"]
    gc()
  } else {
    # for (h in min(model_dt$ref_onset_time):max(model_dt$ref_onset_time)) {
    #   for (r in setdiff(min(model_dt[ref_onset_time == h]$ref_event_time):max(model_dt[ref_onset_time == h]$ref_event_time), omitted_event_time)) {
    #     var <- sprintf("ref_onset_time%s_et%s", h, r)
    #     model_dt[, (var) := as.integer(ref_onset_time == h & ref_event_time == r)]
    #   }
    # }

    for (r in setdiff(min(model_dt$ref_event_time):max(model_dt$ref_event_time), omitted_event_time)) {
      # var = sprintf("et%s", r)
      # model_dt[, (var) := as.integer(ref_event_time == r)]
      var <- sprintf("att%s", r)
      model_dt[, (var) := as.integer(ref_event_time == r & ref_onset_time == get(onset_time_var))]
    }

    new_cols <- setdiff(colnames(model_dt), start_cols)

    new_cols_used <- gsub("\\-", "lead", new_cols) # "-" wreaks havoc otherwise, syntactically
    setnames(model_dt, new_cols, new_cols_used)

    # Should now be able to combine all of the above in a natural way
    felm_formula_input <- paste(c(new_cols_used), collapse = "+")

    est <- felm(as.formula(paste0(eval(outcomevar), " ~ ", felm_formula_input, " | unit_sample + ref_onset_ref_event_time | 0 | cluster_on_this")),
      data = model_dt,
      nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
    )
    gc()

    all_added_cols <- setdiff(colnames(model_dt), start_cols)

    model_dt[, (all_added_cols) := NULL]
    gc()

    results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
    gc()

    results <- results[grep("att", rn)]
    results[grep("lead", rn), rn := gsub("lead", "-", rn)]
    results[, e := "Pooled"]
    results[, event_time := as.integer(gsub("att", "", rn))]
    results[, rn := "att"]

    setnames(results, c("e"), onset_time_var)
    setnames(results, c("Estimate", "Cluster s.e."), c("estimate", "cluster_se"))
  }

  setnames(results,"Pr(>|t|)","pval")
  results[, pval := round(pval,8)]

  return(results)
}


ES_estimate_std_did <- function(long_data,
                                outcomevar,
                                unit_var,
                                cal_time_var,
                                onset_time_var,
                                cohort_specific_trends = FALSE,
                                omitted_event_time = -2,
                                cluster_vars = NULL) {

  # Just in case, we immediately make a copy of the input long_data and run everything on the full copy
  # Can revisit this to remove the copy for memory efficiency at a later point.

  model_dt <- copy(long_data)
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

  if (cohort_specific_trends == TRUE) {

    # Further omit min event_time + 1
    addl_event_time_to_omit <- abs(min(model_dt$event_time) + 1)
    event_time_form <- gsub(sprintf("event_timelead%s\\+", addl_event_time_to_omit), "", event_time_form)

    est <- felm(as.formula(sprintf("%s ~ %s + factor(%s)*%s | %s | 0 | %s", outcomevar, event_time_form, onset_time_var, cal_time_var, fe_form, cluster_on_this)),
                data = model_dt,
      nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
    )
    gc()
  } else {
    est <- felm(as.formula(sprintf("%s ~ %s | %s | 0 | %s", outcomevar, event_time_form, fe_form, cluster_on_this)),
                data = model_dt,
      nostats = FALSE, keepX = FALSE, keepCX = FALSE, psdef = FALSE, kclass = FALSE
    )
    gc()
  }


  all_added_cols <- setdiff(colnames(model_dt), start_cols)

  model_dt[, (all_added_cols) := NULL]
  gc()

  results <- as.data.table(summary(est, robust = TRUE)$coefficients, keep.rownames = TRUE)
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
