
#' @export
ES_simulate_data <- function(units = 1e4,
                             min_cal_time = 1999,
                             max_cal_time = 2005,
                             min_onset_time = 2001,
                             max_onset_time = 2006,
                             epsilon_sd = 0.2,
                             homogeneous_ATT = TRUE,
                             oversample_one_year = TRUE,
                             cohort_specific_trends = FALSE,
                             post_treat_dynamics = TRUE,
                             anticipation = FALSE,
                             cohort_specific_anticipation = FALSE) {

  Units <- 1:units
  Times <- min_cal_time:max_cal_time
  W <- min_onset_time:max_onset_time

  # Make balanced panel in calendar time
  sim_data <- data.table(tin = rep(Units, length(Times)))
  setorderv(sim_data, "tin")
  sim_data[, tax_yr := rep(Times, length(Units))]

  # Different cohorts determined by time of treatment
  win_yr_data <- data.table(tin = Units, win_yr = sample(W, size = length(Units), replace = TRUE))
  sim_data <- merge(sim_data, win_yr_data, by = "tin")
  win_yr_data <- NULL

  if (oversample_one_year == TRUE) {

    # $@$ Randomly make the 2003 cohort much larger
    sim_data[, recode := runif(1, min = 0, max = 1), list(tin)]
    sim_data[recode > 0.50, win_yr := 2003]
  }

  # Different initial condition by cohort (just to generate level differences)
  cohort_fes <- data.table(win_yr = sort(unique(sim_data$win_yr)), alpha_e = runif(uniqueN(sim_data$win_yr), min = 0.6, max = 0.9))

  # Generic (but parallel) trend
  tax_yr_fes <- data.table(tax_yr = sort(unique(sim_data$tax_yr)), delta_t = runif(uniqueN(sim_data$tax_yr), min = -0.1, max = 0.1))

  sim_data <- merge(sim_data, cohort_fes, by = "win_yr")
  sim_data <- merge(sim_data, tax_yr_fes, by = "tax_yr")

  if (homogeneous_ATT == FALSE) {

    # Different on-impact treatment effect by cohort
    cohort_specific_te_0 <- data.table(
      win_yr = sort(unique(sim_data$win_yr)),
      # te_0 = c(-0.2431251, -0.2135971, -0.2690795, -0.2816714, -0.2741780, -0.2002287)
      te_0 = (-2:3 - 3.5) / 10
    )

    # No anticipation; nonlinear cohort-specific dynamics (growing TE) after "on-impact"
    # will do a quadratic in event_time post-treatment with cohort-specific loadings
    cohort_specific_post_linear <- data.table(
      win_yr = sort(unique(sim_data$win_yr)),
      # trend_load = c(-0.01788569, -0.01061920, -0.01784453, -0.01432203, -0.01355335, -0.02377261)
      trend_load = (-2:3 - 3.5) / 100
    )
    cohort_specific_post_quad <- data.table(
      win_yr = sort(unique(sim_data$win_yr)),
      # quad_load = c(-0.0002752961, -0.0001757367, -0.0002231508, -0.0001444998, -0.0002591208, -0.0001080127)
      quad_load = (-2:3 - 3.5) / 10000
    )

    params <- merge(cohort_specific_te_0, cohort_specific_post_linear, by = "win_yr", all = TRUE)
    params <- merge(params, cohort_specific_post_quad, by = "win_yr", all = TRUE)
    params[, win_yr := as.character(win_yr)]

    sim_data <- merge(sim_data, cohort_specific_te_0, by = "win_yr")
    sim_data <- merge(sim_data, cohort_specific_post_linear, by = "win_yr")
    sim_data <- merge(sim_data, cohort_specific_post_quad, by = "win_yr")
  } else {

    # Homogeneous on-impact treatment effect by cohort
    te_0 <- -0.2431251

    # No anticipation; nonlinear cohort-specific dynamics (growing TE) after "on-impact"
    # will do a quadratic in event_time post-treatment with cohort-specific loadings
    trend_load <- -0.01784453
    quad_load <- -0.0001432203

    params <- c(te_0, trend_load, quad_load)
    names(params) <- c("te_0", "trend_load", "quad_load")
  }

  if (cohort_specific_trends == TRUE) {

    # Different calendar time trends by cohort
    cohort_specific_cal_time_trend <- data.table(
      win_yr = sort(unique(sim_data$win_yr)),
      cal_linear_trend = (1:6 - 3.5) / 100
    )
    sim_data <- merge(sim_data, cohort_specific_cal_time_trend, by = "win_yr")
    sim_data[, delta_t := delta_t + (tax_yr - min_cal_time) * cal_linear_trend]
  }


  setorderv(sim_data, c("tin", "tax_yr"))
  sim_data[, event_time := tax_yr - win_yr]

  # Have all the ingredients in place to determine the observed outcome, but for structural error
  sim_data[, epsilon := rnorm(dim(sim_data)[1], mean = 0, sd = epsilon_sd)]

  if (post_treat_dynamics == FALSE & homogeneous_ATT == TRUE) {
    trend_load <- 0
    quad_load <- 0

    params[names(params) %in% c("trend_load", "quad_load")] <- 0
  } else if (post_treat_dynamics == FALSE & homogeneous_ATT == FALSE) {
    sim_data[, trend_load := 0]
    sim_data[, quad_load := 0]

    params[, trend_load := 0]
    params[, quad_load := 0]
  }

  if(homogeneous_ATT == T & anticipation == T & cohort_specific_anticipation == F){
    antic_params <- c(-.1,-.1)
    params <- c(params,antic_params[1],antic_params[2])
  } else if(homogeneous_ATT == T & anticipation == T & cohort_specific_anticipation == T){

    # Different two-period anticipations effects by cohort
    cohort_specific_antic_params <- data.table(
      win_yr = sort(unique(sim_data$win_yr)),
      # te_0 = c(-0.2431251, -0.2135971, -0.2690795, -0.2816714, -0.2741780, -0.2002287)
      antic_1 = (-3:2) / 10,
      antic_2 = (-3:2) / 10
    )

    sim_data <- merge(sim_data, cohort_specific_antic_params, by = "win_yr")
    cohort_specific_antic_params[, win_yr := as.character(win_yr)]

  }

  sim_data[, outcome := alpha_e + delta_t + epsilon]
  sim_data[event_time >= 0, outcome := outcome + te_0 + (event_time * trend_load) + ((event_time)^2 * quad_load)]
  if(homogeneous_ATT == T & anticipation == T & cohort_specific_anticipation == F){
    sim_data[event_time < 0, outcome := outcome + (event_time==(-2))*antic_params[1] + (event_time==(-1))*antic_params[2]]
  } else if(homogeneous_ATT == T & anticipation == T & cohort_specific_anticipation == T){
    sim_data[event_time < 0, outcome := outcome + (event_time==(-2))*antic_1 + (event_time==(-1))*antic_2]
  }

  output <- list()
  output[[1]] <- sim_data
  output[[2]] <- params

  if(homogeneous_ATT == T & anticipation == T & cohort_specific_anticipation == T){
    output[[3]] <- cohort_specific_antic_params
  }


  return(output)
}

ES_simulate_estimator_comparison <- function(units = 1e4,
                                             seed = 1,
                                             oversample_one_year = FALSE,
                                             omitted_event_time = -2,
                                             cohort_specific_trends = FALSE,
                                             correct_pre_trends = FALSE,
                                             anticipation = FALSE,
                                             max_control_gap = Inf,
                                             min_control_gap = 1,
                                             homogeneous_ATT = TRUE,
                                             cohort_specific_anticipation = FALSE) {
  set.seed(seed)

  sim_result <- ES_simulate_data(units,
                                 oversample_one_year = oversample_one_year,
                                 cohort_specific_trends = cohort_specific_trends,
                                 anticipation = anticipation,
                                 homogeneous_ATT = homogeneous_ATT,
                                 cohort_specific_anticipation = cohort_specific_anticipation)

  if (correct_pre_trends == TRUE) {
    long_data <- ES_parallelize_trends(long_data = sim_result[[1]],outcomevar = "outcome",cal_time_var = "tax_yr",onset_time_var = "win_yr"
    )
  } else {
    long_data <- sim_result[[1]]
  }

  params <- sim_result[[2]]

  if(homogeneous_ATT == T & anticipation == T & cohort_specific_anticipation == T){
    cohort_specific_antic_params <- sim_result[[3]]
  }


  ES_data <- ES_clean_data(
    long_data = long_data,
    outcomevar = "outcome",
    unit_var = "tin",
    cal_time_var = "tax_yr",
    onset_time_var = "win_yr",
    min_control_gap = min_control_gap,
    max_control_gap = max_control_gap,
    omitted_event_time = omitted_event_time
  )

  mean2003_2002 = mean(ES_data[win_yr == 2002 & tax_yr == 2003]$outcome)
  mean1999_2002 = mean(ES_data[win_yr == 2002 & tax_yr == 1999]$outcome)
  mean2003_2006 = mean(ES_data[win_yr == 2006 & tax_yr == 2003]$outcome)
  mean1999_2006 = mean(ES_data[win_yr == 2006 & tax_yr == 1999]$outcome)

  event_time_1_coh_2002 = (mean2003_2002 - mean1999_2002) - (mean2003_2006 - mean1999_2006)

  ref_cohort_check = list()
  z = 0
  for(w in sort(unique(ES_data$ref_onset_time))){
    z = z + 1
    ref_cohort_check[[z]] = setorderv(ES_data[ref_onset_time == w, .N, by = c("win_yr", "ref_onset_time", "ref_event_time")], c("win_yr", "ref_onset_time", "ref_event_time"))
  }
  ref_cohort_check = rbindlist(ref_cohort_check, use.names = TRUE)

  # check how we are doing against true values

  ES_results_hetero <- ES_estimate_ATT(
    ES_data = ES_data,
    onset_time_var = "win_yr",
    cluster_vars = c("tin", "tax_yr"),
    homogeneous_ATT = FALSE,
    omitted_event_time = omitted_event_time
  )

  ES_results_homo <- ES_estimate_ATT(
    ES_data = ES_data,
    onset_time_var = "win_yr",
    cluster_vars = c("tin", "tax_yr"),
    homogeneous_ATT = TRUE,
    omitted_event_time = omitted_event_time
  )

  es_results <- ES_estimate_std_did(
    long_data = sim_result[[1]],
    outcomevar = "outcome",
    unit_var = "tin",
    cal_time_var = "tax_yr",
    onset_time_var = "win_yr",
    cohort_specific_trends = correct_pre_trends,
    omitted_event_time = omitted_event_time
  )

  figdata <- rbindlist(list(ES_results_hetero, ES_results_homo, es_results), use.names = TRUE)

  if(homogeneous_ATT==T & cohort_specific_anticipation==F){

    # getting true value in
    figdata_temp <- copy(figdata)
    figdata_temp[, te_0 := params[1]]
    figdata_temp[, trend_load := params[2]]
    figdata_temp[, quad_load := params[3]]
    figdata_temp[event_time < 0, true_value := 0]
    if(anticipation){
      figdata_temp[event_time < 0, true_value := (event_time==(-2))*params[4] + (event_time==(-1))*params[5]]
    }
    figdata_temp[event_time >= 0, true_value := te_0 + (event_time * trend_load) + ((event_time)^2 * quad_load)]

    figdata_temp <- figdata_temp[win_yr == "Standard DiD"]
    figdata_temp <- figdata_temp[, list(rn, win_yr, event_time, true_value)]
    figdata_temp[win_yr == "Standard DiD", win_yr := "True Value"]
    setnames(figdata_temp, c("true_value"), c("estimate"))

    figdata <- rbindlist(list(figdata, figdata_temp), use.names = TRUE, fill = TRUE)
    figdata[is.na(cluster_se), cluster_se := 0]

    figdata[, jitter := .GRP, by = win_yr]

    fig <- ggplot(aes(
      x = event_time + (jitter - 3.5) / 14,
      y = estimate, colour = factor(win_yr)
    ), data = figdata) + geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - 1.96 * cluster_se, ymax = estimate + 1.96 * cluster_se)) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time", y = "ATT", color = "Estimator") +
      annotate(geom = "text", x = omitted_event_time, y = 0, label = "(Omitted)", size = 4)

    event_time_1_coh_2002 = c(event_time_1_coh_2002, figdata[win_yr == "True Value" & event_time == 1]$estimate)
    names(event_time_1_coh_2002) <- c("by_hand_DiD", "true_value")

    results = list()
    results[[1]] = fig
    results[[2]] = ref_cohort_check
    results[[3]] = event_time_1_coh_2002

  } else if(homogeneous_ATT==T & anticipation==T & cohort_specific_anticipation==T){

    # getting true value in
    figdata_temp <- copy(figdata)
    figdata_temp[, te_0 := params[1]]
    figdata_temp[, trend_load := params[2]]
    figdata_temp[, quad_load := params[3]]
    figdata_temp[event_time < 0, true_value := 0]
    if(anticipation == T & cohort_specific_anticipation == F){
      figdata_temp[event_time < 0, true_value := (event_time==(-2))*params[4] + (event_time==(-1))*params[5]]
    }
    figdata_temp[event_time >= 0, true_value := te_0 + (event_time * trend_load) + ((event_time)^2 * quad_load)]

    figdata_temp = merge(figdata_temp, cohort_specific_antic_params, by = "win_yr", all = TRUE)
    figdata_temp[between(event_time, omitted_event_time, 0, incbounds = FALSE), true_value := (event_time==(-2))*antic_1 + (event_time==(-1))*antic_2]
    figdata_temp[event_time <= omitted_event_time, true_value := 0]

    figdata_temp <- figdata_temp[win_yr == "Standard DiD"]
    figdata_temp <- figdata_temp[, list(rn, win_yr, event_time, true_value)]
    figdata_temp[win_yr == "Standard DiD", win_yr := "True Value"]
    setnames(figdata_temp, c("true_value"), c("estimate"))

    figdata <- rbindlist(list(figdata, figdata_temp), use.names = TRUE, fill = TRUE)
    figdata[is.na(cluster_se), cluster_se := 0]

    figdata = merge(figdata, cohort_specific_antic_params, by = "win_yr", all = TRUE)
    figdata[between(event_time, omitted_event_time, 0, incbounds = FALSE), true_value := (event_time==(-2))*antic_1 + (event_time==(-1))*antic_2]
    figdata[event_time <= omitted_event_time & !(win_yr %in% c("Pooled", "Standard DiD", "True Value")), true_value := NA]

    figdata[, jitter := .GRP, by = win_yr]

    fig <- ggplot(aes(
      x = event_time + (jitter - 3.5) / 14,
      y = estimate, colour = factor(win_yr)
    ), data = figdata) + geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - 1.96 * cluster_se, ymax = estimate + 1.96 * cluster_se)) +
      geom_point(aes(y = true_value), color = "#000000", shape = 4)+
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time", y = "ATT", color = "Win Year") +
      annotate(geom = "text", x = omitted_event_time, y = 0, label = "(Omitted)", size = 4)

    event_time_1_coh_2002 = c(event_time_1_coh_2002, figdata[win_yr == "True Value" & event_time == 1]$estimate)
    names(event_time_1_coh_2002) <- c("by_hand_DiD", "true_value")

    results = list()
    results[[1]] = fig
    results[[2]] = ref_cohort_check
    results[[3]] = event_time_1_coh_2002

  } else{

    # getting true value in
    figdata_temp <- copy(figdata)
    figdata_temp = merge(figdata_temp, params, by = "win_yr", all = TRUE)
    figdata_temp[event_time < 0 & !(win_yr %in% c("Pooled", "Standard DiD")), true_value := 0]
    figdata_temp[event_time >= 0, true_value := te_0 + (event_time * trend_load) + ((event_time)^2 * quad_load)]

    figdata_temp[, jitter := .GRP, by = win_yr]

    fig <- ggplot(aes(
      x = event_time + (jitter - 3.5) / 14,
      y = estimate, colour = factor(win_yr)
    ), data = figdata_temp) + geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - 1.96 * cluster_se, ymax = estimate + 1.96 * cluster_se)) +
      geom_point(aes(y = true_value), color = "#000000", shape = 4)+
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time", y = "ATT", color = "Estimator") +
      annotate(geom = "text", x = omitted_event_time, y = 0, label = "(Omitted)", size = 4)

    results = list()
    results[[1]] = fig
    results[[2]] = ref_cohort_check

  }

  return(results)
}
