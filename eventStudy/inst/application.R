rm(list = ls(all = TRUE))
gc()

options(scipen = 10000)
options(show.error.locations=TRUE)
options(warn=1)


library(data.table)
library(testthat)
library(futile.logger)
library(lfe)
library(multiwayvcov)
library(lmtest)
library(ggplot2)
library(scales)
library(eventStudy)




main_figs <- function(){

  for(s in c(1,2)){
    for(u in c(1e3, 1e4)){
      fig = ES_simulate_estimator_comparison(units = u, seed = s, homogeneous_ATT = TRUE)[[1]]
      ggsave(sprintf("inst/figures/event_time_all_seed%s_size%s_homogeneous.png", s, u), width = 12, height = 5)
    }
  }

  for(s in c(1,2)){
    for(u in c(1e3, 1e4)){
      fig = ES_simulate_estimator_comparison(units = u, seed = s, homogeneous_ATT = FALSE)[[1]]
      ggsave(sprintf("inst/figures/event_time_all_seed%s_size%s_heterogeneous.png", s, u), width = 12, height = 5)
    }
  }


  u = 1e4
  s = 1

  fig = ES_simulate_estimator_comparison(units = u, seed = s, cohort_specific_trends=T, correct_pre_trends = F)[[1]]
  ggsave(sprintf("inst/figures/event_time_all_seed%s_size%s_trends.png", s, u), width = 12, height = 5)
  fig = ES_simulate_estimator_comparison(units = u, seed = s, cohort_specific_trends=T, correct_pre_trends = T)[[1]]
  ggsave(sprintf("inst/figures/event_time_all_seed%s_size%s_trends_corrected.png", s, u), width = 12, height = 5)

}


antic_figs <- function(){

  u = 1e4
  s = 1

  results_anticipation = ES_simulate_estimator_comparison(units = u, seed = s, anticipation=T, omitted_event_time = -3, cohort_specific_anticipation = T)
  fig = results_anticipation[[1]]
    ggsave(sprintf("inst/figures/event_time_all_seed%s_size%s_anticipation.png", s, u), width = 12, height = 5)
  sample_check_anticipation = results_anticipation[[2]]
  sample_check_anticipation[, anticipation_corrected := 0]
  event_time_1_coh_2002 = results_anticipation[[3]]

  results_anticipation_corrected = ES_simulate_estimator_comparison(units = u, seed = s, oversample_one_year=F, anticipation=T, min_control_gap = 3, max_control_gap = Inf,omitted_event_time = -3, cohort_specific_anticipation = T)
  fig = results_anticipation_corrected[[1]]
    ggsave(sprintf("inst/figures/event_time_all_seed%s_size%s_anticipation_corrected.png", s, u), width = 12, height = 5)
  sample_check_anticipation_corrected = results_anticipation_corrected[[2]]
  sample_check_anticipation_corrected[, anticipation_corrected := 1]
  event_time_1_coh_2002_corrected = results_anticipation_corrected[[3]]

  # view_cohorttrends <- ES_simulate_data(units=1e4, cohort_specific_trends=T)[[1]][,list(outcome=mean(outcome)),list(tax_yr,win_yr)]
  # ggplot(aes(x=tax_yr,y=outcome,color=factor(win_yr)),data=view_cohorttrends) + geom_line() + theme_bw(base_size=16) + labs(x="Year",y="Outcome",colour="Cohort")
  #
  # view_anticipation <- ES_simulate_data(units=1e4, anticipation = T)[[1]][,list(outcome=mean(outcome)),list(tax_yr,win_yr)]
  # ggplot(aes(x=tax_yr,y=outcome,color=factor(win_yr)),data=view_anticipation) + geom_line() + theme_bw(base_size=16) + labs(x="Year",y="Outcome",colour="Cohort")

  results = list()
  results[[1]] = rbindlist(list(sample_check_anticipation, sample_check_anticipation_corrected), use.names = TRUE)
  results[[2]] = c(event_time_1_coh_2002, event_time_1_coh_2002_corrected)

  return(results)

}




test_restrictions <- function(){

  units = 1e4
  seed = 1
  oversample_one_year = FALSE
  omitted_event_time = -2
  cohort_specific_trends = FALSE
  correct_pre_trends = FALSE
  anticipation = FALSE
  max_control_gap = 1
  min_control_gap = 1

  set.seed(seed)

  sim_result <- ES_simulate_data(units,
                                 oversample_one_year = oversample_one_year,
                                 cohort_specific_trends = cohort_specific_trends,
                                 anticipation = anticipation)

  if (correct_pre_trends == TRUE) {
    long_data <- ES_parallelize_trends(long_data = sim_result[[1]],outcomevar = "outcome",cal_time_var = "tax_yr",onset_time_var = "win_yr"
    )
  } else {
    long_data <- sim_result[[1]]
  }

  params <- sim_result[[2]]


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

  # check how we are doing against true values

  ES_results_hetero <- ES_estimate_ATT(
    ES_data = ES_data,
    onset_time_var = "win_yr",
    cluster_vars = c("tin", "tax_yr"),
    homogeneous_ATT = FALSE
  )

  ES_results_homo <- ES_estimate_ATT(
    ES_data = ES_data,
    onset_time_var = "win_yr",
    cluster_vars = c("tin", "tax_yr"),
    homogeneous_ATT = TRUE
  )

  check1 <- ES_results_hetero[,list(estimate_heterog=mean(estimate,na.rm=T)),list(event_time)][order(event_time)][,.(event_time,estimate_heterog)]
  check2 <- ES_results_homo[,list(event_time,estimate_homog=estimate)]
  check <- merge(check1,check2,by='event_time')
  print(check[])

}



