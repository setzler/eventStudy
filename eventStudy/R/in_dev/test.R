# Things to check:

rm(list = ls(all = TRUE))
gc()

options(scipen = 999)
options(show.error.locations=TRUE)
options(warn=1)

library(data.table)
library(futile.logger)
library(checkmate)
library(lfe)
library(ggplot2)
library(scales)
library(parallel)
library(eventStudy)

my_dt <- readRDS("~/Downloads/old_downloads/ES_testing/testdt.rds")
my_dt[, keep_case := between(age, 21, 64, incbounds = TRUE)]

# ###
# # TESTING THAT vanilla event study gives same estimates for ES and by_cohort_ES
# ###

collapse_table <- data.table(a = c("pre","on_impact", "post_avg"),b = c(list(-7:-1), list(1), list(1:5)))

# vanilla run
set.seed(1)
res0a <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0,
  calculate_collapse_estimates = TRUE,
  collapse_inputs = copy(collapse_table))

res0b <- by_cohort_ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0,
  calculate_collapse_estimates = TRUE,
  collapse_inputs = copy(collapse_table), cohor)
rm(collapse_table)

# ###
# # TESTING THAT new ES_check_inputs() works correctly
# ###

collapse_table <- data.table(a = c("pre","on_impact", "post_avg"),b = c(list(-7:-1), list(1), list(1:5)))

# vanilla run
set.seed(1)
res1a1 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0,
  discrete_covars = "age",
  calculate_collapse_estimates = TRUE,
  collapse_inputs = copy(collapse_table))
rm(collapse_table)

# no keep_case2, so should break there
set.seed(1)
res1a1 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case2",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case2",
  control_subset_event_time = 0)

# haven't defined any collapsed estimate inputs
set.seed(1)
res1a1 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0,
  calculate_collapse_estimates = TRUE)


# ###
# # TESTING THAT ref_discrete_covars is working correctly
# ###
#
# # should break as 'age2' isn't in the data
# set.seed(1)
# res1a3 <- ES(
#   long_data = copy(my_dt),
#   outcomevar = "db_w2_wages",
#   unit_var = "tin",
#   cluster_vars = "tin",
#   cal_time_var = "tax_yr",
#   onset_time_var = "win_yr",
#   omitted_event_time = -2,
#   homogeneous_ATT = FALSE,
#   treated_subset_var = "keep_case",
#   treated_subset_event_time = 0,
#   control_subset_var = "keep_case",
#   control_subset_event_time = 0,
#   ref_discrete_covars = "age2",
#   ref_discrete_covar_event_time = 0)
#
# # should work
# set.seed(1)
# res1a3 <- ES(
#   long_data = copy(my_dt),
#   outcomevar = "db_w2_wages",
#   unit_var = "tin",
#   cluster_vars = "tin",
#   cal_time_var = "tax_yr",
#   onset_time_var = "win_yr",
#   omitted_event_time = -2,
#   homogeneous_ATT = FALSE,
#   treated_subset_var = "keep_case",
#   treated_subset_event_time = 0,
#   control_subset_var = "keep_case",
#   control_subset_event_time = 0,
#   ref_discrete_covars = "age",
#   ref_discrete_covar_event_time = 0)
#
# ###
# # TESTING THAT ref_cont_covars is working correctly
# ###
#
# # should break as 'age2' isn't in the data
# set.seed(1)
# res1a4 <- ES(
#   long_data = copy(my_dt),
#   outcomevar = "db_w2_wages",
#   unit_var = "tin",
#   cluster_vars = "tin",
#   cal_time_var = "tax_yr",
#   onset_time_var = "win_yr",
#   omitted_event_time = -2,
#   homogeneous_ATT = FALSE,
#   treated_subset_var = "keep_case",
#   treated_subset_event_time = 0,
#   control_subset_var = "keep_case",
#   control_subset_event_time = 0,
#   ref_cont_covars = "age2",
#   ref_cont_covar_event_time = 0)
#
# # should work
# set.seed(1)
# res1a4 <- ES(
#   long_data = copy(my_dt),
#   outcomevar = "db_w2_wages",
#   unit_var = "tin",
#   cluster_vars = "tin",
#   cal_time_var = "tax_yr",
#   onset_time_var = "win_yr",
#   omitted_event_time = -2,
#   homogeneous_ATT = FALSE,
#   treated_subset_var = "keep_case",
#   treated_subset_event_time = 0,
#   control_subset_var = "keep_case",
#   control_subset_event_time = 0,
#   ref_cont_covars = "age",
#   ref_cont_covar_event_time = 0)
#
# identical(res1a3, res1a4)

###
# TESTING THAT linearDiD works correctly
###

# linearDiD needs to be logical
set.seed(1)
res1a1 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0,
  linearDiD = "gross_winnings")

# linearDiD needs to be set to TRUE explicitly
set.seed(1)
res1a1 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0,
  linearDiD_treat_var = "gross_winnings")

# linearDiD_treat_var needs to be in the data
set.seed(1)
res1a1 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0,
  linearDiD = TRUE,
  linearDiD_treat_var = "gross_winnings")

# this should work
collapse_table <- data.table(a = c("pre","on_impact", "post_avg"),b = c(list(-7:-1), list(1), list(1:5)))
set.seed(1)
res1a2 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0,
  linearDiD = TRUE,
  linearDiD_treat_var = "gross_winnings_adj_equiv",
  discrete_covars = "age",
  calculate_collapse_estimates = TRUE,
  collapse_inputs = copy(collapse_table))

###
# TESTING THAT THE treated_subset_var / control_subset_var is working correctly
###

set.seed(1)
res1a1 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0)

set.seed(1)
res1a2 <- ES(
  long_data = copy(my_dt),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = -2,
  control_subset_var = "keep_case",
  control_subset_event_time = -2)

# Now let's make it missing in a single event time for a person where ref_event_time matches true event_time
# Done for a non-repeating event-time
my_dt2 <- copy(my_dt)
my_dt2[tin == 11586195 & tax_yr == 2001, age := NA]
my_dt2[, keep_case := between(age, 21, 64, incbounds = TRUE)]
set.seed(1)
res1b1 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0)

set.seed(1)
res1b2 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = -2,
  control_subset_var = "keep_case",
  control_subset_event_time = -2)

# Now let's make it missing in a single event time for a person where ref_event_time matches true event_time
# Done for a repeating event time
my_dt2 <- copy(my_dt)
my_dt2[tin == 11586195 & tax_yr == 1999, age := NA]
my_dt2[, keep_case := between(age, 21, 64, incbounds = TRUE)]
set.seed(1)
res1c1 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0)

set.seed(1)
res1c2 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = -2,
  control_subset_var = "keep_case",
  control_subset_event_time = -2)

# Now let's make it missing in all but a single event time for a person where ref_event_time matches true event_time
# Done for a non-repeating event-time
my_dt2 <- copy(my_dt)
my_dt2[tin == 11586195 & !(tax_yr == 2001), age := NA]
my_dt2[, keep_case := between(age, 21, 64, incbounds = TRUE)]
set.seed(1)
res1d1 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0)

set.seed(1)
res1d2 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = -2,
  control_subset_var = "keep_case",
  control_subset_event_time = -2)

# Now let's make it missing in all but a single event time for a person where ref_event_time matches true event_time
# Done for a repeating event time
my_dt2 <- copy(my_dt)
my_dt2[tin == 11586195 & !(tax_yr == 1999), age := NA]
my_dt2[, keep_case := between(age, 21, 64, incbounds = TRUE)]
set.seed(1)
res1e1 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0)

set.seed(1)
res1e2 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = -2,
  control_subset_var = "keep_case",
  control_subset_event_time = -2)

# Now let's make it all missing for a person where ref_event_time matches true event_time
my_dt2 <- copy(my_dt)
my_dt2[tin == 11586195, age := NA]
my_dt2[, keep_case := between(age, 21, 64, incbounds = TRUE)]
set.seed(1)
res1e1 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = 0,
  control_subset_var = "keep_case",
  control_subset_event_time = 0)

set.seed(1)
res1e2 <- ES(
  long_data = copy(my_dt2),
  outcomevar = "db_w2_wages",
  unit_var = "tin",
  cluster_vars = "tin",
  cal_time_var = "tax_yr",
  onset_time_var = "win_yr",
  omitted_event_time = -2,
  homogeneous_ATT = FALSE,
  treated_subset_var = "keep_case",
  treated_subset_event_time = -2,
  control_subset_var = "keep_case",
  control_subset_event_time = -2)














# # TO DELETE
#
# long_data = copy(my_dt); outcomevar="db_w2_wages"; unit_var="tin"; cal_time_var="tax_yr"; onset_time_var="win_yr"; cluster_vars="tin";
# omitted_event_time= -2; anticipation = 0; min_control_gap=1; max_control_gap=Inf;
# linearize_pretrends=FALSE; fill_zeros=FALSE; residualize_covariates = FALSE;
# control_subset_var=NA; control_subset_event_time=0;
# treated_subset_var=NA; treated_subset_event_time=0;
# control_subset_var2=NA; control_subset_event_time2=0;
# treated_subset_var2=NA; treated_subset_event_time2=0;
# discrete_covars = "age"; cont_covars = NULL; never_treat_action = 'none';
# homogeneous_ATT = FALSE; reg_weights = NULL; add_unit_fes = FALSE;
# bootstrapES = FALSE; bootstrap_iters = 1; bootstrap_num_cores = 1;
# ipw = FALSE; ipw_model = 'linear'; ipw_composition_change = FALSE; ipw_keep_data = FALSE; ipw_ps_lower_bound = 0; ipw_ps_upper_bound = 1;
# event_vs_noevent = FALSE;
# ref_discrete_covars = "age"; ref_discrete_covar_event_time=0; ref_cont_covars = NULL; ref_cont_covar_event_time=0;
# calculate_collapse_estimates = FALSE; collapse_inputs = NULL;
# ref_reg_weights = NULL; ref_reg_weights_event_time=0;
# ntile_var = NULL; ntile_event_time = -2; ntiles = NA; ntile_var_value = NA; ntile_avg = FALSE
#
# # just for ES_clean_data
# endog_var = NULL
#
# # END TO DELETE



#
# ES_clean <- copy(stack_across_cohorts_balanced_treated_control)
# ES_clean_testing <- copy(ES_clean)
# ES_clean_testing[tin == 7955188 | tin == 11586195]
#
# treated_subset_var <- "keep_case"
# control_subset_var <- "keep_case"
#
# # case 1 - no missings
# summary(ES_clean_testing$age)
#
# ES_clean_testing[, keep_case := between(age, 21, 64, incbounds = TRUE)]
# treated_subset_event_time <- 0
# control_subset_event_time <- 0
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
#
# treated_subset_event_time <- -2
# control_subset_event_time <- -2
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
#
# # case 2 - introduce missing in a non-repeated event_time (ref_event_time == 0)
# ES_clean_testing[, age2 := age]
# ES_clean_testing[(tin == 7955188| tin == 11586195) & ref_event_time == 0, age2 := NA]
# ES_clean_testing[, keep_case := between(age2, 21, 64, incbounds = TRUE)]
# treated_subset_event_time <- 0
# control_subset_event_time <- 0
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is missing at the desired event time (but not always missing within c(unit_var, "ref_onset_time"), we get all 0s, as we want)
#
# # now what if given the above missing we were actually interested in -2
# treated_subset_event_time <- -2
# control_subset_event_time <- -2
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # also works correctly
#
# # case 3 - introduce missing in a repeated event time
# ES_clean_testing[, age2 := age]
# ES_clean_testing[(tin == 7955188| tin == 11586195) & ref_event_time == -2, age2 := NA]
# ES_clean_testing[, keep_case := between(age2, 21, 64, incbounds = TRUE)]
# treated_subset_event_time <- 0
# control_subset_event_time <- 0
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is NOT missing at the desired event time, we get the desired output
#
# # now what if given the above missing we were actually interested in -2
# treated_subset_event_time <- -2
# control_subset_event_time <- -2
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is missing at the desired event time (but not always missing within c(unit_var, "ref_onset_time"), we get all 0s, as we want)
#
# # case 4 -- missing everywhere excluding the desired event time, and that is a non-repeated event time
# ES_clean_testing[, age2 := age]
# ES_clean_testing[(tin == 7955188| tin == 11586195) & !(ref_event_time == 0), age2 := NA]
# ES_clean_testing[, keep_case := between(age2, 21, 64, incbounds = TRUE)]
# treated_subset_event_time <- 0
# control_subset_event_time <- 0
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is NOT missing at the desired event time, we get the desired output
#
# # now what if given the above missing we were actually interested in -2
# treated_subset_event_time <- -2
# control_subset_event_time <- -2
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is missing at the desired event time (but not always missing within c(unit_var, "ref_onset_time"), we get all 0s, as we want)
#
# # case 5 -- missing everywhere but the desired event time, and that is a repeated event time
# ES_clean_testing[, age2 := age]
# ES_clean_testing[(tin == 7955188| tin == 11586195) & !(ref_event_time == -2), age2 := NA]
# ES_clean_testing[, keep_case := between(age2, 21, 64, incbounds = TRUE)]
# treated_subset_event_time <- 0
# control_subset_event_time <- 0
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is missing at the desired event time (but not always missing within c(unit_var, "ref_onset_time"), we get all 0s, as we want)
#
# # now what if given the above missing we were actually interested in -2
# treated_subset_event_time <- -2
# control_subset_event_time <- -2
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is NOT missing at the desired event time, we get the desired output
#
# # case 6 -- missing everywhere
# ES_clean_testing[, age2 := age]
# ES_clean_testing[(tin == 7955188| tin == 11586195), age2 := NA]
# ES_clean_testing[, keep_case := between(age2, 21, 64, incbounds = TRUE)]
# treated_subset_event_time <- 0
# control_subset_event_time <- 0
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is NOT missing at the desired event time, we get the desired output
#
# # now what if given the above missing we were actually interested in -2
# treated_subset_event_time <- -2
# control_subset_event_time <- -2
# ES_clean_testing[, c("valid_treated_group","valid_control_group") := NULL]
# ES_clean_testing[, valid_treated_group := as.integer(max2((get(treated_subset_var)*(ref_event_time==treated_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[, valid_control_group := as.integer(max2((get(control_subset_var)*(ref_event_time==control_subset_event_time)), na.rm = TRUE)), by=c(unit_var, "ref_onset_time")]
# ES_clean_testing[tin == 7955188]
# ES_clean_testing[tin == 11586195]
# # ^ as expected, since the variable of interest is missing at the desired event time (but not always missing within c(unit_var, "ref_onset_time"), we get all 0s, as we want)
#
#
