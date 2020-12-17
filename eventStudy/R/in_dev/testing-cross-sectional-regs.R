library(eventStudy)

options(scipen = 999)
options(show.error.locations=TRUE)
options(warn=1)

##########################################################
# Testing; MAKE SURE TO CLEAN AND REBUILD BEFORE RUNNING!
##########################################################

RNGkind("L'Ecuyer-CMRG")
seed = 52649583L
set.seed(seed)

# simulate data
# mydt <- ES_simulate_data(units = 5000,
#                          time_vary_confounds_high_dim = TRUE)[[1]]
# mydt <- mydt[, list(tin,tax_yr,win_yr,outcome,time_vary_var_high_dim)]
# setnames(mydt,c("time_vary_var_high_dim"), c("age"))
mydt <- ES_simulate_data(units = 5000)[[1]]
mydt <- mydt[, list(tin,tax_yr,win_yr,outcome)]

setorderv(mydt, c("tin", "tax_yr"))
summary(mydt)

#########################
# DiD with Later Winners
#########################

# mydt[, age_case := as.logical(between(age, 21, 64, incbounds = TRUE))]
omitted_event_time  <- -2
anticipation <- 0
collapse_table <- data.table(a = c("pre_3",
                                   "one_to_two",
                                   "one_to_three",
                                   "one_to_four",
                                   "two_to_three",
                                   "two_to_four",
                                   "three_to_four"),
                             b = c(list(-3:(-1 - anticipation)),
                                   list(1:2),
                                   list(1:3),
                                   list(1:4),
                                   list(2:3),
                                   list(2:4),
                                   list(3:4)))

res1_new <- ES(long_data = copy(mydt),
               outcomevar = "outcome",
               unit_var = "tin",
               cal_time_var = "tax_yr",
               onset_time_var = "win_yr",
               cluster_vars = "tin",
               omitted_event_time = omitted_event_time,
               calculate_collapse_estimates = TRUE,
               collapse_inputs = collapse_table,
               cross_sectional_spec = TRUE)

omitted_diffs_t <- res1_new[rn == "treatment_means" & ref_event_time == -2 & treated == 1,
                            list(ref_onset_time,
                                 ref_event_time,
                                 estimate)]
setnames(omitted_diffs_t, "estimate", "t_estimate")
omitted_diffs_c <- res1_new[rn == "treatment_means" & ref_event_time == -2 & treated == 0,
                            list(ref_onset_time,
                                 ref_event_time,
                                 estimate)]
setnames(omitted_diffs_c, "estimate", "c_estimate")
omitted_diffs <- merge(omitted_diffs_t, omitted_diffs_c,
                       by = c("ref_onset_time","ref_event_time"),
                       sort = FALSE)
omitted_diffs[, estimate := (t_estimate - c_estimate)]

res1_new[rn == "catt" & ref_event_time == -2, list(ref_onset_time, ref_event_time, estimate)]

omitted_diffs_t1 <- res1_new[rn == "treatment_means" & ref_event_time == 1 & treated == 1,
                            list(ref_onset_time,
                                 ref_event_time,
                                 estimate)]
setnames(omitted_diffs_t1, "estimate", "t_estimate")
omitted_diffs_c1 <- res1_new[rn == "treatment_means" & ref_event_time == 1 & treated == 0,
                            list(ref_onset_time,
                                 ref_event_time,
                                 estimate)]
setnames(omitted_diffs_c1, "estimate", "c_estimate")
omitted_diffs1 <- merge(omitted_diffs_t1, omitted_diffs_c1,
                       by = c("ref_onset_time","ref_event_time"),
                       sort = FALSE)
omitted_diffs1[, estimate := (t_estimate - c_estimate)]

res1_new[rn == "catt" & ref_event_time == 1, list(ref_onset_time, ref_event_time, estimate)]

# res1_old <- ES(long_data = copy(mydt),
#                outcomevar = "outcome",
#                unit_var = "tin",
#                cal_time_var = "tax_yr",
#                onset_time_var = "win_yr",
#                cluster_vars = "tin",
#                discrete_covars = "age",
#                omitted_event_time = omitted_event_time,
#                control_subset_var = "age_case",
#                control_subset_event_time = 0,
#                treated_subset_var = "age_case",
#                treated_subset_event_time = 0,
#                calculate_collapse_estimates = TRUE,
#                collapse_inputs = collapse_table)
#
# res1_new <- ES(long_data = copy(mydt),
#                outcomevar = "outcome",
#                unit_var = "tin",
#                cal_time_var = "tax_yr",
#                onset_time_var = "win_yr",
#                cluster_vars = "tin",
#                discrete_covars = "age",
#                omitted_event_time = omitted_event_time,
#                control_subset_var = "age_case",
#                control_subset_event_time = 0,
#                treated_subset_var = "age_case",
#                treated_subset_event_time = 0,
#                calculate_collapse_estimates = TRUE,
#                collapse_inputs = collapse_table,
#                cross_sectional_spec = TRUE)
#
# # Do they agree on all of the treatment_control_means (excluding ref_event_time == omitted_event_time)
#
# res1_old_tc_means <- res1_old[rn == "treatment_means" & ref_event_time != omitted_event_time]
# res1_new_tc_means <- res1_new[rn == "treatment_means" & ref_event_time != omitted_event_time]
#
# identical(setorderv(res1_old_tc_means, c("ref_onset_time", "ref_event_time")),
#           setorderv(res1_new_tc_means, c("ref_onset_time", "ref_event_time"))) # TRUE
