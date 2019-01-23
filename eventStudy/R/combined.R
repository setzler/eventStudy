

#' @export
ES <- function(long_data, outcomevar, unit_var, cal_time_var, onset_time_var, cluster_vars,
               omitted_event_time= -2, min_control_gap=1, max_control_gap=Inf,
               control_subset_var=NA, control_subset_event_time=0, fill_zeros=F, vars_to_fill=NULL){

  # fill with zeros
  if(fill_zeros){
    long_data <- ES_expand_to_balance(long_data, vars_to_fill = vars_to_fill, unit_var, cal_time_var)
  }

  # process data
  ES_data <- ES_clean_data(long_data = long_data, outcomevar = outcomevar,
                           unit_var = unit_var, cal_time_var = cal_time_var, onset_time_var = onset_time_var,
                           min_control_gap = min_control_gap, max_control_gap = max_control_gap, omitted_event_time = omitted_event_time,
                           control_subset_var = control_subset_var, control_subset_event_time = control_subset_event_time)

  # collect ATT estimates
  ES_results_hetero <- ES_estimate_ATT(ES_data = ES_data, outcomevar=outcomevar, onset_time_var = onset_time_var, cluster_vars = cluster_vars, homogeneous_ATT = FALSE)
  ES_results_homo <- ES_estimate_ATT(ES_data = ES_data, outcomevar=outcomevar, onset_time_var = onset_time_var, cluster_vars = cluster_vars, homogeneous_ATT = TRUE)
  setnames(ES_results_hetero,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))
  setnames(ES_results_homo,c(onset_time_var,"event_time"),c("ref_onset_time","ref_event_time"))

  # collect levels by treatment/control
  ES_treatcontrol_means <- ES_data[,list(rn="treatment_means", estimate = mean(get(outcomevar)), cluster_se = sd(get(outcomevar))/sqrt(.N)),list(ref_onset_time,ref_event_time,treated)][order(ref_onset_time,ref_event_time,treated)]

  # export
  figdata <- rbindlist(list(ES_results_hetero, ES_results_homo, ES_treatcontrol_means), use.names = TRUE, fill=TRUE)
  figdata[is.na(cluster_se), cluster_se := 0]
  return(figdata)
}

#' @export
ES_plot_ATTs <- function(figdata, lower_event = -3, upper_event = 5, ci_factor = 1.96, homogeneous_only = FALSE){

  figdata <- figdata[rn %in% c("att","catt")]
  figdata[, ref_event_time := as.numeric(ref_event_time)]
  figdata <- figdata[ ref_event_time >= lower_event & ref_event_time <= upper_event]
  figdata[, jitter := .GRP, by = ref_onset_time]
  jitter_center <- figdata[,median(unique(jitter))]
  jitter_scale <- figdata[,length(unique(jitter))]
  figdata[, jitter_event_time := ref_event_time + (jitter - jitter_center) / jitter_scale]

  if(homogeneous_only){
    fig <- ggplot(aes( x = ref_event_time, y = estimate), data = figdata[ref_onset_time=="Pooled"]) +
      geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - ci_factor * cluster_se, ymax = estimate + ci_factor * cluster_se)) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time") +
      annotate(geom = "text", x = -2, y = 0, label = "(Omitted)", size = 4)
  } else {
    fig <- ggplot(aes( x = jitter_event_time, y = estimate, colour = factor(ref_onset_time)), data = figdata) +
      geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - ci_factor * cluster_se, ymax = estimate + ci_factor * cluster_se)) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time", color = "Cohort") +
      annotate(geom = "text", x = -2, y = 0, label = "(Omitted)", size = 4)
  }

  return(fig)

}

#' @export
ES_plot_levels <- function(figdata, cohort_subset = NA, lower_event = -3, upper_event = 5, omitted_event_time = -2){

  figdata <- figdata[rn %in% c("treatment_means")]
  figdata[, estimate := estimate - estimate[ref_event_time == omitted_event_time], list(ref_onset_time,treated)]
  figdata[, ref_onset_time :=  as.integer(as.character(ref_onset_time))]
  figdata[, ref_event_time := as.integer(as.character(ref_event_time))]
  figdata[treated==1, treatment := "Treatment"]
  figdata[treated==0, treatment := "Control"]
  figdata[, treatment := factor(treatment,levels=c("Treatment","Control"))]
  figdata[, year := ref_onset_time + ref_event_time]
  figdata <- figdata[ref_event_time >= lower_event & ref_event_time <= upper_event]

  if(!is.na(cohort_subset)){
    figdata <- figdata[ref_onset_time %in% cohort_subset]
  }

  gg <- ggplot(aes(x=year,y=estimate,colour=factor(ref_onset_time),linetype=treatment),data=figdata) +
    geom_line() + geom_point() + geom_vline(aes(xintercept=ref_onset_time,colour=factor(ref_onset_time)),linetype=3) +
    theme_bw(base_size = 16) + labs(colour="Cohort",linetype="",x="Year") +
    scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

  return(gg)

}




