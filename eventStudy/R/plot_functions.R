#' @export
ES_plot_ATTs <- function(figdata, lower_event = -5, upper_event = 5, ci_factor = 1.96, homogeneous_ATT = FALSE, omitted_event_time = -2, bootstrap_ses = FALSE){

  # As this function changes the underlying results data and figdata tends to be a small file, make a copy now
  figdata_temp <- copy(figdata)

  figdata_temp <- figdata_temp[rn %in% c("att","catt")]
  figdata_temp[, ref_event_time := as.numeric(ref_event_time)]
  figdata_temp <- figdata_temp[ ref_event_time >= lower_event & ref_event_time <= upper_event]
  figdata_temp[, jitter := .GRP, by = ref_onset_time]
  jitter_center <- figdata_temp[,median(unique(jitter))]
  jitter_scale <- figdata_temp[,length(unique(jitter))]
  figdata_temp[, jitter_event_time := ref_event_time + (jitter - jitter_center) / jitter_scale]

  if(bootstrap_ses == TRUE){
    if(!("bootstrap_se" %in% colnames(figdata_temp))){
      stop(print("Variable 'bootstrap_se' is not found in results data you provided. Let me suggest bootstrap_ses = FALSE."))
    } else{
      figdata_temp[, cluster_se := NULL]
      setnames(figdata_temp, "bootstrap_se", "cluster_se")
    }
  }

  if(homogeneous_ATT){
    fig <- ggplot(aes( x = ref_event_time, y = estimate), data = figdata_temp[ref_onset_time=="Pooled"]) +
      geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - ci_factor * cluster_se, ymax = estimate + ci_factor * cluster_se)) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time") +
      annotate(geom = "text", x = omitted_event_time, y = 0, label = "(Omitted)", size = 4)
  } else {
    fig <- ggplot(aes( x = jitter_event_time, y = estimate, colour = factor(ref_onset_time)), data = figdata_temp) +
      geom_point() + theme_bw(base_size = 16) +
      geom_errorbar(aes(ymin = estimate - ci_factor * cluster_se, ymax = estimate + ci_factor * cluster_se)) +
      scale_x_continuous(breaks = pretty_breaks()) +
      labs(x = "Event Time", color = "Cohort") +
      annotate(geom = "text", x = omitted_event_time, y = 0, label = "(Omitted)", size = 4)
  }

  return(fig)

}

#' @export
ES_plot_levels <- function(figdata, cohort_subset = NA, lower_event = -5, upper_event = 5, omitted_event_time = -2){

  # As this function changes the underlying results data and figdata tends to be a small file, make a copy now
  figdata_temp <- copy(figdata)

  figdata_temp <- figdata_temp[rn %in% c("treatment_means")]
  figdata_temp[, estimate := estimate - estimate[ref_event_time == omitted_event_time], list(ref_onset_time,treated)]
  figdata_temp[, ref_onset_time :=  as.integer(as.character(ref_onset_time))]
  figdata_temp[, ref_event_time := as.integer(as.character(ref_event_time))]
  figdata_temp[treated==1, treatment := "Treatment"]
  figdata_temp[treated==0, treatment := "Control"]
  figdata_temp[, treatment := factor(treatment,levels=c("Treatment","Control"))]
  figdata_temp[, year := ref_onset_time + ref_event_time]
  figdata_temp <- figdata_temp[ref_event_time >= lower_event & ref_event_time <= upper_event]

  if(length(cohort_subset) > 0 & !is.na(cohort_subset)){
    figdata_temp <- figdata_temp[ref_onset_time %in% cohort_subset]
  }

  gg <- ggplot(aes(x=year,y=estimate,colour=factor(ref_onset_time),linetype=treatment),data=figdata_temp) +
    geom_line() + geom_point() + geom_vline(aes(xintercept=ref_onset_time,colour=factor(ref_onset_time)),linetype=3) +
    theme_bw(base_size = 16) + labs(colour="Cohort",linetype="",x="Year") +
    scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

  return(gg)

}

#' @export
ES_plot_raw_levels <- function(figdata, cohort_subset = NA, lower_event = -5, upper_event = 5, omitted_event_time = -2){

  # As this function changes the underlying results data and figdata tends to be a small file, make a copy now
  figdata_temp <- copy(figdata)

  figdata_temp <- figdata_temp[rn %in% c("treatment_means")]
  figdata_temp[, ref_onset_time :=  as.integer(as.character(ref_onset_time))]
  figdata_temp[, ref_event_time := as.integer(as.character(ref_event_time))]
  figdata_temp[treated==1, treatment := "Treatment"]
  figdata_temp[treated==0, treatment := "Control"]
  figdata_temp[, treatment := factor(treatment,levels=c("Treatment","Control"))]
  figdata_temp[, year := ref_onset_time + ref_event_time]
  figdata_temp <- figdata_temp[ref_event_time >= lower_event & ref_event_time <= upper_event]

  if(length(cohort_subset) > 0 & !is.na(cohort_subset)){
    figdata_temp <- figdata_temp[ref_onset_time %in% cohort_subset]
  }

  gg <- ggplot(aes(x=year,y=estimate,colour=factor(ref_onset_time),linetype=treatment),data=figdata_temp) +
    geom_line() + geom_point() + geom_vline(aes(xintercept=ref_onset_time,colour=factor(ref_onset_time)),linetype=3) +
    theme_bw(base_size = 16) + labs(colour="Cohort",linetype="",x="Year") +
    scale_x_continuous(breaks = pretty_breaks()) + scale_y_continuous(breaks = pretty_breaks())

  return(gg)

}

#' @export
ES_plot_ipw <- function(figdata,
                        pooled = FALSE,
                        cohort_subset = NA,
                        event_time_subset = NA,
                        type = 'density',
                        binwidth = 0.01,
                        align = "vertical",
                        omitted_event_time = -2){

  figdata_temp <- figdata[!is.na(pr) & ref_event_time != omitted_event_time, list(ref_onset_time, ref_event_time, catt_specific_sample, treated, pr)]
  figdata_temp[, ref_onset_time :=  as.integer(as.character(ref_onset_time))]
  figdata_temp[, ref_event_time := as.integer(as.character(ref_event_time))]
  figdata_temp[treated==1, treatment := "Treatment"]
  figdata_temp[treated==0, treatment := "Control"]
  figdata_temp[, treatment := factor(treatment,levels=c("Control","Treatment"))]

  if(pooled == TRUE){
    if(align == "overlap"){
      if(type == 'density'){
        gg <- ggplot(figdata_temp, aes(x=pr)) +
          geom_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 0], color = 'black', show.legend = FALSE) +
          stat_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 0], color = 'black', geom = "line", position = "identity", size = 0) +
          geom_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 1], color = 'black', show.legend = FALSE) +
          stat_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 1], color = 'black', geom = "line", position = "identity", size = 0) +
          theme_bw(base_size = 12) +
          labs(x = "Estimated Propensity", y = "Density", linetype = "Group") +
          scale_x_continuous(limits = c(0,1), breaks = pretty_breaks()) +
          guides(linetype = guide_legend(override.aes = list(size = 1)))
      } else if(type == "histogram"){
        gg <- ggplot(figdata_temp, aes(x=pr)) +
          scale_fill_grey(start = 0.5, end = 1) +
          geom_histogram(aes(y=..density.., fill=factor(treatment)), binwidth = binwidth, data = figdata_temp[treated == 0], color = NA) +
          geom_histogram(aes(y=..density.., fill=factor(treatment)), binwidth = binwidth, data = figdata_temp[treated == 1], color = 'black', alpha = 0) +
          theme_bw(base_size = 12) +
          labs(x = "Estimated Propensity", y = "Density", fill = "Group") +
          scale_x_continuous(breaks = pretty_breaks()) +
          coord_cartesian(xlim = c(0,1))
      }
    } else if(align == "vertical"){
      if(type == 'density'){
        gg <- ggplot(figdata_temp, aes(x=pr)) +
          geom_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 0], color = 'black', show.legend = FALSE) +
          stat_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 0], color = 'black', geom = "line", position = "identity", size = 0) +
          geom_density(aes(y = -..density.., linetype = factor(treatment)), data = figdata_temp[treated == 1], color = 'black', show.legend = FALSE) +
          stat_density(aes(y = -..density.., linetype = factor(treatment)), data = figdata_temp[treated == 1], color = 'black', geom = "line", position = "identity", size = 0) +
          theme_bw(base_size = 12) +
          labs(x = "Estimated Propensity", y = "Density", linetype = "Group") +
          scale_x_continuous(limits = c(0,1), breaks = pretty_breaks()) +
          guides(linetype = guide_legend(override.aes = list(size = 1)))
      } else if(type == "histogram"){
        gg <- ggplot(figdata_temp, aes(x=pr)) +
          scale_fill_grey(start = 0.5, end = 1) +
          geom_histogram(aes(y=..density.., fill=factor(treatment)), binwidth = binwidth, data = figdata_temp[treated == 0], color = 'white') +
          geom_histogram(aes(y= - ..density.., fill=factor(treatment)), binwidth = binwidth, data = figdata_temp[treated == 1], color = 'black') +
          theme_bw(base_size = 12) +
          labs(x = "Estimated Propensity", y = "Density", fill = "Group") +
          scale_x_continuous(breaks = pretty_breaks()) +
          coord_cartesian(xlim = c(0,1))
      }
    }
  } else{

    if(length(cohort_subset) > 0 & !is.na(cohort_subset)){
      figdata_temp <- figdata_temp[ref_onset_time %in% cohort_subset]
    }

    if(length(event_time_subset) > 0 & !is.na(event_time_subset)){
      figdata_temp <- figdata_temp[ref_event_time %in% event_time_subset]
    }

    if(align == "overlap"){
      if(type == 'density'){
        gg <- ggplot(figdata_temp, aes(x=pr)) +
          geom_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 0], color = 'black', show.legend = FALSE) +
          stat_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 0], color = 'black', geom = "line", position = "identity", size = 0) +
          geom_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 1], color = 'black', show.legend = FALSE) +
          stat_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 1], color = 'black', geom = "line", position = "identity", size = 0) +
          theme_bw(base_size = 12) +
          labs(x = "Estimated Propensity", y = "Density", linetype = "Group") +
          scale_x_continuous(limits = c(0,1), breaks = pretty_breaks()) +
          facet_wrap(~ref_onset_time, scales = c("free")) +
          guides(linetype = guide_legend(override.aes = list(size = 1)))
      } else if(type == "histogram"){
        gg <- ggplot(figdata_temp, aes(x=pr)) +
          scale_fill_grey(start = 0.5, end = 1) +
          geom_histogram(aes(y=..density.., fill=factor(treatment)), binwidth = binwidth, data = figdata_temp[treated == 0], color = NA) +
          geom_histogram(aes(y=..density.., fill=factor(treatment)), binwidth = binwidth, data = figdata_temp[treated == 1], color = 'black', alpha = 0) +
          theme_bw(base_size = 12) +
          labs(x = "Estimated Propensity", y = "Density", fill = "Group") +
          scale_x_continuous(breaks = pretty_breaks()) +
          facet_wrap(~ref_onset_time, scales = c("free")) +
          coord_cartesian(xlim = c(0,1))
      }
    } else if(align == "vertical"){
      if(type == 'density'){
        gg <- ggplot(figdata_temp, aes(x=pr)) +
          geom_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 0], color = 'black', show.legend = FALSE) +
          stat_density(aes(y = ..density.., linetype = factor(treatment)), data = figdata_temp[treated == 0], color = 'black', geom = "line", position = "identity", size = 0) +
          geom_density(aes(y = -..density.., linetype = factor(treatment)), data = figdata_temp[treated == 1], color = 'black', show.legend = FALSE) +
          stat_density(aes(y = -..density.., linetype = factor(treatment)), data = figdata_temp[treated == 1], color = 'black', geom = "line", position = "identity", size = 0) +
          theme_bw(base_size = 12) +
          labs(x = "Estimated Propensity", y = "Density", linetype = "Group") +
          scale_x_continuous(limits = c(0,1), breaks = pretty_breaks()) +
          facet_wrap(~ref_onset_time, scales = c("free")) +
          guides(linetype = guide_legend(override.aes = list(size = 1)))
      } else if(type == "histogram"){
        gg <- ggplot(figdata_temp, aes(x=pr)) +
          scale_fill_grey(start = 0.5, end = 1) +
          geom_histogram(aes(y=..density.., fill=factor(treatment)), binwidth = binwidth, data = figdata_temp[treated == 0], color = 'white') +
          geom_histogram(aes(y= - ..density.., fill=factor(treatment)), binwidth = binwidth, data = figdata_temp[treated == 1], color = 'black') +
          theme_bw(base_size = 12) +
          labs(x = "Estimated Propensity", y = "Density", fill = "Group") +
          scale_x_continuous(breaks = pretty_breaks()) +
          facet_wrap(~ref_onset_time, scales = c("free")) +
          coord_cartesian(xlim = c(0,1))
      }
    }

  }

  return(gg)

}
