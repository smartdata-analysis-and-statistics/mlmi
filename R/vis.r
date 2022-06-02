library(ggplot2)

plot_example_trajectory <- function(dat,
                                    sel_patid = 1,
                                    outcome_var = "edss",
                                    treat_var = "trt",
                                    time_var = "time") {
  
  dat <- dat %>% filter(patid == sel_patid)
  
  dat$y <- dat[,outcome_var]
  dat$time <- dat[,time_var]
  
  ggplot(dat, aes(x = time, y = y)) +
    geom_line(col = "red") +
    geom_point(aes(x = time, y = y), size = 2) +
    xlab("Time") +
    ylab("EDSS")
}

plot_dens_x <- function(dat,
                        x_var = "age",
                        x_label = "Age",
                        outcome_var = "edss",
                        treat_var = "trt",
                        time_var = "time") {
  
  dat$time <- dat[,time_var]
  dat$trt <- dat[,treat_var]
  dat$x <- dat[,x_var]
  
  dat <- dat %>% group_by(patid,centerid) %>% summarize(x = first(x),
                                                        trt = first(trt))
  
  dat$Treatment <- factor(dat$trt, levels = c(0,1), labels = c("DMT A", "DMT B"))
  
  mean_age <- data.frame(Treatment = c("DMT A", "DMT B"), 
                         wt = c(mean(subset(dat, trt == 0)$x), 
                                mean(subset(dat, trt == 1)$x)))
  
  ggplot(dat, aes(x = x, color = Treatment, fill = Treatment)) + 
    geom_histogram(aes(y=..density..), position="identity", alpha = 0.5, fill = "white") +
    geom_density(alpha = 0.6) +
    facet_grid(Treatment ~ .) +
    geom_vline(aes(xintercept = wt), mean_age, lty = 2) +
    xlab(x_label) + 
    scale_color_brewer(palette = "Paired") + 
    scale_fill_brewer(palette = "Paired") + 
    theme(legend.position = "none")
}

plot_max_fup <- function(dat,
                         x_label = "Age",
                         outcome_var = "edss",
                         treat_var = "trt",
                         time_var = "time") {
  
  dat$time <- dat[,time_var]
  dat$trt <- dat[,treat_var]
  
  dat$Treatment <- factor(dat$trt, levels = c(0,1), labels = c("DMT A", "DMT B"))
  
  plotdat <- dat %>%
    group_by(patid, trt, centerid) %>%
    summarize(max_fup = max(time), Treatment = unique(Treatment))
  
  # Estimate median of max_fup for each center
  mfup <- plotdat %>%
    group_by(centerid) %>%
    summarize(median_mfup = median(max_fup)) %>% 
    arrange(median_mfup)
  
  plotdat$centerid <- factor(plotdat$centerid, labels = mfup$centerid, levels = mfup$centerid)
  
  ggplot(plotdat, aes(x = centerid, y = max_fup, fill = Treatment)) + 
    geom_boxplot() + 
    ylab("Maximum follow-up (months)") +
    xlab("Center") +
    scale_fill_brewer(palette = "Blues")
}


plot_imputed_trajectories <- function(fits, 
                                      sel_patid = 1,
                                      outcome_var = "edss",
                                      treat_var = "trt",
                                      time_var = "time",
                                      col.imputed = "#999999", 
                                      plot.imp = 4) {
  
  dat_pat_imputed <- NULL
  dat_pat_orig <- subset(complete(fits, 0), patid == sel_patid)
  dat_pat_orig <- rename(dat_pat_orig, time = time_var)
  
  for (i in seq(min(plot.imp, fits$m))) {
    dat_pat <- subset(complete(fits, i), patid == sel_patid)
    dat_pat <- rename(dat_pat, y_imputed = outcome_var)
    dat_pat$source <- ifelse(!is.na(dat_pat_orig[,outcome_var]),"Observed", "Imputed")
    dat_pat$imputation <- i
    dat_pat_imputed <- rbind(dat_pat_imputed, dat_pat)
  }
  dat_pat_imputed <- rename(dat_pat_imputed, time = time_var)
  

  g <- ggplot(dat_pat_imputed, aes(x = time, y = y_imputed, group = source)) +
    geom_point(aes(shape = source, color = source), size = 4) +
    scale_color_manual(values = c('#000000', col.imputed )) +
    ylim(0,8) +
    xlab("Days elapsed since treatment start") +
    ylab("EDSS") +
    theme(legend.title = element_blank()) +
    theme(legend.position = "bottom") +
    facet_wrap(~imputation, ncol = 2)
  return(g)
  
  
}
  