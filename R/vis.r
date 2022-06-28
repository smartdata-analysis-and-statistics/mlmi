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

plot_md_pattern <- function(simpars, censorFUN) {
  
  # Generate a dataset
  simpars$npatients <- 2000
  simpars$ncenters <- 50
  sim_data <- sim_data(simpars)
  
  # Introduce missing values
  sim_data <- censorFUN(data = sim_data, save_prob_y_obs = TRUE)
    
  ggdat <- sim_data %>% group_by(time, trt) %>% 
    summarise(pr_yobs = mean(prob_yobs),
              pr_yobs_iqrl = quantile(prob_yobs, probs = 0.25),
              pr_yobs_iqrm = quantile(prob_yobs, probs = 0.50),
              pr_yobs_iqru = quantile(prob_yobs, probs = 0.75))

  
  # Remove rows where time=0
  ggdat <- subset(ggdat, time > 0)
  
  # Generate factors for relevant variables
  ggdat$Treatment <- factor(ggdat$trt, levels = c(0,1), labels = c("DMT A", "DMT B"))
  
  
  ggplot(ggdat, aes(x = time, y = pr_yobs_iqrm, fill = Treatment, color = Treatment)) + 
    geom_errorbar(aes(group = time, ymax = pr_yobs_iqru, ymin = pr_yobs_iqrl)) +
    geom_point(aes(group = time)) + 
    ylab("Probability of observing the EDSS score") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    xlab("Time (months)") + 
    scale_x_continuous(breaks = c(1, 9, 18, 27, 36, 45, 54, 60)) +
    scale_color_manual(values = c("#2171b5", "#6baed6", "hotpink1", "hotpink4")) +
    theme(legend.position = "bottom") +
    facet_wrap(~Treatment)

}

  