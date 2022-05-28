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
  