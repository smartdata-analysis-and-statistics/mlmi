library(ggplot2)
library(grid)
library(gridExtra)
library(kableExtra)

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

# Visualize data from the case study
plot_baseline_characteristics <- function(data) {
  
  data <- data %>% group_by(mpi) %>% summarize(
    nvisit = max(visit),
    agemspt = first(baseval_age),
    sex = first(baseval_sex),
    educ = first(baseval_educ),
    msdur = first(baseval_msdur),
    mstype = first(baseval_mstype),
    relapses = first(baseval_relapses),
    prior_dmt_effic = first(baseval_prior_dmt_effic),
    cardio = first(baseval_cardio),
    diabetes = first(baseval_diabetes),
    baseval_pdds = first(baseval_pdds),
    dmt_nm = first(dmt_nm)
  )
  
  dat_complete <- subset(data, !is.na(agemspt) & !is.na(sex) & !is.na(educ) & !is.na(msdur) & !is.na(mstype) & !is.na(relapses) & !is.na(prior_dmt_effic) &
                           !is.na(cardio) & !is.na(diabetes) & !is.na(baseval_pdds))
  
  
  d <- data.frame(Parameter = character(),
                  DMF = character(), 
                  FTY = character(), 
                  missing = character())
  
  dhpack <- data.frame("Title" = character(),
                       "index" = numeric())
  
  
  # N
  d <- d %>% add_row(data.frame(Parameter = "Total sample size",
                                DMF = paste(nrow(subset(data, dmt_nm == "Tecfidera")), sep = ""),
                                FTY = paste(nrow(subset(data, dmt_nm == "Gilenya")), sep = ""),
                                missing = ""))
  
  # Complete cases
  d <- d %>% add_row(data.frame(Parameter = "Complete cases", 
                                DMF = paste(nrow(subset(dat_complete, dmt_nm == "Tecfidera")), sep = ""),
                                FTY = paste(nrow(subset(dat_complete, dmt_nm == "Gilenya")), sep = ""),
                                missing = ""))
  dhpack <- dhpack %>% add_row(data.frame("Title" = "Sample size, n", "index" = 2))
  
  
  # Age
  age_t0 <- subset(data, dmt_nm == "Gilenya")$agemspt
  age_t1 <- subset(data, dmt_nm == "Tecfidera")$agemspt
  d <- d %>% add_row(data.frame(Parameter = "Median, years (IQR)",
                                DMF = paste(round(median(age_t1, na.rm = T),0), " [", round(quantile(age_t1, na.rm = T)["25%"],0), " - ", round(quantile(age_t1, na.rm = T)["75%"],0), "]", sep = ""),
                                FTY = paste(round(median(age_t0, na.rm = T),0), " [", round(quantile(age_t0, na.rm = T)["25%"],0), " - ", round(quantile(age_t0, na.rm = T)["75%"],0), "]", sep = ""),
                                missing = paste(sum(is.na(data$agemspt)), " [", round(sum(is.na(data$agemspt))*100/nrow(data),0), "%]")))
  dhpack <- dhpack %>% add_row(data.frame("Title" = "Age", "index" = 1))
  
  # Gender (0 = female, 1 = male)
  sex_t0 <- subset(data, dmt_nm == "Gilenya")$sex
  sex_t1 <- subset(data, dmt_nm == "Tecfidera")$sex
  d <- d %>% add_row(data.frame(Parameter = "Male, n (%)",
                                DMF = paste(sum(sex_t1 == 1), " [", round(sum(sex_t1 == 1)*100/sum(!is.na(sex_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(sex_t0 == 1), " [", round(sum(sex_t0 == 1)*100/sum(!is.na(sex_t0)), 0), "%]", sep = ""),
                           missing = paste(sum(is.na(data$sex)), " [", round(sum(is.na(data$sex))*100/nrow(data),0), "%]")))
  d <- d %>% add_row(data.frame(Parameter = "Female, n (%)",
                                DMF = paste(sum(sex_t1 == 0), " [", round(sum(sex_t1 == 0)*100/sum(!is.na(sex_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(sex_t0 == 0), " [", round(sum(sex_t0 == 0)*100/sum(!is.na(sex_t0)), 0), "%]", sep = ""),
                           missing = paste(sum(is.na(data$sex)), " [", round(sum(is.na(data$sex))*100/nrow(data),0), "%]")))
  dhpack <- dhpack %>% add_row(data.frame("Title" = "Gender", "index" = 2))
  
  # Years of education
  educ_t0 <- subset(data, dmt_nm == "Gilenya")$educ
  educ_t1 <- subset(data, dmt_nm == "Tecfidera")$educ
  d <- d %>% add_row(data.frame(Parameter = "Median, years (IQR)",
                                DMF = paste(round(median(educ_t1, na.rm = T),0), " [", round(quantile(educ_t1, na.rm = T)["25%"],0), " - ", round(quantile(educ_t1, na.rm = T)["75%"],0), "]", sep = ""),
                           FTY = paste(round(median(educ_t0, na.rm = T),0), " [", round(quantile(educ_t0, na.rm = T)["25%"],0), " - ", round(quantile(educ_t0, na.rm = T)["75%"],0), "]", sep = ""),
                           missing = paste(sum(is.na(data$educ)), " [", round(sum(is.na(data$educ))*100/nrow(data),0), "%]")))
  dhpack <- dhpack %>% add_row(data.frame("Title" = "Years of education", "index" = 1))
  
  
  # Disease duration
  msdur_t0 <- subset(data, dmt_nm == "Gilenya")$msdur
  msdur_t1 <- subset(data, dmt_nm == "Tecfidera")$msdur
  d <- d %>% add_row(data.frame(Parameter = "Median, years (IQR)",
                                DMF = paste(round(median(msdur_t1, na.rm = T),0), " [", round(quantile(msdur_t1, na.rm = T)["25%"],0), " - ", round(quantile(msdur_t1, na.rm = T)["75%"],0), "]", sep = ""),
                           FTY = paste(round(median(msdur_t0, na.rm = T),0), " [", round(quantile(msdur_t0, na.rm = T)["25%"],0), " - ", round(quantile(msdur_t0, na.rm = T)["75%"],0), "]", sep = ""),
                           missing = paste(sum(is.na(data$msdur)), " [", round(sum(is.na(data$msdur))*100/nrow(data),0), "%]")))
  dhpack <- dhpack %>% add_row(data.frame("Title" = "Disease duration", "index" = 1))
  
  # MS type
  mstype_t0 <- subset(data, dmt_nm == "Gilenya")$mstype
  mstype_t1 <- subset(data, dmt_nm == "Tecfidera")$mstype
  d <- d %>% add_row(data.frame(Parameter = "Relapsing MS (remitting / progressive)",
                                DMF = paste(sum(mstype_t1 %in% c("Relapsing Remitting MS", "Progressive Relapsing MS"), na.rm = T), " [", round(sum(mstype_t1 %in% c("Relapsing Remitting MS", "Progressive Relapsing MS"), na.rm = T)*100/sum(!is.na(mstype_t1)), 0), "%]", sep = ""),
                                FTY = paste(sum(mstype_t0 %in% c("Relapsing Remitting MS", "Progressive Relapsing MS"), na.rm = T), " [", round(sum(mstype_t0 %in% c("Relapsing Remitting MS", "Progressive Relapsing MS"), na.rm = T)*100/sum(!is.na(mstype_t0)), 0), "%]", sep = ""),
                                missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "Primary Progressive MS",
                                DMF = paste(sum(mstype_t1 == "Primary Progressive MS", na.rm = T), " [", round(sum(mstype_t1 == "Primary Progressive MS", na.rm = T)*100/sum(!is.na(mstype_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(mstype_t0 == "Primary Progressive MS", na.rm = T), " [", round(sum(mstype_t0 == "Primary Progressive MS", na.rm = T)*100/sum(!is.na(mstype_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "Secondary Progressive MS",
                                DMF = paste(sum(mstype_t1 == "Secondary Progressive MS", na.rm = T), " [", round(sum(mstype_t1 == "Secondary Progressive MS", na.rm = T)*100/sum(!is.na(mstype_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(mstype_t0 == "Secondary Progressive MS", na.rm = T), " [", round(sum(mstype_t0 == "Secondary Progressive MS", na.rm = T)*100/sum(!is.na(mstype_t0)), 0), "%]", sep = ""),
                           missing = ""))
  dhpack <- dhpack %>% add_row(data.frame("Title" = "MS type, n (%)", "index" = 3))
  
  # Number of relapses in previous year
  relapses_t0 <- subset(data, dmt_nm == "Gilenya")$relapses
  relapses_t1 <- subset(data, dmt_nm == "Tecfidera")$relapses

  d <- d %>% add_row(data.frame(Parameter = "0",
                                DMF = paste(sum(relapses_t1 == "0", na.rm = T), " [", round(sum(relapses_t1 == "0", na.rm = T)*100/sum(!is.na(relapses_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(relapses_t0 == "0", na.rm = T), " [", round(sum(relapses_t0 == "0", na.rm = T)*100/sum(!is.na(relapses_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "1",
                                DMF = paste(sum(relapses_t1 == "1", na.rm = T), " [", round(sum(relapses_t1 == "1", na.rm = T)*100/sum(!is.na(relapses_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(relapses_t0 == "1", na.rm = T), " [", round(sum(relapses_t0 == "1", na.rm = T)*100/sum(!is.na(relapses_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "2",
                                DMF = paste(sum(relapses_t1 == "2", na.rm = T), " [", round(sum(relapses_t1 == "2", na.rm = T)*100/sum(!is.na(relapses_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(relapses_t0 == "2", na.rm = T), " [", round(sum(relapses_t0 == "2", na.rm = T)*100/sum(!is.na(relapses_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = ">=3",
                                DMF = paste(sum(relapses_t1 == "3_plus", na.rm = T), " [", round(sum(relapses_t1 == "3_plus", na.rm = T)*100/sum(!is.na(relapses_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(relapses_t0 == "3_plus", na.rm = T), " [", round(sum(relapses_t0 == "3_plus", na.rm = T)*100/sum(!is.na(relapses_t0)), 0), "%]", sep = ""),
                           missing = ""))
  
  dhpack <- dhpack %>% add_row(data.frame("Title" = "Number of relapses in previous year, n (%)", "index" = 4))
  
  # Primary DMT efficacy in previous year, n (%)
  prior_dmt_effic_t0 <- subset(data, dmt_nm == "Gilenya")$prior_dmt_effic
  prior_dmt_effic_t1 <- subset(data, dmt_nm == "Tecfidera")$prior_dmt_effic
  d <- d %>% add_row(data.frame(Parameter = "High",
                                DMF = paste(sum(prior_dmt_effic_t1 == "High efficacy", na.rm = T), " [", round(sum(prior_dmt_effic_t1 == "High efficacy", na.rm = T)*100/sum(!is.na(prior_dmt_effic_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(prior_dmt_effic_t0 == "High efficacy", na.rm = T), " [", round(sum(prior_dmt_effic_t0 == "High efficacy", na.rm = T)*100/sum(!is.na(prior_dmt_effic_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "Medium",
                                DMF = paste(sum(prior_dmt_effic_t1 == "Medium efficacy", na.rm = T), " [", round(sum(prior_dmt_effic_t1 == "Medium efficacy", na.rm = T)*100/sum(!is.na(prior_dmt_effic_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(prior_dmt_effic_t0 == "Medium efficacy", na.rm = T), " [", round(sum(prior_dmt_effic_t0 == "Medium efficacy", na.rm = T)*100/sum(!is.na(prior_dmt_effic_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "Low",
                                DMF = paste(sum(prior_dmt_effic_t1 == "Low efficacy", na.rm = T), " [", round(sum(prior_dmt_effic_t1 == "Low efficacy", na.rm = T)*100/sum(!is.na(prior_dmt_effic_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(prior_dmt_effic_t0 == "Low efficacy", na.rm = T), " [", round(sum(prior_dmt_effic_t0 == "Low efficacy", na.rm = T)*100/sum(!is.na(prior_dmt_effic_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "None",
                                DMF = paste(sum(prior_dmt_effic_t1 == "No prior DMT", na.rm = T), " [", round(sum(prior_dmt_effic_t1 == "No prior DMT", na.rm = T)*100/sum(!is.na(prior_dmt_effic_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(prior_dmt_effic_t0 == "No prior DMT", na.rm = T), " [", round(sum(prior_dmt_effic_t0 == "No prior DMT", na.rm = T)*100/sum(!is.na(prior_dmt_effic_t0)), 0), "%]", sep = ""),
                           missing = ""))
  
  dhpack <- dhpack %>% add_row(data.frame("Title" = "Primary DMT efficacy in previous year, n (%)", "index" = 4))
  
  
  # History of cardiovascular disease, n (%)
  cvd_t0 <- subset(data, dmt_nm == "Gilenya")$cardio
  cvd_t1 <- subset(data, dmt_nm == "Tecfidera")$cardio
  d <- d %>% add_row(data.frame(Parameter = "Cardiovascular disease",
                                DMF = paste(sum(cvd_t1 == 1, na.rm = T), " [", round(sum(cvd_t1 == 1, na.rm = T)*100/sum(!is.na(cvd_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(cvd_t0 == 1, na.rm = T), " [", round(sum(cvd_t0 == 1, na.rm = T)*100/sum(!is.na(cvd_t0)), 0), "%]", sep = ""),
                           missing = paste(sum(is.na(data$cardio)), " [", round(sum(is.na(data$cardio))*100/nrow(data),0), "%]")))
  
  # History of diabetes, n (%)
  dia_t0 <- subset(data, dmt_nm == "Gilenya")$diabetes
  dia_t1 <- subset(data, dmt_nm == "Tecfidera")$diabetes
  d <- d %>% add_row(data.frame(Parameter = "Diabetes",
                                DMF = paste(sum(dia_t1 == 1, na.rm = T), " [", round(sum(dia_t1 == 1, na.rm = T)*100/sum(!is.na(dia_t1)), 0), "%]", sep = ""),
                           FTY = paste(sum(dia_t0 == 1, na.rm = T), " [", round(sum(dia_t0 == 1, na.rm = T)*100/sum(!is.na(dia_t0)), 0), "%]", sep = ""),
                           missing = paste(sum(is.na(data$diabetes)), " [", round(sum(is.na(data$diabetes))*100/nrow(data),0), "%]")))
  
  dhpack <- dhpack %>% add_row(data.frame("Title" = "Medical history, n (%)", "index" = 2))
  
  # PDDS score, n (%)
  pdds_t0 <- subset(data, dmt_nm == "Gilenya")$baseval_pdds
  pdds_t1 <- subset(data, dmt_nm == "Tecfidera")$baseval_pdds
  d <- d %>% add_row(data.frame(Parameter = "0 - 1",
                                DMF = paste(length(which(pdds_t1 <= 1)), " [", round(length(which(pdds_t1 <= 1))*100/sum(!is.na(pdds_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(pdds_t0 <= 1)), " [", round(length(which(pdds_t0 <= 1))*100/sum(!is.na(pdds_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "2 - 3",
                                DMF = paste(length(which(pdds_t1 >= 2 & pdds_t1 <= 3)), " [", round(length(which(pdds_t1 >= 2 & pdds_t1 <= 3))*100/sum(!is.na(pdds_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(pdds_t0 >= 2 & pdds_t0 <= 3)), " [", round(length(which(pdds_t0 >= 2 & pdds_t0 <= 3))*100/sum(!is.na(pdds_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "4 - 5",
                                DMF = paste(length(which(pdds_t1 >= 4 & pdds_t1 <= 5)), " [", round(length(which(pdds_t1 >= 4 & pdds_t1 <= 5))*100/sum(!is.na(pdds_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(pdds_t0 >= 4 & pdds_t0 <= 5)), " [", round(length(which(pdds_t0 >=4 & pdds_t0 <= 5))*100/sum(!is.na(pdds_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = ">= 6",
                                DMF = paste(length(which(pdds_t1 >= 6 )), " [", round(length(which(pdds_t1 >= 6))*100/sum(!is.na(pdds_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(pdds_t0 >= 6)), " [", round(length(which(pdds_t0 >= 6 ))*100/sum(!is.na(pdds_t0)), 0), "%]", sep = ""),
                           missing = ""))
  
  dhpack <- dhpack %>% add_row(data.frame("Title" = "PDDS score, n (%)", "index" = 4))

  # Number of visits
  nvisit_t0 <- subset(data, dmt_nm == "Gilenya")$nvisit
  nvisit_t1 <- subset(data, dmt_nm == "Tecfidera")$nvisit
  
  d <- d %>% add_row(data.frame(Parameter = "1 visit",
                                DMF = paste(length(which(nvisit_t1 == 1)), " [", round(length(which(nvisit_t1 == 1))*100/sum(!is.na(nvisit_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(nvisit_t0 == 1)), " [", round(length(which(nvisit_t0 == 1))*100/sum(!is.na(nvisit_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "2 visits",
                                DMF = paste(length(which(nvisit_t1 == 2)), " [", round(length(which(nvisit_t1 == 2))*100/sum(!is.na(nvisit_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(nvisit_t0 == 2)), " [", round(length(which(nvisit_t0 == 2))*100/sum(!is.na(nvisit_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "3 visits",
                                DMF = paste(length(which(nvisit_t1 == 3)), " [", round(length(which(nvisit_t1 == 3))*100/sum(!is.na(nvisit_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(nvisit_t0 == 3)), " [", round(length(which(nvisit_t0 == 3))*100/sum(!is.na(nvisit_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = "4 visits",
                                DMF = paste(length(which(nvisit_t1 == 4)), " [", round(length(which(nvisit_t1 == 4))*100/sum(!is.na(nvisit_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(nvisit_t0 == 4)), " [", round(length(which(nvisit_t0 == 4))*100/sum(!is.na(nvisit_t0)), 0), "%]", sep = ""),
                           missing = ""))
  d <- d %>% add_row(data.frame(Parameter = ">= 5 visits",
                                DMF = paste(length(which(nvisit_t1 >= 5)), " [", round(length(which(nvisit_t1 >= 5))*100/sum(!is.na(nvisit_t1)), 0), "%]", sep = ""),
                           FTY = paste(length(which(nvisit_t0 >= 5)), " [", round(length(which(nvisit_t0 >= 5))*100/sum(!is.na(nvisit_t0)), 0), "%]", sep = ""),
                           missing = ""))
  

  dhpack <- dhpack %>% add_row(data.frame("Title" = "Visit count, n (%)", "index" = 5))
  
  dat_groups <- dhpack %>% pull("index")
  names(dat_groups) <- dhpack %>% pull("Title")
  
  d %>% dplyr::select(Parameter, DMF, FTY) %>% 
    knitr::kable(
      format = "latex",
      align = "lcc",
      booktabs = TRUE,
      longtable = TRUE,
      linesep = ""
    ) %>%
    kableExtra::kable_styling(
      latex_options = c("repeat_header"),
    ) %>%
    pack_rows(index = dat_groups)
}

plot_interval_visits <- function(data) {
  
  # Remove patients with only one visit
  pdat <- data %>%
    group_by(mpi, dmt_nm) %>%
    summarize(numvisits = max(visit))
  
  plotdat <- subset(data, (mpi %in% subset(pdat, numvisits > 1)$mpi))
  plotdat$difftime[2:nrow(plotdat)] <- plotdat$time[-1] - plotdat$time[-nrow(plotdat)]
  plotdat$difftime[plotdat$difftime < 0] <- NA
  plotdat <- subset(plotdat, !is.na(difftime))
  
  ggplot(plotdat, aes(x = difftime, color = treatment)) +
    geom_histogram(position = "dodge", aes(fill = treatment)) +
    xlab("Intervisit time (#days)") +
    scale_color_brewer(palette = "Paired") +
    scale_fill_brewer(palette = "Paired") +
    ggtitle(paste("Interval between consecutive visits; N =", length(unique(plotdat$mpi)), "patients"))  + 
    theme(legend.position = "bottom")
 
  
}
  