library(ipwCoxCSV)
library(stringr)

derive_max_fup <- function(data,
                           outcome_var = "edss") {
  
  data$y <- data[,outcome_var]
  
  data <- data %>%
    group_by(patid) %>%
    mutate(max_fup = max(time[!is.na(y)]))
  
  data <- dplyr::select(data, -y)
  
  data.frame(data)
}




evaluate_rmse <- function(data,
                          outcome_var = "edss",  # Variable that contains the incomplete data
                          orig_outcome_var = "edss_orig", # Variable that contains the original data
                          y_label, # Variable that contains the imputed data
                          window_size = 1) {
  
  data <- rename(data, y = outcome_var)
  
  # Restrict the dataset to the follow-up and grid
  impdata <- subset(data, time <= max_fup & time %% window_size == 0 & is.na(y))
  
  rmse <- sqrt(mean((impdata[,y_label] - impdata[,orig_outcome_var])**2))
  data.frame(nobs_imputed = nrow(impdata), rmse = rmse)
}

evaluate_reference <- function(data,
                               window_size = 3,
                               confirmation_window = 1) {
  
  # Derive TTE set for the original data
  prog_data <- derive_cdp(data = data, 
                          outcome_var = "edss_orig",
                          window_size = window_size,
                          confirmation_window = confirmation_window,
                          extrapolate = "all") 
  
  # Derive treatment effect for the original data
  fit <- ipwCoxCluster(data.frame(prog_data), 
                       indID = "centerid", 
                       indA = "x", 
                       indX = c("age", "edss_baseval"), 
                       indStatus = "event", 
                       indTime = "tte", 
                       ties = "breslow",
                       confidence = 0.95)
  
  
  # Return the effect size and standard error
  data.frame("method" = "Reference",
             "est_logHR" = fit["conventional weights","log HR Estimate"], # Estimated logHR
             "est_HR" = fit["conventional weights","HR Estimate"], # Estimated HR
             "est_HR_CIl" = fit["conventional weights","HR 95% CI-low"], 
             "est_HR_CIu" = fit["conventional weights","HR 95% CI-up"],  
             "window" = window_size,
             "confirmation" = confirmation_window,
             "nobs_imputed" = 0,
             "rmse" = 0)
}

evaluate_locf <- function(data, 
                          window_size, 
                          confirmation_window,
                          extrapolate_cdp,
                          method_lbl = "LOCF") 
  {
  
  # Generate imputed dataset according to LOCF
  impdata <- impute_y_locf(data) 
  
  # Identify number of imputed values and their RMSE
  q_imp <- evaluate_rmse(impdata, y_label =  "y_imp_locf", window_size = window_size)
  
  # Convert longitudinal data to TTE dataset
  prog_data <- derive_cdp(data = impdata,
                          outcome_var = "y_imp_locf",
                          window_size = window_size, 
                          confirmation_window = confirmation_window,
                          extrapolate = extrapolate_cdp) 
  
  # Estimate the marginal treatment effect
  fit <- ipwCoxCluster(data.frame(prog_data), 
                       indID = "centerid", 
                       indA = "x", 
                       indX = c("age", "edss_baseval"), 
                       indStatus = "event", 
                       indTime = "tte", 
                       ties = "breslow",
                       confidence = 0.95)

  
  # Return the effect size and standard error
  data.frame("method" = method_lbl,
             "est_logHR" = fit["conventional weights","log HR Estimate"], # Estimated logHR
             "est_HR" = fit["conventional weights","HR Estimate"], # Estimated logHR
             "est_HR_CIl" = fit["conventional weights","HR 95% CI-low"], 
             "est_HR_CIu" = fit["conventional weights","HR 95% CI-up"],  
             "window" = window_size, 
             "confirmation" = confirmation_window,
             "nobs_imputed" = q_imp$nobs_imputed,
             "rmse" = q_imp$rmse)
}

evaluate_rounding <- function(data, 
                              window_size, 
                              confirmation_window,
                              extrapolate_cdp,
                              method_lbl = "Rounding"
) {
  
  # Generate imputed dataset according to Rounding
  impdata <- impute_y_rounding(data)
  
  # Identify number of imputed values and their RMSE
  q_imp <- evaluate_rmse(impdata, y_label =  "y_imp_rounding", window_size = window_size)
  
  # Convert longitudinal data to TTE dataset
  prog_data <- derive_cdp(data = impdata, 
                          outcome_var = "y_imp_rounding",
                          window_size = window_size,
                          confirmation_window = confirmation_window,
                          extrapolate = extrapolate_cdp) 
  
  # Estimate the marginal treatment effect
  fit <- ipwCoxCluster(data.frame(prog_data), 
                       indID = "centerid", 
                       indA = "x", 
                       indX = c("age", "edss_baseval"), 
                       indStatus = "event", 
                       indTime = "tte", 
                       ties = "breslow",
                       confidence = 0.95)
  

  # Return the effect size and standard error
  data.frame("method" = method_lbl,
             "est_logHR" = fit["conventional weights","log HR Estimate"], # Estimated logHR
             "est_HR" = fit["conventional weights","HR Estimate"], # Estimated logHR
             "est_HR_CIl" = fit["conventional weights","HR 95% CI-low"], 
             "est_HR_CIu" = fit["conventional weights","HR 95% CI-up"],  
             "window" = window_size, 
             "confirmation" = confirmation_window,
             "nobs_imputed" = q_imp$nobs_imputed,
             "rmse" = q_imp$rmse)
}

evaluate_lme_3l <- function(data,
                            n.imp = 10,
                            window_size, 
                            confirmation_window,
                            extrapolate_cdp,
                            sampling_resid = "condMVN",
                            method_lbl = "LME-CE (PMM)"
) {
  
  if (missing(n.imp) | n.imp < 1) {
    stop("Invalid number of imputations")
  }
  
  impdata <- impute_y_mice_3l(data = data, n.imp = n.imp, sampling_resid = sampling_resid)
  
  if (is.null(impdata$fit)) {
    out <- data.frame("method" = c(paste(method_lbl, "(PMM)", paste(method_lbl, "(EDSS conversion)"))),
                      "est_logHR" = c(NA, NA), # Estimated logHR
                      "est_HR" = c(NA, NA), # Estimated logHR
                      "est_HR_CIl" = c(NA, NA), 
                      "est_HR_CIu" = c(NA, NA),  
                      "window" = rep(window_size,2), # Size of the imputation grid
                      "confirmation" = rep(confirmation_window,2), # Time needed to confirm a progression (months)
                      "nobs_imputed" = c(NA, NA), # Total number of imputed values
                      "rmse" = c(NA, NA))
    return(out)
  }
  
  # Identify number of imputed values and their RMSE for imputing the expected value
  
  q_imp_pmm <- evaluate_rmse(impdata$data, y_label = "y_imp_lme3l_mi_pmm", window_size = window_size)
  q_imp_cnv <- evaluate_rmse(impdata$data, y_label = "y_imp_lme3l_mi_conv", window_size = window_size)
  
  y_imp_pmm <-  colnames(impdata$data)[which(str_detect(colnames(impdata$data),"y_imp_lme3l_mi_pmm_"))]
  y_imp_cnv <-  colnames(impdata$data)[which(str_detect(colnames(impdata$data),"y_imp_lme3l_mi_conv_"))]
  
  
  results_pmm <- results_cnv <- data.frame(est_logHR = numeric(), est_HR = numeric(), est_HR_CIl = numeric(), est_HR_CIu = numeric())
  
  
  for (imp in seq(length(y_imp_pmm))) {
    
    # Predictive mean matching analyses
    prog_data <- derive_cdp(data = impdata$data, 
                            outcome_var = y_imp_pmm[imp],
                            window_size = window_size,
                            confirmation_window = confirmation_window,
                            extrapolate = extrapolate_cdp) 
    
    fit <- ipwCoxCluster(data.frame(prog_data), 
                         indID = "centerid", 
                         indA = "x", 
                         indX = c("age", "edss_baseval"), 
                         indStatus = "event", 
                         indTime = "tte", 
                         ties = "breslow",
                         confidence = 0.95)
    
    results_pmm[imp,] <- fit["conventional weights", c("log HR Estimate", "HR Estimate", "HR 95% CI-low", "HR 95% CI-up")] 
    
    # Conversion based analyses
    prog_data <- derive_cdp(data = impdata$data, 
                            outcome_var = y_imp_cnv[imp],
                            window_size = window_size,
                            confirmation_window = confirmation_window,
                            extrapolate = extrapolate_cdp) 
    
    fit <- ipwCoxCluster(data.frame(prog_data), 
                         indID = "centerid", 
                         indA = "x", 
                         indX = c("age", "edss_baseval"), 
                         indStatus = "event", 
                         indTime = "tte", 
                         ties = "breslow",
                         confidence = 0.95)
    
    results_cnv[imp,] <- fit["conventional weights", c("log HR Estimate", "HR Estimate", "HR 95% CI-low", "HR 95% CI-up")] 
    
  }
  
  df_pmm <- data.frame("method" = paste(method_lbl, "(PMM)") ,
                       "est_logHR" = mean(results_pmm$est_logHR), # Estimated logHR
                       "est_HR" = mean(results_pmm$est_HR), # Estimated logHR
                       "est_HR_CIl" = mean(results_pmm$est_HR_CIl), 
                       "est_HR_CIu" = mean(results_pmm$est_HR_CIu),  
                       "window" = window_size, # Size of the imputation grid
                       "confirmation" = confirmation_window, # Time needed to confirm a progression (months)
                       "nobs_imputed" = q_imp_pmm$nobs_imputed, # Total number of imputed values
                       "rmse" = q_imp_pmm$rmse)
  df_cnv <- data.frame("method" = paste(method_lbl, "(EDSS conversion)") ,
                       "est_logHR" = mean(results_cnv$est_logHR), # Estimated logHR
                       "est_HR" = mean(results_cnv$est_HR), # Estimated logHR
                       "est_HR_CIl" = mean(results_cnv$est_HR_CIl), 
                       "est_HR_CIu" = mean(results_cnv$est_HR_CIu),  
                       "window" = window_size, # Size of the imputation grid
                       "confirmation" = confirmation_window, # Time needed to confirm a progression (months)
                       "nobs_imputed" = q_imp_cnv$nobs_imputed, # Total number of imputed values
                       "rmse" = q_imp_cnv$rmse)
  
  return(rbind(df_pmm, df_cnv))
}

run_sim <- function(simpars, 
                    censorFUN,
                    window_size = 3, 
                    confirmation_window = 1, 
                    extrapolate_cdp = "progression",
                    seed) {
  
  
  
  if (!missing(seed)) {
    set.seed(seed)
  } else (
    seed <- NA
  )
  
  # Generate a dataset
  dat <- sim_data(simpars)
  
  # Introduce missing values
  misdat <- censorFUN(dat, keep_baseline_visit = TRUE, remove_missing_visits = FALSE)

  # Save original values for y as reference
  misdat$edss_orig <- dat$edss 
  
  # Derive total follow-up for each patient
  misdat <- derive_max_fup(data = misdat)
  
  #########################################################################################
  # Evaluate the methods for a window of 3 months
  #########################################################################################
  results <- data.frame(method = character(),
                        est_logHR = numeric(), 
                        est_HR = numeric(),
                        est_HR_CIl = numeric(),
                        est_HR_CIu = numeric(),
                        window = numeric(), 
                        confirmation = numeric(),
                        nobs_imputed = numeric(),
                        rmse = numeric()) 
  
  # Reference (no missing data)
  results <- results %>% add_row(evaluate_reference(data = misdat, 
                                         window_size = window_size, 
                                         confirmation_window = confirmation_window))
  
  # LOCF
  results <- results %>% add_row(evaluate_locf(data = misdat, 
                                               window_size = window_size, 
                                               confirmation_window = confirmation_window, 
                                               extrapolate_cdp = extrapolate_cdp,
                                               method_lbl = "LOCF"))
  
  # Rounding
  results <- results %>% add_row(evaluate_rounding(data = misdat, 
                                               window_size = window_size, 
                                               confirmation_window = confirmation_window, 
                                               extrapolate_cdp = extrapolate_cdp,
                                               method_lbl = "Rounding"))
  
  # Multilevel imputation
  results <- results %>% add_row(evaluate_lme_3l(data = misdat, 
                  n.imp = 10, 
                  window_size = window_size, 
                  confirmation_window = confirmation_window, 
                  extrapolate_cdp = extrapolate_cdp,
                  sampling_resid = "condMVN", 
                  method_lbl = "LME-CE"))
  
 
  # Get system info
  sys_info <- Sys.info()
  
  # Output the results
  data.frame("delta_xt" = simpars$delta_xt,
             "sim_id" = seed, 
             "system_name" = sys_info["sysname"],
             "system_machine" = sys_info["machine"],
             "R_version" = R.version$version.string,
             results, 
             row.names = NULL)
}



derive_cdp  <- function(data, 
                        outcome_var = "edss",
                        window_size = 3, 
                        confirmation_window = 1, #how many windows ahead do we need to confirm?
                        extrapolate = "progression") #Options: none, progression, all
  {
  
  data$y <- data[,outcome_var]
  
  if (extrapolate == "all") {
      progdat <- data %>%
        group_by(patid) %>%
        summarize(
          centerid = unique(centerid),
          patid = unique(patid),
          x = unique(trt),
          age = unique(age),
          edss_baseval = first(y),
          min_increase = 0,
          max_fup = max(time),
          tte = max(time),
          event = 0)
  } else if (extrapolate %in% c("none", "progression")) {
    progdat <- subset(data, !is.na(edss)) %>%
      group_by(patid) %>%
      summarize(
        centerid = unique(centerid),
        patid = unique(patid),
        x = unique(trt),
        age = unique(age),
        edss_baseval = first(y),
        min_increase = 0,
        max_fup = max(time),
        tte = max(time),
        event = 0)
  }
  
  # Count number of patients
  n_total <- nrow(progdat)
  
  # Get rid of all (imp) observations that are not on the desired grid
  data <- subset(data, time %% window_size == 0)
  
  # For each method, look at the CDSS at time 't' and compare it to baseline
  ytimes <- unique(data$time)
  
  # Identify column last grid position for each patient
  a_maxfup <- matrix(progdat$max_fup, nrow = n_total, ncol = length(ytimes), byrow = TRUE)
  a_ytimes <- matrix(ytimes, nrow = n_total, ncol = length(ytimes), byrow = TRUE)
  a_tdiff <- a_maxfup - a_ytimes
  a_tdiff[which(a_tdiff < 0)] <- NA
  maxfup_index <- apply( a_tdiff, 1, which.min)
  #maxfup_grid <- ytimes[maxfup_index]
  
  npat <- length(unique(data$patid))
  
  # edss_baseval >= 1 & edss_baseval < 6
  patid_low <- data[which(data[,"time"] == 0 & data[,outcome_var] < 1), "patid"] # Minimum increase 1.5
  patid_med <- data[which(data[,"time"] == 0 & data[,outcome_var] >= 1 & data[,outcome_var] < 6), "patid"] # Minimum increase 1.0
  patid_hig <- data[which(data[,"time"] == 0 & data[,outcome_var] >= 6), "patid"]  # Minimum increase 0.5
  
  # Evaluate each subset based on the baseline value
  data$min_increase <- NA 
  data$min_increase[data$patid %in% patid_low] <- 1.5
  data$min_increase[data$patid %in% patid_med] <- 1.0
  data$min_increase[data$patid %in% patid_hig] <- 0.5
  
  pdat <- matrix(data[,outcome_var], nrow = npat, ncol = length(ytimes), byrow = TRUE)
  refdat <- matrix(subset(data, time == 0)[outcome_var][,1] , nrow = npat, ncol = length(ytimes), byrow = FALSE)
  minchange <- matrix(subset(data, time == 0)$min_increase, nrow = npat, ncol = length(ytimes), byrow = FALSE)
  
  # Minimum change depends on baseline EDSS which is always observed
  progression <- pdat - refdat >= minchange
  
  #edss_progression <- rep(NA, npat) # Store the EDSS scores at progression time
  
  ############################################################################################################
  # Calculate CDP directly for all methods
  ############################################################################################################
  cdp <- array(NA, dim = c(npat, length(ytimes)))
  cdp[,1] <- FALSE # We can never have confirmed progression at baseline
  for (j in 2:(length(ytimes) - confirmation_window)) {
    cdp[,j] <- progression[,j] & progression[,(j + confirmation_window)]
  }
  
  # identify the first TRUE statement for cdp
  first_max <- (apply(cdp, 1, which.max))
  
  # Omit CDP for patients where all rows are FALSE
  rows_no_cdp <- which(rowSums(cdp, na.rm = T) == 0)
  
  if (length(rows_no_cdp) > 0) {
    first_max[rows_no_cdp] <- NA
  }
  
  # Omit CDP that occur after the last follow-up
  mat <- data.frame(col_cdp = first_max, col_maxfup = maxfup_index)
  
  if (extrapolate %in% c("all", "progression")) {
    rows_invalid_cdp <- which(mat$col_cdp > mat$col_maxfup)
  } else if (extrapolate == "none") {
    rows_invalid_cdp <- which(mat$col_cdp >= mat$col_maxfup) # The confirmationm point needs to be at the follow-up point
  }
  
  
  if (length(rows_invalid_cdp > 0)) {
    first_max[rows_invalid_cdp] <- NA
  }
  
  array_ind <- data.frame(patid = unique(data$patid), col_cdp = first_max)
  array_ind <- array_ind[-which(is.na(array_ind$col_cdp)),]
  array_ind$tte_cdp <- ytimes[array_ind$col_cdp]
  
  progdat[which(progdat$patid %in% patid_low), "min_increase"] <- 1.5
  progdat[which(progdat$patid %in% patid_med), "min_increase"] <- 1.0
  progdat[which(progdat$patid %in% patid_hig), "min_increase"] <- 0.5
  
  progdat[which(progdat$patid %in% array_ind$patid), "tte"] <- array_ind$tte_cdp
  progdat[which(progdat$patid %in% array_ind$patid), "event"] <- 1
  
  
  progdat
}