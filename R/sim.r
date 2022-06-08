require(mice)
require(nlme)
require(ipwCoxCSV)

derive_max_fup <- function(data,
                           outcome_var = "edss") {
  
  data$y <- data[,outcome_var]
  
  data <- data %>%
    group_by(patid) %>%
    mutate(max_fup = max(time[!is.na(y)]))
  
  data <- dplyr::select(data, -y)
  
  data.frame(data)
}

impute_y_locf <- function(data,
                          outcome_var = "edss") {
  
  ytimes <- 0:max(data$time)
  y_cens <- array(NA, dim = c(length(unique(data$patid)), length(ytimes)))
  
  for (i in 1:length(ytimes)) {
    y_cens[,i] <- subset(data, time == ytimes[i])[,outcome_var]
    
    if (i > 1) {
      # Replace missing values by observations from previous time
      y_cens[is.na(y_cens[,i]),i] <- y_cens[is.na(y_cens[,i]),(i - 1)] 
    }
  }
  
  # Fill in imputations
  data$y_imp_locf <- as.vector(t(y_cens))
  
  data
}

impute_y_rounding <- function(data,
                              outcome_var = "edss") {

  ytimes <- 0:max(data$time)
  y_cens <- array(NA, dim = c(length(unique(data$patid)), length(ytimes)))
  
  for (i in 1:length(ytimes)) {
    y_cens[,i] <- subset(data, time == ytimes[i])[,outcome_var]
  }
  
  y_imputed <- y_cens
  
  y_obs <- array(0, dim = dim(y_cens))
  y_obs[!is.na(y_cens)] <- 1
  
  
  to_impute <- which(is.na(y_cens), arr.ind = TRUE)
  pb <- txtProgressBar(min = 0, max = nrow(to_impute), style = 3)
  for (i in 1:nrow(to_impute)) {
    # Retrieve column and row
    col <- to_impute[i, "col"]
    row <- to_impute[i, "row"]
    
    # Identify imputation sources
    imp_source <- which(y_obs[row,] == 1)
    
    # Identify imputation distance
    imp_cols <- abs(col - imp_source)
    
    # Impute!
    y_imputed[row, col] <- mean(y_cens[row,imp_source[which(imp_cols == min(imp_cols))]])
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Fill in imputations
  data$y_imp_rounding <- as.vector(t(y_imputed))
  
  data
}


evaluate_rmse <- function(data,
                          outcome_var = "edss", 
                          orig_outcome_var = "edss_orig", 
                          y_label, 
                          window_size = 1) {
  
  data <- rename(data, y = outcome_var)
  
  # Restrict the dataset to the follow-up and grid
  impdata <- subset(data, time <= max_fup & time %% window_size == 0 & is.na(y))
  
  rmse <- sqrt(mean((impdata[,y_label] - impdata[,orig_outcome_var])**2))
  data.frame(nobs_imputed = nrow(impdata), rmse = rmse)
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
  
  results <- results %>% add_row(evaluate_reference(data = misdat, 
                                         window_size = window_size, 
                                         confirmation_window = confirmation_window))
  
  results <- results %>% add_row(evaluate_locf(data = misdat, 
                                               window_size = window_size, 
                                               confirmation_window = confirmation_window, 
                                               extrapolate_cdp = extrapolate_cdp,
                                               method_lbl = "LOCF"))
  
  results <- results %>% add_row(evaluate_rounding(data = misdat, 
                                               window_size = window_size, 
                                               confirmation_window = confirmation_window, 
                                               extrapolate_cdp = extrapolate_cdp,
                                               method_lbl = "Rounding"))
  
 
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
  data.frame("method" = "Reference (OBS)",
             "est_logHR" = fit["conventional weights","log HR Estimate"], # Estimated logHR
             "est_HR" = fit["conventional weights","HR Estimate"], # Estimated HR
             "est_HR_CIl" = fit["conventional weights","HR 95% CI-low"], 
             "est_HR_CIu" = fit["conventional weights","HR 95% CI-up"],  
             "window" = window_size,
             "confirmation" = confirmation_window,
             "nobs_imputed" = 0,
             "rmse" = 0)
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