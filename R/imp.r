library(mice)

source("R/mlmi.r")

impute <- function(dat, 
                   times = seq(0, 60, by = 3),
                   n.imp = 10,
                   maxit = 5,
                   outcome_var = "edss",
                   treat_var = "trt",
                   time_var = "time") {
  
  dat$time <- dat[,time_var]
  dat$trt <- dat[,treat_var]
  
  impdat <- data.frame(patid = rep(unique(dat$patid), each = length(times)),
                       time = rep(times, length(unique(dat$patid))))
  
  impdat <- bind_rows(dat, impdat)

  impdat <- impdat %>% group_by(patid) %>% arrange(time, .by_group = TRUE) %>% summarize(
    centerid = first(na.omit(centerid)),
    patid = first(patid),
    time = unique(time),
    trt = first(na.omit(trt)),
    age = first(na.omit(age)))
  
  # add EDSS scores
  impdat <- full_join(impdat, dat, by = c("patid" = "patid", "time" = "time")) %>% 
    dplyr::select(centerid.x, patid, trt.x, age.x, time, edss)
  colnames(impdat) <- c("centerid", "patid", "trt", "age", "time", "edss")
  
  # Add interaction term between time and treatment
  impdat$trttime <- impdat$time * impdat$trt
  
  # Prepare imputation method
  imp0 <- mice(impdat, maxit = 0)

  #imputation Method
  imeth <- imp0$method
  imeth[outcome_var] <- "mlmi"
  
  predM <- imp0$predictorMatrix
  predM[outcome_var, "centerid"] <- -2
  predM[outcome_var, "patid"] <- -3
  predM[outcome_var, time_var] <- 6
  diag(predM) <- 0
  
  # Run the imputation
  imp <-  mice(impdat, predictorMatrix = predM, method = imeth, maxit = maxit, m = n.imp)
  
  return(imp)
}

