library(mice)
library(nlme)
library(Matrix)

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
  imp <-  mice(impdat, predictorMatrix = predM, method = imeth, maxit = maxit, m = n.imp, printFlag = FALSE)
  
  return(imp)
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



### Fast implementation of MICE algorithm
impute_y_mice_3l <- function(data, 
                             n.imp = 10, 
                             sampling_resid = "condMVN",
                             outcome_var = "edss",
                             treat_var = "trt",
                             time_var = "time"
) {
  
  data$time <- data[,time_var]
  data$trt <- data[,treat_var]
  data$y <- data[,outcome_var]
  data$time_trt <- data$time * data$trt
  
  ####################### PREP
  y <- data$y
  yimp_raw <- yimp_pmm <- replicate(n.imp, y)
  ry <- !(is.na(data$y))
  x <- data[,c("patid", "centerid", "time", "time_trt", "age")]
  type <- c("patid" = -3, "centerid" = -2, "time"  = 6, "time_trt" = 1, "age" = 1)
  
  ####################### START
  wy <- !ry
  x <- cbind(1, as.matrix(x))
  type <- c(2, type)
  names(type)[1] <- colnames(x)[1] <- "(Intercept)"
  
  clust <- names(type[type == -3]) # Patient
  rande <- names(type[type == 2])
  group <- names(type[type == -2]) # Center
  time <- names(type[type == 6]) # Time
  
  fixe <- names(type[type > 0])
  lev <- unique(x[, clust])
  X <- x[, fixe, drop = FALSE]
  Z <- x[, rande, drop = FALSE]
  xobs <- x[ry, , drop = FALSE]
  yobs <- y[ry]
  Xobs <- X[ry, , drop = FALSE]
  Zobs <- Z[ry, , drop = FALSE]
  
  fixmodel <- paste("yobs ~ ", paste(fixe[-1L], collapse = "+"))
  
  fit <- try(lme(formula(fixmodel), 
                 random =  formula(paste("~1|", group, "/", clust)), 
                 correlation = corExp(form = formula(paste("~", time, "|", group, "/", clust))), 
                 data = data.frame(yobs, xobs),
                 control = list(returnObject = TRUE)))
  
  
  if (("try-error" %in% class(fit))) {
    warning("Estimation of multilevel model failed!")
    return(list(data = data, fit = NULL))
  }
  
  vnames <- names(fixef(fit))
  
  # Store the observed predictions
  yhatobs <- predict(fit, level = 2)
  
  # Derive standard error for key model parameters
  fit_intervals <- try(intervals(fit))
  sigmahat <- fit$sigma # fit_intervals$sigma["est."]
  df <- nrow(fit$data) - length(fixef(fit))
  
  if (!("try-error" %in% class(fit_intervals))) {
    se_log_tau_center <- (log(((fit_intervals$reStruct[[group]])["sd((Intercept))","upper"])) - log(((fit_intervals$reStruct[[group]])["sd((Intercept))","lower"])))/(2*qnorm(0.975))
    se_log_tau_patient <- (log(((fit_intervals$reStruct[[clust]])["sd((Intercept))","upper"])) - log(((fit_intervals$reStruct[[clust]])["sd((Intercept))","lower"])))/(2*qnorm(0.975))
    
    # Draw random sample for the cluster effects
    tau_center_star <- exp(rnorm(n = n.imp, mean = log((fit_intervals$reStruct[[group]])["sd((Intercept))","est."]), sd = se_log_tau_center))
    
    # Draw random sample for the patient effects
    tau_patient_star <- exp(rnorm(n = n.imp, mean = log((fit_intervals$reStruct[[clust]])["sd((Intercept))","est."]), sd = se_log_tau_patient))
    
    psi_star <- (tau_center_star**2) + (tau_patient_star**2)
    
    # Draw random sample for the range of the autocorrelation matrix
    se_log_range = (log(fit_intervals$corStruct[, "upper"]) - log(fit_intervals$corStruct[, "lower"])) / (2*qnorm(0.975))
    range_star <- exp(rnorm(n.imp, mean = log(fit_intervals$corStruct[, "est."]), sd = se_log_range))
    
  } else {
    warning("Error when estimating the confidence intervals for the multilevel model")
    
    rancoef_group <- as.matrix(ranef(fit)[[group]])
    lambda_group <- t(rancoef_group) %*% rancoef_group
    psi_group_star <- rep(lambda_group, n.imp)/rchisq(n.imp, df = nrow(rancoef_group)) 
    
    rancoef_clust <- as.matrix(ranef(fit)[[clust]])
    lambda_clust <- t(rancoef_clust) %*% rancoef_clust
    psi_clust_star <- rep(lambda_clust, n.imp)/rchisq(n.imp, df = nrow(rancoef_clust)) 
    
    psi_star <- psi_group_star + psi_clust_star
    
    # 
    cs <- fit$modelStruct$corStruct
    range_est <- as.numeric(coef(cs, unconstrained = FALSE))
    range_star <- rep(range_est, n.imp) #Ignore uncertainty in the estimated range
  }
  
  spatDat <- data.frame(time = 0:max(x[,time]))
  
  
  
  pbimp <- txtProgressBar(min = 0, max = n.imp, style = 3)
  for (imp in 1:n.imp) {
    # Draw a random sample for the residual error
    sigma2star <- df * sigmahat^2/rchisq(n = 1, df = df)
    
    # Rescale the covariance matrix to the new draw of sigma
    covmat <- sigma2star * (vcov(fit)/sigmahat^2)
    rv <- t(chol(covmat))
    
    # Draw random sample for the beta coefficients
    beta_star <- fixef(fit) + rv %*% rnorm(ncol(rv))
    rownames(beta_star) <- vnames
    
    
    cs1Exp <- corExp(value = range_star[imp], form = ~ time)
    cs1Exp <- Initialize(cs1Exp, spatDat)
    sigma_full <- cor2cov(sd = rep(sqrt(sigma2star), nrow(spatDat)), rho = corMatrix(cs1Exp))
    
    # Population-level predictions
    blup_pop_beta <- X[,vnames] %*% beta_star[vnames,]
    
    # Patient-level predictions
    blup_clust_beta <- y_pred_resid <- rep(0, nrow(blup_pop_beta))
    
    # Set up sigma for all clusters
    sigma_all_pat <- list()
    zi_psistar_zi <- list()
    psistar_zi <- list()
    resid_pop_blup <- list()
    zi_all_pat <- list()
    
    for (jj in lev) {
      Xi <- as.matrix(Xobs[xobs[, clust] == jj, ])
      if (sum(xobs[,clust] == jj)  <= 1) {
        # If we only have one observation, Xi is a single-column matrix but should be a single-row
        Xi <- t(Xi)
      }
      
      yi <- yobs[xobs[, clust] == jj]
      ti <- as.matrix(Xobs[xobs[, clust] == jj, time])
      Zi <- as.matrix(Zobs[xobs[, clust] == jj, ])
      idx <- spatDat$time %in% ti
      
      zi_all_pat[[jj]] <- Zi
      sigma_all_pat[[jj]] <- sigma_full[idx, idx]
      
      psistar_zi[[jj]] <- psi_star[imp] %*% t(Zi)
      zi_psistar_zi[[jj]] <- Zi %*% psistar_zi[[jj]]
      
      resid_pop_blup[[jj]] <- (yi - Xi[,vnames] %*% beta_star[vnames,])
    }
    
    sigma_all_pat <- do.call(bdiag, sigma_all_pat)
    psistar_zi <- do.call(c, psistar_zi)
    zi_psistar_zi <- do.call(bdiag, zi_psistar_zi)
    resid_pop_blup <- do.call(bdiag, resid_pop_blup)
    zi_all_pat <- do.call(bdiag, zi_all_pat)
    
    Mi <- psistar_zi %*% chol2inv(chol(zi_psistar_zi + sigma_all_pat))
    myi <- Mi %*% resid_pop_blup # Mean of the conditional random effect for each cluster
    vyi <- psi_star[imp] - Mi %*% zi_all_pat * psi_star[imp]
    
    # Draw random draws for the random effects
    bi_star <- myi + sqrt(vyi) * rnorm(length(myi))
    
    # Derive the linear predictor on the cluster level
    blup_clust_beta <- blup_pop_beta + rep(bi_star, each = nrow(spatDat))
    
    y_pred_resid <- y - blup_clust_beta
    y_imp_resid <- y_pred_resid
    
    # Identify observed and missing residuals
    dependent.ind <- which(is.na(y_pred_resid))
    given.ind <- which(!is.na(y_pred_resid))
    
    if (sampling_resid == "condMVN") {
      sigma_all_pat <- do.call(bdiag, replicate(length(lev), sigma_full, simplify = F))
      
      B <- sigma_all_pat[dependent.ind, dependent.ind]
      C <- sigma_all_pat[dependent.ind, given.ind, drop = FALSE]
      D <- sigma_all_pat[given.ind, given.ind]
      
      CDinv <-  C %*% chol2inv(chol(D)) 
      
      # Derive the conditional residual errors
      cMu <- as.vector(CDinv %*% y_pred_resid[given.ind])
      CH_cVar <- Cholesky(B - CDinv %*% t(C))
      
      y_pred_resid[dependent.ind] <- cMu
      y_imp_resid[dependent.ind] <-  sparseMVN::rmvn.sparse(1, mu = cMu, CH = CH_cVar, prec = FALSE)[1,]
    } else if (sampling_resid == "MVN") {
      # Generate a sample for each patient
      cMu <- rep(0, length(spatDat$time))
      random_resid <- rmvnorm(n = length(lev), mean = cMu, sigma = sigma_full)
      
      y_pred_resid[dependent.ind] <- rep(0, length(dependent.ind))
      y_imp_resid[dependent.ind] <-  as.vector(t(random_resid))[dependent.ind]
    } else {
      y_pred_resid[dependent.ind] <- rep(0, length(dependent.ind))
      y_imp_resid[dependent.ind] <-  rnorm(length(dependent.ind)) * sqrt(sigma2star)
    }
    
    ## Store the raw imputations
    yimp_raw[wy, imp] <-  blup_clust_beta[wy] + y_imp_resid[dependent.ind] # We add "random" conditional residuals
    
    ## Apply predictive mean matching
    yhatmis <- blup_clust_beta[wy] + y_pred_resid[wy] #Note: we here use the expected conditional residuals
    idx_pmm <- mice:::matchindex(d = yhatobs, t = yhatmis, k = 5L)
    yimp_pmm[wy, imp] <- y[ry][idx_pmm]
    
    
    # Update progress bar
    setTxtProgressBar(pbimp, imp)
  }
  close(pbimp)
  
  data[, paste0("y_imp_lme3l_mi_conv_", seq(n.imp))] <- convert_to_EDSS_scale(yimp_raw[,seq(n.imp)])
  data$y_imp_lme3l_mi_conv <- convert_to_EDSS_scale(rowMeans(yimp_raw))
  
  data[, paste0("y_imp_lme3l_mi_pmm_", seq(n.imp))] <- yimp_pmm[,seq(n.imp)]
  data$y_imp_lme3l_mi_pmm <- rowMeans(yimp_pmm)
  
  # ONly keep original columns
  data <- data %>% dplyr::select(centerid, patid, treat_var, age, time_var, outcome_var, edss_orig, 
                                 max_fup, 
                                 paste0("y_imp_lme3l_mi_conv_", seq(n.imp)), 
                                 paste0("y_imp_lme3l_mi_pmm_", seq(n.imp)),
                                 y_imp_lme3l_mi_conv,
                                 y_imp_lme3l_mi_pmm)
  
  
  list(data = data, fit = fit)
}

