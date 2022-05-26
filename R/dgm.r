# Version Dec 22, 15h59
library(dplyr, warn.conflicts = FALSE)
#library(Matrix)
#library(sparseMVN)
#library(stringr)
library(truncnorm)
library(MASS)
library(mvtnorm)
library(nlme)

options(dplyr.summarise.inform = FALSE)

#Generate  EDSS outcomes without any noise (simplified DGM by Gabrielle)
sim_data <- function(simpars) {
  
  # Identify total number of patients
  n_total <- simpars$ncenters * simpars$npatients
  
  # Create the grid
  ytimes <- seq(from = 0, to = simpars$follow_up, by = 1)
  
  # Construct Sigma
  cs1 <- simpars$corFUN(value = simpars$rho)
  cor_matrix <- corMatrix(cs1, covariate = ytimes)
  sd_alpha <- rep(simpars$sd_a_t, (length(ytimes)))
  sigma_alpha <- cor2cov(sd = sd_alpha, rho = cor_matrix)
  
  # Draw a prognostic factor
  age <- rtruncnorm(n = n_total, 
                    a = simpars$min_age, 
                    b = Inf, 
                    mean = simpars$mean_age, 
                    sd = simpars$sd_age)
  
  # Allocate treatment for all patients
  ptreat <- simpars$tx_alloc_FUN(age = age) 
  xtreat <- rbinom(n = n_total, size = 1, prob = ptreat)
  
  # Identify the centers
  centerid <- rep(c(1:simpars$ncenters), each = simpars$npatients)
  
  # Draw the patient effects
  alpha_ij <- rnorm(n = n_total, mean = 0, sd = simpars$sd_alpha_ij)
  
  # Draw the center effects
  beta_1j <- rep(rnorm(simpars$ncenters, mean = 0, sd = simpars$sd_beta1_j), each = simpars$npatients)
  
  # Draw epsilon
  epsilon_tij_x0 <-  mvrnorm(n = n_total, mu = rep(0, length(ytimes)), Sigma = sigma_alpha)
  epsilon_tij_x1 <-  mvrnorm(n = n_total, mu = rep(0, length(ytimes)), Sigma = sigma_alpha)
  
  # Generate matrices with disease trajectory for eacht patient
  delta_baseline <- matrix(simpars$intercept, nrow = n_total, ncol = length(ytimes))
  
  # Patient-specific baseline risk (constant  over time)
  delta_patient <- matrix(alpha_ij, nrow = n_total, ncol = length(ytimes), byrow = FALSE) 
  
  # Cluster-specific basline risk (constant over time)
  delta_cluster <- matrix(beta_1j, nrow = n_total, ncol = length(ytimes), byrow = FALSE) 
  
  # Time effect is identical for all patients, but varies over time
  delta_time <- matrix(simpars$beta_t * ytimes + simpars$beta_t2 * (ytimes**2), 
                       nrow = n_total, ncol = length(ytimes), byrow = TRUE) 
  
  # Treatment effect for received treatment
  
  delta_x1 <- matrix(simpars$delta_xt * ytimes + simpars$delta_xt2 * (ytimes**2), 
                     nrow = n_total, ncol = length(ytimes), byrow = TRUE)
  
  # Age effect
  delta_age <- matrix(simpars$beta_age * age, nrow = n_total, ncol = length(ytimes), byrow = FALSE) 
  
  latent_y_x0 <- delta_baseline + delta_patient + delta_cluster + delta_time + delta_age + epsilon_tij_x0
  latent_y_x1 <- delta_baseline + delta_patient + delta_cluster + delta_time + delta_x1 + delta_age + epsilon_tij_x1
  dsx <- cbind(l_x0 = as.vector(t(latent_y_x0)), 
               l_x1 = as.vector(t(latent_y_x1)), 
               x0 = rep(xtreat, each = length(ytimes)) == 0, 
               x1 = rep(xtreat, each = length(ytimes)) == 1,
               l_dr = NA)
  dsx[,"l_dr"] <- ifelse(dsx[,"x0"] == 1, dsx[,"l_x0"], dsx[,"l_x1"])
  
  
  mat <- matrix(NA, nrow = (n_total * length(ytimes)), ncol = 13)
  colnames(mat) <- c("centerid", 
                     "patid", 
                     "dgm_alpha_ij", 
                     "dgm_beta_1j", 
                     "x",
                     "age", 
                     "time", "l_x0", "l_x1", "l_dr", "edss_x0", "edss_x1", "edss")
  
  mat[,"centerid"] <- rep(centerid, each = length(ytimes))
  mat[, "patid"] <- rep(seq(n_total), each = length(ytimes))
  mat[, "dgm_alpha_ij"] <- as.vector(t(delta_patient))
  mat[, "dgm_beta_1j"] <- as.vector(t(delta_cluster))
  mat[, "x"] <- rep(xtreat, each = length(ytimes)) # Received treatment
  mat[, "age"] <- rep(age, each = length(ytimes))
  mat[, "time"] <- rep(ytimes, n_total)
  mat[, "l_x0"] <- dsx[,"l_x0"]  # Latent outcome under x=0
  mat[, "l_x1"] <- dsx[,"l_x1"]  # Latent outcome under x=1
  mat[, "l_dr"] <- dsx[,"l_dr"] # Latent outcome under received treatment
  mat[, "edss_x0"] <- convert_to_EDSS_scale(mat[, "l_x0"]) #Observed outcome
  mat[, "edss_x1"] <- convert_to_EDSS_scale(mat[, "l_x1"]) #Observed outcome
  mat[, "edss"] <- convert_to_EDSS_scale(mat[, "l_dr"]) #Observed outcome
  
  data.frame(mat)
}

setup <- function(tx_alloc_FUN = treatment_alloc_confounding, # Function for treatment allocation
                  delta_xt = 0, # DGM - interaction treatment time
                  delta_xt2 = 0) {
  
  # Default scenario
  npatients <- 500
  ncenters <- 20
  follow_up <- 12*5 # Total follow-up (number of months)
  sd_a_t <- 0.5          # DGM - Within-visit variation in EDSS scores
  intercept <- 1.3295    # DGM - Mean baseline EDDS score
  sd_alpha_ij <- 1.46    # DGM - Between-subject variation in baseline EDSS
  sd_beta1_j <- 0.20     # DGM - Between-site variation in baseline EDSS
  
  # Distribution of age
  mean_age <- 42.41
  sd_age <- 10.53
  min_age <- 18
  
  # Prognostic effects
  beta_age <- 0.05 # DGM - prognostic effect of age
  beta_t <- 0.014  # DGM - prognostic effect of time
  beta_t2 <- 0     # DGM - prognostic effect of time squared
  rho <- 0.8       # DGM - autocorrelation of between alpha_tij
  corFUN <- corAR1 # DGM - correlation structure of the latent EDSS scores
  
  # Report the theoretical mean EDSS at baseline
  mean_edss_t0 <- intercept + beta_age * mean_age

  # Report the theoretical standard deviation of the baseline EDSS within each center
  sd_within_edss_t0 <- sqrt(sd_alpha_ij**2 + (beta_age**2) * (sd_age**2) +  sd_a_t**2)

  # Report the theoretical standard deviation of the mean EDSS across centers
  sd_between_edss_t0 <- sqrt(sd_within_edss_t0**2 + sd_beta1_j**2)
  
  list(npatients = npatients, 
       ncenters = ncenters, 
       follow_up = follow_up,
       sd_a_t = sd_a_t, 
       intercept = intercept,
       sd_alpha_ij = sd_alpha_ij,
       sd_beta1_j = sd_beta1_j,
       mean_age = mean_age,
       sd_age = sd_age,
       min_age = min_age,
       beta_age = beta_age,
       beta_t = beta_t,
       beta_t2 = beta_t2,
       delta_xt = delta_xt,
       delta_xt2 = delta_xt2,
       rho = rho,
       corFUN = corFUN,
       tx_alloc_FUN = tx_alloc_FUN)
}

logit <- function(x) { 
  log(x) - log(1 - x)
}

expit <- function(x) { 
  1/(1 + exp(-x))
}


# Exponential decay function for within person correlation
cor_decay_exp <- function(time1, time2, range = 70) {
  exp(-abs(time1 - time2)/range)
}

cor_comp_symm <- function(time1, time2, rho = 0.87) {
  if (length(time1) > 1 && length(time1) != length(time2)) {
    stop("Invalid input for 'time1'")
  }
  cor <- rep(rho, length(time2)) 
  cor[time1 == time2] <- 1
  cor
}

cor_none <- function(time1, time2) {
  cor_comp_symm(time1, time2, rho = 0)
}

# convert correlation matrix to covariance matrix 
cor2cov <- function(sd, rho) {
  if (length(sd) != nrow(rho) & length(sd) != ncol(rho)) {
    stop("Invalid dimensions of 'sigma' and/or 'rho'")
  }
  D <- diag(sd, nrow = length(sd), ncol = length(sd))
  vmat <- D %*% rho %*% D
  vmat
}

treatment_alloc_randomized <- function(age) {
  0.5
}

treatment_alloc_confounding <- function(age) {
  1/(1 + exp(-(0.7 - 0.032*age - 0.0001*(age**2))))
}

censor_visits_1 <- function(data, 
                            outcome_var = "edss", # The outcome variable 
                            time_var = "time") {
  
  ncenters <- length(unique(data$centerid))
  
  u <- rnorm(ncenters, mean = 0, sd = 0.15)
  
  prob_yobs <- expit(-1.94 + u[data$centerid])
  
  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}

censor_visits_2 <- function(data, 
                            outcome_var = "edss", 
                            time_var = "time") {

  prob_yobs <- rep(expit(-1.94), nrow(data))

  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}


# Visits missing according to center and time
censor_visits_3 <- function(data,
                            outcome_var = "edss", 
                            time_var = "time") {
  
  ncenters <- length(unique(data$centerid))
  
  # Draw the center effects
  u <- rnorm(ncenters, mean = 0, sd = 0.15)
  
  prob_yobs <- expit(-2.70 + u[data$centerid] - 0.7 * log(data[,time_var]/24))
  
  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}

censor_visits_4 <- function(data, 
                            outcome_var = "edss", 
                            time_var = "time") {
  
  ncenters <- length(unique(data$centerid))
  
  # Draw the center effects
  u <- rnorm(ncenters, mean = 0, sd = 0.15)
  
  prob_yobs <- expit(-2.31 + u[data$centerid] - 0.5 * log(data[,time_var]/36) - data$x*0.8)
  
  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}


# Visits missing according to center and time
censor_visits_5 <- function(data, 
                            outcome_var = "edss", 
                            time_var = "time") {
  
  ncenters <- length(unique(data$centerid))
  
  # Draw the center effects
  u <- rnorm(ncenters, mean = 0, sd = 0.15)
  
  prob_yobs <- expit(-2.31 + u[data$centerid] - 1.1 * log(data[,time_var]/36) - data$x*0.8)
  
  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}

censor_visits_6 <- function(data,
                            outcome_var = "edss", 
                            time_var = "time") {
  
  ncenters <- length(unique(data$centerid))
  
  # Draw the center effects
  u <- rnorm(ncenters, mean = 0, sd = 0.15)
  
  prob_yobs <- expit(-1.6 + u[data$centerid] - data$x*0.7)

  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}

censor_visits_7 <- function(data, 
                            outcome_var = "edss", 
                            time_var = "time") {
  
  prob_yobs <- rep(0.03, nrow(data)) 
  prob_yobs[data$x == 0 & data[,time_var] %% 6 == 0] <- 0.85
  prob_yobs[data$x == 1 & data[,time_var] %% 9 == 0] <- 0.67

  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}


# 3 vs 9 months schedule
censor_visits_8 <- function(data, 
                            outcome_var = "edss", 
                            time_var = "time") {
  
  prob_yobs <- rep(0.03, nrow(data))
  prob_yobs[data$x == 0 & data[,time_var] %% 3 == 0] <- 0.35
  prob_yobs[data$x == 1 & data[,time_var] %% 9 == 0] <- 0.55
  
  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}

censor_visits_9 <- function(data, 
                            outcome_var = "edss", 
                            time_var = "time") {
  
  # Define MNAR model
  prob_yobs <- expit(-0.5  -  0.5 * data[,outcome_var] + 0.5 * data$x)

  # Set outcome to NA where missing
  data[rbinom(nrow(data), size = 1, prob = prob_yobs) == 0, outcome_var] <- NA
  
  data %>% na.omit()
}


convert_to_EDSS_scale <- function(x) {
  x <- round(x*2)/2
  x[which(x < 0)] <- 0
  x[which(x > 9.5)] <- 9.5
  x
}

