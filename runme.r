source("R/dgm.r")
source("R/vis.r")
source("R/imp.r")

# Set up the data-generating mechanism
# Introduce a modest treatment effect of 0.007 change in EDSS per day
simpars <- setup(delta_xt = -0.007) 

# Generate a complete dataset
dat <- sim_data(simpars)

# Introduce informative visit patterns
misdat <- censor_visits_7(dat)

# Visualize the observed trajectory of one patient
plot_example_trajectory(misdat, sel_patid = 1)

# Visualize the age distribution for both treatment groups
plot_dens_x(misdat, x_var = "age", x_label = "Age")

# Visualize the total follow-up for each center
plot_max_fup(misdat)

imptimes <- seq(0, 60, by = 3) # Which time points do we need imputations?
imp <- impute(misdat, times = imptimes, maxit = 1)

# Visualize the imputed trajectory of one patient