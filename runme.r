source("R/dgm.r")

# Set up the data-generating mechanism
simpars <- setup() 

# Generate a complete dataset
dat <- sim_data(simpars)

# Introduce missing visits
misdat <- censor_visits_7(data = dat, outcome = "y_dr")
