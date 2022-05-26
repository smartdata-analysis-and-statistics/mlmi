source("R/dgm.r")

# Set up the data-generating mechanism
simpars <- setup() 

# Generate a complete dataset
dat <- sim_data(simpars)

# Introduce informative visit patterns
misdat <- censor_visits_7(data = dat)

# Visualize trajectories
