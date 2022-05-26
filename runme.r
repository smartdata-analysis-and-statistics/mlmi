#####################################################
# Initialize libararies
#####################################################
source("R/dgm.r")

# Set up the data-generating mechanism
simpars <- setup() 

# Generate a dataset without missings
dat <- simdata(simpars)

#censorFUN = censor_visits_7,