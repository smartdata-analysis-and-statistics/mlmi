#####################################################
# Initialize libararies
#####################################################
source("R/dgm.r")

#####################################################
# Initialize background processes
#####################################################
logger <- flog.namespace() # Initiate Logger

#####################################################
# Generate a dataset
#####################################################

simpars <- load_scenario(censorFUN = censor_visits_7,
                         tx_alloc_FUN = treatment_alloc_confounding,
                         logger = logger) 
dat <- simdata(simpars, logger)
