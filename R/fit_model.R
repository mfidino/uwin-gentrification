library(sf)
library(nimble)
library(dplyr)
library(cli)
library(MCMCvis)

# the path to the data
data_path <- "../uwin-dataset/data_2021-07-27/cleaned_data/full_capture_history.csv"

# ALL OF THESE WILL BE REMOVED AT SOME POINT. 
# The city abbreviations
cities <- c(
  "ahga", "autx", "boma","buny","chil","deco","inin","ioio"#,
  #"jams","lrar","lbca","mawi","naca","phaz","phaz2","poor",
  #"rony","scut","sewa","slmo","tawa","toon","uril","wide"
)
# The spceies
species <- c(
  "cottontail_sp", "coyote", "fox_squirrel",# "gray_fox",
  "gray_squirrel_sp"#, "mule_deer", "raccoon", "red_fox",
  #"striped_skunk", "virginia_opossum", "white_tailed_deer", "woodchuck"
)

# the years of data to collect (using for now)
years <- c(
  "18|19|20|21"
)

# Test the autologistic indexing term is specified right. Works on
#  current dataset, but want to be able to test it with new data.
test_autologistic <- FALSE

# Prep data for the model
source("./R/format_data_for_analysis.R")

# Get the model ready
source("./nimble/msom.R")
# Make the model

# Currently setting this up for 1 chain to test,
#  will abstract out to running in parallel later...
model <- nimble::nimbleModel(
  msom_code,
  data = data_list,
  constants = constant_list,
  inits = inits(1)
)

to_monitor <- names(my_inits())
to_monitor <- to_monitor[
  -which(
  to_monitor %in% c("z", "x", "ssc_psi", "ssc_rho")
  )
]
mcmc_out <- nimble::nimbleMCMC(
  msom_code,
  constants = constant_list,
  data = data_list,
  inits = my_inits(),
  dimensions = list(psi = constant_list$ndata),
  nchains = 1,
  niter = 1000,
  nburnin = 500,
  monitors = to_monitor 
)

model_configured <- nimble::configureMCMC(
  model,
  print = TRUE,
  monitors = to_monitor
)

model_configured$addMonitors("b_species")
model_mcmc <- nimble::buildMCMC(
  model_configured
)

model_compiled <- nimble::compileNimble(
  model_mcmc,
  project = model
)

results <- nimble::runMCMC(
  model_compiled,
  niter = 50
)


dl <- constant_list

dl$y <- data_list$y
dl$J <- data_list$J
dl$psicov <- data_list$psicov
dl$omegacov <- data_list$omegacov
dl$rhocov <- data_list$rhocov

ini <- my_inits()
m1 <- run.jags(
  "./jags/msom.R",
  monitor = c("a_species", "b_species_city", "theta_psi",
              "c_species_city", "a_community", "tau_community_omega",
              "b_community", "tau_community_psi", "tau_species_psi",
              "tau_shape_psi", "tau_rate_psi", "b_species", 
              "c_community", "tau_community_rho", "tau_shape_rho",
              "tau_rate_rho", "c_species", "tau_species_rho",
              "city_shape_psi", "city_rate_psi", "city_shape_rho",
              "city_rate_rho", "city_tau_psi", "city_tau_rho",
              "theta_community", "tau_community_theta",
              "tau_shape_theta", "tau_rate_theta", "theta_species",
              "tau_species_theta"),
  data = dl,
  n.chains = 3,
  burnin = 20000,
  adapt = 1000,
  sample = 10000,
  method = "parallel",
  inits = inits_longshot
)

msum <- summary(m1)
  
yo <- round(msum,2)
