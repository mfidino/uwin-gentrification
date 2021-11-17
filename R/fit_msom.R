library(sf)
library(runjags)
library(dplyr)
library(cli)


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
              "tau_species_theta","z"),
  data = jags_list,
  n.chains = 3,
  burnin = 50,
  adapt = 50,
  sample = 50,
  method = "parallel",
  inits = jags_inits
)

# save it
saveRDS(m1, "./mcmc_output/trial_msom.RDS")

