#########################################
#
# Fit initial occupancy model
#
# Written by M. Fidino
#
#########################################




library(dplyr)
library(runjags)
library(googledrive)

source("./R/prep_data_occupancy.R")


my_start <- Sys.time()

cl <- parallel::makeCluster(4)

# Note: Takes about a day and a half to run.
m1 <- runjags::run.jags(
  "./JAGS/multi_scale_occupancy.R",
  monitor = c(
    # among-city regression
    "a_among", "tau_among", "b_among",
    # within-city occupancy
    "b_within", "tau_within", "tau_shape", "tau_rate",
    "b_species", "tau_species","b_species_city",
    # auto-logistic term
    "theta_mu", "tau_theta", "theta_shape", "theta_rate",
    "theta_species", "tau_theta_species", "theta",
    # seasonal variation occupancy
    "c_shape_psi", "c_rate_psi", "c_tau_psi",
    # within-city detection
    "c_det", "tau_det", "tau_shape_rho", "tau_rate_rho",
    "c_species_det","tau_species_det", "c_species_city",
    # seasonal variation detection
    "c_shape_rho","c_rate_rho", "c_tau_rho",
    # latent state of whether species is available for sampling in city
    "x"
  ),
  inits = inits,
  data = data_list,
  n.chains = 4,
  adapt =  1000,
  burnin = 100000,
  sample = 30000,
  thin = 3,
  modules = "glm",
  #method = "parallel"
  method = "rjparallel",
  cl = cl
)

parallel::stopCluster(cl)



saveRDS(m1, "./results/occupancy_model_fit.RDS")
# and then move it to the cloud
googledrive::drive_upload(
  "./results/occupancy_model_fit.RDS",
  "~/gentrification_analysis/occupancy_model_fit.rds"
)


my_end <- Sys.time()

# summarise different parts of the model

my_params <- c(
  # among-city regression
  "a_among", "tau_among", "b_among",
  # within-city occupancy
  "b_within", "tau_within", "tau_shape", "tau_rate",
  "b_species", "tau_species","b_species_city",
  # auto-logistic term
  "theta_mu", "tau_theta", "theta_shape", "theta_rate",
  "theta_species", "tau_theta_species", "theta",
  # seasonal variation occupancy
  "c_shape_psi", "c_rate_psi", "c_tau_psi",
  # within-city detection
  "c_det", "tau_det", "tau_shape_rho", "tau_rate_rho",
  "c_species_det","tau_species_det", "c_species_city",
  # seasonal variation detection
  "c_shape_rho","c_rate_rho", "c_tau_rho"
)

my_sum <- vector("list", length = length(my_params))
for(i in 1:length(my_params)){
  my_sum[[i]] <- summary(
    m1,
    vars = my_params[i]
  )
}

saveRDS(my_sum, "./results/occupancy_model_fit_summary.RDS")


googledrive::drive_upload(
  "./results/occupancy_model_fit_summary.RDS",
  "~/gentrification_analysis/occupancy_model_fit_summary.RDS"
)


