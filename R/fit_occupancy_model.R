#########################################
#
# Fit initial occupancy model
#
# Written by M. Fidino
#
#########################################




library(dplyr)
library(runjags)
#library(googledrive)

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
  burnin = 125000,
  sample = 30000,
  thin = 3,
  modules = "glm",
  #method = "parallel"
  method = "rjparallel",
  cl = cl
)

parallel::stopCluster(cl)



saveRDS(m1, "./results/occupancy_model_fit_simpler2.RDS")
# and then move it to the cloud
#googledrive::drive_upload(
#  "./results/occupancy_model_fit.RDS",
#  "~/gentrification_analysis/occupancy_model_fit.rds"
#)

#m1 <- readRDS("./results/occupancy_model_fit_simpler.RDS")

my_end <- Sys.time()

# summarise different parts of the model

my_pars <- colnames(m1$mcmc[[1]])

# split into batches of 100
my_pars <- split(
  my_pars,
  factor(floor(1:length(my_pars) / 100))
)
summary_list <- vector("list", length = length(my_pars))

for(i in 1:length(my_pars)){
  print(i)
  summary_list[[i]] <- summary(m1, vars = my_pars[[i]])
}

summary_list <- do.call("rbind", summary_list)


saveRDS(summary_list, "./results/occupancy_model_fit_simpler_summary2.RDS")

# sample z


my_start <- Sys.time()

cl <- parallel::makeCluster(4)



mz <- runjags::extend.jags(
  m1,
  add.monitor = "z",
  drop.monitor = c(
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
  modules = "glm",
  adapt =  1000,
  sample = 5000,
  thin = 3,
  method = "rjparallel",
  cl = cl
)

parallel::stopCluster(cl)





cl <- parallel::makeCluster(4)

# Note: Takes about a day and a half to run.
mz <- runjags::run.jags(
  "./JAGS/multi_scale_occupancy.R",
  monitor = c(
    "z"
  ),
  inits = inits,
  data = data_list,
  n.chains = 4,
  adapt =  1000,
  burnin = 100000,
  sample = 5000,
  thin = 3,
  modules = "glm",
  #method = "parallel"
  method = "rjparallel",
  cl = cl
)

parallel::stopCluster(cl)


saveRDS(mz, "./results/occupancy_z.RDS")