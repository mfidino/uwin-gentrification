#####################################################
#
# Simulate latent state from occupancy model results
#
#
#
####################################################

analysis <- "beta"

# placeholder for now until I run the model without cougar
species_to_drop <- "cougar"

library(runjags)
library(dplyr)

# Prep the data for the model
source("./R/prep_data_occupancy.R")

source("./R/alpha_beta_functions.R")

cat("loading in run.jags file...\n")
# Load in the occupancy model results
mout <- readRDS(
  "./results/occupancy_model_fit_simpler2.RDS"
)

cat("binding posterior simulations...\n")
# Compile the posterior
mcmc <- do.call(
  "rbind",
  mout$mcmc
)

# take a random sample to iterate through
set.seed(11556644)
my_samples <- ifelse(analysis == "beta", 5000, 10000)
mcsamp <- mcmc[sample(1:nrow(mcmc), my_samples),]

rm(mout, mcmc)
gc()
# make pieces of this sample because we cannot iterate through
#  the whole thing
mcsamp_list <- vector(
  "list",
  length = 10
)

my_groups <- rep(1:10, each = my_samples/10)
ngroup <- 10
for(i in 1:10){
  mcsamp_list[[i]] <- split_mcmc(
    mcsamp[which(my_groups == i),]
  )
}
# and then just get all of it too for some other
#  calculations
mcsamp <- split_mcmc(
  mcsamp
)

# simulate the probability of:
#  1) Species presence in a city
#  2) Species occupancy at a site
#  3) Species probability of detection

# determine if species was detected in a city. This will
#  modify our z matrix posterior calculation.
sp_observed <- apply(
  mcsamp$x,
  c(2,3),
  sum
)
sp_observed[sp_observed < nrow(mcsamp$x)] <- 0
sp_observed[sp_observed > 0] <- 1

# the z vector for one mcmc sample, plus the other objects
#  that are of similar length. The object tmp_covs has all
#  the site info and is the same length. 
site_info <- tmp_covs

# figure out how many unique site, city, seasons we have
unq_site_samps <- site_info[,1:3]
unq_site_samps <- unq_site_samps[!duplicated(unq_site_samps),]

unq_site_samps$Season <- order_seasons(unq_site_samps$Season)

# order by city, season, and then site
unq_site_samps <- unq_site_samps[
  order(unq_site_samps$City, unq_site_samps$Season, unq_site_samps$Site),
]

# and now make a matrix to store the sp_rich results for
if(analysis == "alpha"){
  sp_rich_mcmc <- matrix(
    NA,
    ncol = nrow(unq_site_samps),
    nrow = nrow(mcsamp$a_among)
  )
  sp_rich_mcmc <- matrix(
    NA,
    ncol = 999,
    nrow = nrow(mcsamp$a_among)
  )
  
}
if(analysis == "beta"){
  sp_dat <- tmp
  
  beta_results <- vector("list", length = my_samples)
}

# and now we need to split the mcmc into pieces because we cannot store
#  all of the results at once. Splitting into pieces of 1K.
data_list$ncov_within <- 4
data_list$ncov_det <- 4

data_list$psi_covs <- cbind(
  data_list$psi_covs,
  data_list$psi_covs[,2] * data_list$psi_covs[,3]
)
data_list$rho_covs <- cbind(
  data_list$rho_covs,
  data_list$rho_covs[,2] * data_list$rho_covs[,3]
)
for(gr in 1:ngroup){
  cat(
    paste0("\ngroup ", gr, " of ", ngroup,"...\n")
  )
  mcsamp_piece <- mcsamp_list[[gr]]
  # where to store stuff in sp_rich_mcmc
  mc_loc <- which(my_groups == gr)
  
  source("./R/sample_z.R")
}


# calculate mean and sd
if(analysis == "alpha"){
  sp_rich$mu <- apply(sp_rich_mcmc, 2, mean)
  sp_rich$sd <- apply(sp_rich_mcmc, 2, sd)
  sp_rich <- sp_rich[,-which(colnames(sp_rich) == "rich")]
  write.csv(
    sp_rich,
    "./results/alpha_for_stage_two_collapsed.csv",
    row.names = FALSE
  )
}

if(analysis == "beta"){
  gc()
  pb <- txtProgressBar(max = length(beta_results))
  for(i in 1:length(beta_results)){
    setTxtProgressBar(pb, i)
    beta_results[[i]] <- data.frame(
      dissim = beta_results[[i]][,1],
      rich = beta_results[[i]][,2],
      loc = as.character(1:nrow(beta_results[[i]]))
    )
  }
  beta_results <- do.call("rbind", beta_results)
  gc()
  beta_summary <- beta_results %>% 
    dplyr::group_by(loc) %>% 
    dplyr::summarise(
      mu_beta = mean(dissim, na.rm = TRUE),
      var_beta = sd(dissim, na.rm = TRUE)^2,
      mu_rich = mean(rich, na.rm = TRUE),
      var_rich = sd(rich, na.rm = TRUE)^2,
      na_count = sum(is.na(dissim))
    )
  rm(beta_results)
  gc()
  
  write.csv(
    beta_summary,
    "./results/beta_summary_for_analysis_collapsed_vegan.csv",
    row.names = FALSE
  )

}
