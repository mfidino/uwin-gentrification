msom_mcmc <- function(seed, data_list, constant_list){
library(nimble)
  
msom_code <-nimble::nimbleCode(
  {
    for(site in 1:nsite){
      for(species in 1:nspecies){
        logit(psi[site,species]) <- inprod(
          beta_psi[species, 1:npar_psi],
          design_matrix_psi[site,1:npar_psi]
        )
        z[site,species] ~ dbern(
          psi[site,species]
        )
        logit(rho[site,species]) <- inprod(
          beta_rho[species,1:npar_rho],
          design_matrix_rho[site,1:npar_rho]
        )
        y[site,species] ~ dbin(
          rho[site,species] * z[site,species],
          nsurvey
        )
      }
    }
    for(psii in 1:npar_psi){
      beta_psi_mu[psii] ~ dlogis(0,1)
      tau_psi[psii] ~ dgamma(0.001,0.001)
      sd_psi[psii] <- 1 / sqrt(tau_psi[psii])
      for(species in 1:nspecies){
        beta_psi[species,psii] ~ dnorm(
          beta_psi_mu[psii],
          tau_psi[psii]
        )
      }
    }
    for(rhoi in 1:npar_rho){
      beta_rho_mu[rhoi] ~ dlogis(0,1)
      tau_rho[rhoi] ~ dgamma(0.001,0.001)
      sd_rho[rhoi] <- 1 / sqrt(tau_rho[rhoi])
      for(species in 1:nspecies){
        beta_rho[species,rhoi] ~ dnorm(
          beta_rho_mu[rhoi],
          tau_rho[rhoi]
        )
      }
    }
  }
)

my_inits <- function(){
  list(
    z = matrix(
      1,
      nrow = constant_list$nsite,
      ncol = constant_list$nspecies
    ),
    beta_psi = matrix(
      rnorm(constant_list$npar_psi * constant_list$nspecies),
      nrow  = constant_list$nspecies,
      ncol = constant_list$npar_psi),
    beta_rho = matrix(
      rnorm(constant_list$npar_psi * constant_list$nspecies),
      nrow  = constant_list$nspecies,
      ncol = constant_list$npar_psi)
  )
}

msom <- nimble::nimbleModel(
  code = msom_code,
  name = "msom",
  constants = constant_list,
  data = data_list,
  inits = my_inits()
)

tmp_msom <- compileNimble(msom)

msom_configure <- nimble::configureMCMC(
  msom
)
msom_configure$addSampler(
  target = c("beta_psi_mu", "beta_psi"),
  type = "RW_block",
  control = list(adaptInterval = 100)
)
msom_configure$addSampler(
  target = c("beta_rho_mu", "beta_rho"),
  type = "RW_block",
  control = list(adaptInterval = 100)
)

msom_configure$setMonitors(
  c("beta_psi_mu", "beta_psi", "beta_rho_mu","beta_rho", "sd_psi", "sd_rho")
)

msom_configure$setMonitors2("z")
msom_configure$setThin2(5)

msom_build <- nimble::buildMCMC(
  msom_configure
)

# compile the model
msom_compiled <- nimble::compileNimble(
  msom_build,
  project = msom,
  resetFunctions = TRUE
)



results <- nimble::runMCMC(
  msom_compiled,
  niter = 40000,
  setSeed = seed,
  inits = my_inits()
)

return(results)
}

cli::cli_alert_info("Model code created as object 'msom_compiled'")
