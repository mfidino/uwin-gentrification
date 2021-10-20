model{
  #------------------------------
  # multi-species occupancy model
  #------------------------------
  for(site in 1:nsite){
    for(species in 1:nspecies){
      # Linear predictor latent state model
      logit(psi[site,species]) <- inprod(
        beta_psi[species,],
        design_matrix_psi[site,]
      )
      # z is a Bernoulli random variable
      z[site,species] ~ dbern(
        psi[site,species]
      )
      # Linear predictor data model
      logit(rho[site,species]) <- inprod(
        beta_rho[species,],
        design_matrix_rho[site,]
      )
      # y is a binomial process
      y[site,species] ~ dbin(
        rho[site,species] * z[site,species],
        nsurvey
      )
    }
  }
  #-------
  # priors
  #-------
  # Multispecies occupancy
  for(psii in 1:npar_psi){
    # community mu & sd
    beta_psi_mu[psii] ~ dlogis(0,1)
    tau_psi[psii] ~ dgamma(0.001,0.001)
    sd_psi[psii] <- 1 / sqrt(tau_psi[psii])
    beta_rho_mu[psii] ~ dlogis(0,1)
    tau_rho[psii] ~ dgamma(0.001,0.001)
    sd_rho[psii] <- 1 / sqrt(tau_rho[psii])
    # Species specific coefficients
    for(species in 1:nspecies){
      beta_psi[species,psii] ~ dnorm(
        beta_psi_mu[psii],
        tau_psi[psii]
      )
      beta_rho[species,psii] ~ dnorm(
        beta_rho_mu[psii],
        tau_rho[psii]
      )
    }
  }
}