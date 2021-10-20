model{
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
      # sum richness
    }
      rich[site] <- sum(z[site,])
      #linear predictor richness
      log(alpha[site]) <- inprod(
        beta_alpha,
        design_matrix_alpha[site,]
      )
      # code up log likelihood
      alpha_lik[site] <- -alpha[site] + 
        rich[site]*log(alpha[site]) - logfact(rich[site])
      ones[site] ~ dbern(exp(alpha_lik[site])/CONSTANT)
  }
  for(psii in 1:npar_psi){
    # community mu & sd
    beta_psi_mu[psii] ~ dlogis(0,1)
    tau_psi[psii] ~ dgamma(0.001,0.001)
    sd_psi[psii] <- 1 / sqrt(tau_psi[psii])
    beta_rho_mu[psii] ~ dlogis(0,1)
    # Species specific coefficients
    for(species in 1:nspecies){
      beta_psi[species,psii] ~ dnorm(
        beta_psi_mu[psii],
        tau_psi[psii]
      )
    }
  }
  # Alpha diversity
  for(alphai in 1:npar_alpha){
    beta_alpha[alphai] ~ dnorm(0,0.01)
  }
  
}