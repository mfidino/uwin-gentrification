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
  #----------------------
  # Alpha diversity model
  #----------------------
  # for(site in 1:nsite){
  #   # Derive species richness
  #   rich[site] <- sum(z[site,])
  #   # linear predictor
  #   log(mu_alpha[site]) <- inprod(
  #     beta_alpha,
  #     design_matrix_alpha[site,]
  #     )
  #   # Poisson likelihood
  #   alpha_lik[site] <- -mu_alpha[site] + 
  #     rich[site]*log(mu_alpha[site]) - logfact(rich[site])
  #   # rich is a Poisson random variable
  #   alpha_ones[site] ~ dpois(exp(alpha_lik[site])/CONSTANT)
  # }
  #---------------------
  # Beta diversity model
  #---------------------
  # for(i in 1:n){
  #   # Get number of dissimilar species between site pairs
  #   y1[i] <- sum(
  #     (1 - z[siteA_id[i],]) *
  #       z[siteB_id[i],]
  #   )
  #   # Get total richness between site pairs
  #   y2[i] <- nspecies - sum(
  #     (1 - z[siteA_id[i],]) *
  #       (1 - z[siteB_id[i],])
  #   )
  #   # Linear predictor
  #   logit(pi[i]) <- b0 + inprod(
  #     beta_beta,
  #     design_matrix_beta[i,]
  #   )
  #   # Turn Pr(pi[i]) to zero if no species present
  #   pi2[i] <- pi[i] * step(y2[i]-1)
  #   # code up binomial likelihood. JAGS does not have a binomial
  #   # coefficient so it's a bit of a beast. First line
  #   #  is literally just the binomial coefficient.
  #   beta_lik[i] <- (logfact(y2[i]) - (logfact(y1[i]) + logfact(y2[i] - y1[i]))) + 
  #     (y1[i] * log(pi2[i])) + ((y2[i] - y1[i]) * log(1 - pi2[i]))
  #   
  #   # y1 is a binomial process
  #   beta_ones[i] ~ dbern( 
  #     (exp(beta_lik[i])/CONSTANT) 
  #   )
  # }
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
  # Alpha diversity
  for(alphai in 1:npar_alpha){
    beta_alpha[alphai] ~ dnorm(0,0.01)
  }
  # beta diversity
  for(betai in 1:npar_beta){
    beta_log[betai] ~ dnorm(0, 0.01)
    beta_beta[betai] <-  exp(beta_log[betai])
  }
  b0 ~ dlogis(0,1)
}