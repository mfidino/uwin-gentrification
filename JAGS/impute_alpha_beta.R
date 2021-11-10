model{
  #------------------------------
  # multi-species occupancy model
  #------------------------------
  #----------------------
  # Alpha diversity model
  #----------------------
  for(site in 1:nsite){
    log(mu_alpha[site]) <- inprod(
      beta_alpha,
      design_matrix_alpha[site,]
      )
    # Poisson likelihood
    # rich is a Poisson random variable
    alphaz[site] ~ dpois(mu_alpha[site])
  }
  #---------------------
  # Beta diversity model
  #---------------------
  for(i in 1:n){
    # Linear predictor
    logit(pi[i]) <- b0 + inprod(
      beta_beta,
      design_matrix_beta[i,]
    )
    # Turn Pr(pi[i]) to zero if no species present
    pi2[i] <- pi[i] * step(betaz[i,2]-1)
    # code up binomial likelihood. JAGS does not have a binomial
    # coefficient so it's a bit of a beast. First line
    #  is literally just the binomial coefficient.
    # y1 is a binomial process
    betaz[i,1] ~ dbin(
      pi2[i],
      betaz[i,2]
    )
  }
  #-------
  # priors
  #-------
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