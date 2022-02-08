model{
  #---------------------
  # Beta diversity model
  #---------------------
  for(bi in 1:ndata_beta){
    # Linear predictor
    logit(pi[bi]) <- b0[city_vec_beta[bi]] + 
    inprod(
      beta_exp[city_vec_beta[bi],],
      design_matrix_beta[bi,] 
    ) + 
    beta_site_re[siteA_id[bi]] +
    beta_site_re[siteB_id[bi]]
    # Turn Pr(pi[i]) to zero if no species present
    pi2[bi] <- pi[bi] * step(betaz[bi,2]-1)
    betaz[bi,1] ~ dbin(
      pi2[bi],
      betaz[bi,2]
    )
  }
  #-------
  # priors
  #-------
  # beta diversity slope terms
  for(betai in 1:npar_beta){
    beta_mu[betai] ~ dnorm(0, 0.01)
    beta_tau[betai] ~ dgamma(1, 1)
    for(city in 1:ncity){
      beta_log[city,betai] ~ dnorm(
        beta_mu[betai],
        beta_tau[betai]
      )
      beta_exp[city,betai] <- exp(
        beta_log[city, betai]
      )
    }
  }
  # random intercept for beta diversity model
  b0_mu ~ dlogis(0,1)
  b0_tau ~ dgamma(1,1)
  b0_sd <- 1/sqrt(b0_tau)
  for(city in 1:ncity){
    b0[city] ~ dnorm(b0_mu,b0_tau)
  }
  beta_site_tau ~ dgamma(1,1)
  beta_site_sd <- 1/sqrt(beta_site_tau)
  for(site in 1:nsite){
    beta_site_re[site] ~ dnorm(0, beta_site_tau)
  }
}