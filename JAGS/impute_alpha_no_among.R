model{
  #----------------------
  # Alpha diversity model
  #----------------------

  
  for(ai in 1:ndata_alpha){
    log(mu_alpha[ai]) <- inprod(
      alpha[city_vec_alpha[ai],1:npar_alpha],
      design_matrix_alpha[ai,1:npar_alpha]
    ) + alpha_site_re[site_vec_alpha[ai]] + # site random effect
      resid[ai] # residual error
    # Poisson likelihood
    # rich is a Poisson random variable
    alpha_z[ai] ~ dnorm(
      mu_alpha[ai],
      pow(alpha_sd_known[ai], -2)
    )
  }
  #-------
  # priors
  #-------
  # Alpha diversity
  for(alphai in 1:npar_alpha){
    alpha_mu[alphai] ~ dlogis(0,1)
    alpha_tau[alphai] ~ dgamma(1,1)
    alpha_sd[alphai] <- 1/sqrt(alpha_tau[alphai])
    for(city in 1:ncity){
      alpha[city,alphai] ~ dnorm(
        alpha_mu[alphai],
        alpha_tau[alphai]
      )
    }
  }
  alpha_site_tau ~ dgamma(1,1)
  alpha_site_sd <- 1/sqrt(alpha_site_tau)
   for(site in 1:nsite){
     alpha_site_re[site] ~ dnorm(0,alpha_site_tau)
   }
  resid_tau ~ dgamma(1,1)
  resid_sd <- 1 / sqrt(resid_tau)
  for(ai in 1:ndata_alpha){
    resid[ai] ~ dnorm(0, resid_tau)
  }
}
