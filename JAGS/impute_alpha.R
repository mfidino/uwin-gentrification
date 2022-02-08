model{
  #----------------------
  # Alpha diversity model
  #----------------------
  for(ai in 1:ndata_alpha){
    log(mu_alpha[ai]) <- inprod(
      alpha[city_vec_alpha[ai],1:npar_alpha],
      design_matrix_alpha[ai,1:npar_alpha]
    ) + alpha_site_re[site_vec_alpha[ai]]
    # Poisson likelihood
    # rich is a Poisson random variable
    alphaz[ai] ~ dpois(mu_alpha[ai])
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
      # alpha_re[city,alphai] ~ dnorm(0,1)
      # alpha[city,alphai] <- 
      #   alpha_mu[alphai] + 
      #   alpha_tau[alphai] * alpha_re[city,alphai]
    }
  }
  #alpha_site_shape ~ dunif(0.001, 100)
  #alpha_site_rate ~ dunif(0.001, 100)
   #for(city in 1:ncity){
   #  alpha_site_tau[city] ~ dgamma(alpha_site_shape,alpha_site_rate)
   #  alpha_site_sd[city] <- 1/sqrt(alpha_site_tau[city])
   #}
  alpha_site_tau ~ dgamma(1,1)
  alpha_site_sd <- 1/sqrt(alpha_site_tau)
   for(site in 1:nsite){
     alpha_site_re[site] ~ dnorm(0,alpha_site_tau)
   }
}
