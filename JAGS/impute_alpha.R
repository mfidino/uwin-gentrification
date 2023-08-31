model{
  #----------------------
  # Alpha diversity model
  #----------------------
  # fit the model, first time sampling
  for(ai in 1:ndata_alpha){
    log(mu_alpha[ai]) <- inprod(
      alpha[city_vec_alpha[ai],1:npar_alpha],
      design_matrix_alpha[ai,1:npar_alpha]
    ) + resid[ai] # residual error
    alpha_z[ai] ~ dnorm(
      mu_alpha[ai],
      pow(alpha_sd_known[ai], -2)
    )
  }
  #-------
  # priors
  #-------
  # intercept and slope terms
  # for(alphai in 1:npar_alpha){
  #   for(betai in 1:npar_among){
  #     alpha_mu[alphai,betai] ~ dt(0,2.5,1)
  #   }
  #   alpha_tau[alphai] ~ dgamma(1,1)
  #   alpha_sd[alphai] <- 1/sqrt(alpha_tau[alphai])
  #   for(city in 1:ncity){
  #     alpha[city,alphai] ~ dnorm(
  #       inprod(alpha_mu[alphai,], design_matrix_among[city,]),
  #       alpha_tau[alphai]
  #     )
  #   }
  # }
  
  # intercept and slope terms
  for(alphai in 1:npar_alpha){
   alpha_mu[alphai] ~ dt(0, 2.5, 1)
   alpha_tau[alphai] ~ dgamma(1,1)
   alpha_sd[alphai] <- 1/sqrt(alpha_tau[alphai])
   for(city in 1:ncity){
     alpha[city,alphai] ~ dnorm(
       alpha_mu[alphai],
       alpha_tau[alphai]
     )
   }
  }
 # added residual error
  resid_tau ~ dgamma(1,1)
  resid_sd <- 1 / sqrt(resid_tau)
  for(ai in 1:ndata_alpha){
    resid[ai] ~ dnorm(0, resid_tau)
  }
  
}