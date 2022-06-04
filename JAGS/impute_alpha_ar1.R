model{
  #----------------------
  # Alpha diversity model
  #----------------------

  # simulate richness with error,
  #  using log scale to keep positive, as
  #  we are logging the estimated richness
  #  (to scale it a bit).
  for(i in 1:ndata_alpha){
    tmp_rich[i] ~ dlnorm(
      log_mu[i],
      pow(log_sd[i], -2)
      )
    # log the richness to bring the range down
    log_rich[i] <- log(tmp_rich[i])
  }
  # fit the model, first time sampling
  for(ai in 1:alpha_first){
    log(mu_alpha[ai]) <- inprod(
      alpha[city_vec_alpha[ai],1:npar_alpha],
      design_matrix_alpha[ai,1:npar_alpha]
    ) + resid[ai] # residual error
    # Poisson likelihood
    # rich is a Poisson random variable
    alpha_z[ai] ~ dnorm(
      mu_alpha[ai],
      pow(alpha_sd_known[ai], -2)
    )
  }
  # if there was an AR(1) sample collected, add
  #  the auto-regressive term
  for(aii in (alpha_first+1):alpha_resample){
    log(mu_alpha[aii]) <- inprod(
      alpha[city_vec_alpha[aii],1:npar_alpha],
      design_matrix_alpha[aii,1:npar_alpha]
    ) + resid[aii] +# residual error 
    theta[city_vec_alpha[aii]] *
      (log_rich[last_sample_vec[aii]])
    # Poisson likelihood
    # rich is a Poisson random variable
    alpha_z[aii] ~ dnorm(
      mu_alpha[aii],
      pow(alpha_sd_known[aii], -2)
    )
  }
  #-------
  # priors
  #-------
  # Alpha diversity AR(1) term
  theta_mu ~ dt(0,2.5,1)
  theta_tau ~ dgamma(1,1)
  theta_sd <- 1/sqrt(theta_tau)
  for(city in 1:ncity){
    theta[city] ~ dnorm(
      theta_mu,
      theta_tau
    )
  }
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
var mu_alpha[ndata_alpha];