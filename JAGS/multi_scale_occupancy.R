model{
  # priors for among-city variability, from past research
  # we have strong evidence this slope term is positive.
  a_among[1] ~ dt(0, 2.5, 1)
  a_among[2] ~ dt(2, 2.5, 1)
  for(i in 1:ncov_among){
    #a_among[i] ~ dt(0, 2.5, 1)
    tau_among[i] ~ dgamma(1, 1)
  }
  for(i in 1:ncov_among){
    for(j in 1:nspecies){
      b_among[i,j] ~ dnorm(a_among[i], tau_among[i])
    }
  }
  # among-city regression (latent-state)
  for(j in 1:nspecies){
    for(k in 1:ncity){
      logit(among_mu[j,k]) <- inprod(b_among[,j], among_covs[j,k,])
      x[j,k] ~ dbern(among_mu[j,k])
    }
  }
  # prior for random intercepts & slopes within city
  for(l in 1:ncov_within){
    b_within[l] ~ dt(0, 2.5, 1)
    tau_within[l] ~ dgamma(1, 1)
    tau_shape[l] ~ dunif(0.01, 10) # hyper priors for species random effects.
    tau_rate[l] ~ dunif(0.01, 10)
    for(j in 1:nspecies){
      b_species[l,j] ~ dnorm(b_within[l], tau_within[l])
      tau_species[l,j] ~ dgamma(tau_shape[l], tau_rate[l])
      for(k in 1:ncity){
        b_species_city[l,j,k] ~ dnorm(b_species[l,j], tau_species[l,j])
      }
    }
  }
  # prior for auotlogistic term
  theta_mu ~ dt(0, 2.5, 1)
  tau_theta ~ dgamma(1, 1)
  theta_shape ~ dunif(0.1, 10)
  theta_rate ~ dunif(0.1, 10)
  for(j in 1:nspecies){
    theta_species[j] ~ dnorm(theta_mu, tau_theta)
    tau_theta_species[j] ~ dgamma(theta_shape, theta_rate)
    for(k in 1:ncity){
      theta[j,k] ~ dnorm(theta_species[j], tau_theta_species[j])
    }
  }
  # Setting up within-city seasonal variation now that we have
  #  a structure set up for each species in a city.
  c_shape_psi ~ dunif(0.001, 10)
  c_rate_psi ~ dunif(0.001, 10)
  for(k in 1:ncity){
    c_tau_psi[k] ~ dgamma(c_shape_psi, c_rate_psi)
  }
  # Refererence the correct season_city_tau for the seasonal random effect.
  for(m in 1:nseason_params){
    # species, season, and city seasonal variation.
    ssc_psi[m] ~ dnorm(
      b_species_city[
        1,
        combo_species_idx[m],
        combo_city_idx[m]
        ],
      c_tau_psi[combo_city_idx[m]]
    )
  }
  # within-city regression (latent-state) for 'first' surveys
  for(s in 1:nsamples_one){
    logit(psi[s]) <- inprod(
      b_species_city[2:ncov_within,species_idx[s],city_idx[s]],
      psi_covs[s,2:ncov_within]
    ) + 
      ssc_psi[combo_idx[s]]
    z[s] ~ dbern(psi[s] * x[species_idx[s], city_idx[s]])
  }
  # within-city regression (latent-state) with auto-logistic term
  for(ss in (nsamples_one+1):nsamples_two){
    logit(psi[ss]) <- inprod(
      b_species_city[2:ncov_within,species_idx[ss],city_idx[ss]],
      psi_covs[ss,2:ncov_within]
    ) + 
      ssc_psi[combo_idx[ss]] + 
      theta[species_idx[ss], city_idx[ss]] * z[last_sample_vec[ss]]
    z[ss] ~ dbern(psi[ss] * x[species_idx[ss], city_idx[ss]])
    
  }
  # prior for detection
  # prior for random intercepts & slopes within city
  for(l in 1:ncov_det){
    c_det[l] ~ dt(0, 2.5, 1)
    tau_det[l] ~ dgamma(1, 1)
    tau_shape_rho[l] ~ dunif(0.01, 10) # hyper priors for species random effects.
    tau_rate_rho[l] ~ dunif(0.01, 10)
    for(j in 1:nspecies){
      c_species_det[l,j] ~ dnorm(c_det[l], tau_det[l])
      tau_species_det[l,j] ~ dgamma(tau_shape_rho[l], tau_rate_rho[l])
      for(k in 1:ncity){
        c_species_city[l,j,k] ~ dnorm(c_species_det[l,j], tau_species_det[l,j])
      }
    }
  }
  # Setting up within-city seasonal variation now that we have
  #  a structure set up for each species in a city.
  c_shape_rho ~ dunif(0.001, 10)
  c_rate_rho ~ dunif(0.001, 10)
  for(k in 1:ncity){
    c_tau_rho[k] ~ dgamma(c_shape_rho, c_rate_rho)
  }
  # Refererence the correct season_city_tau for the seasonal random effect.
  for(m in 1:nseason_params){
    # species, season, and city seasonal variation.
    ssc_rho[m] ~ dnorm(
      c_species_city[
        1,
        combo_species_idx[m],
        combo_city_idx[m]
        ],
      c_tau_rho[combo_city_idx[m]]
    )
  }
  # within-city regression (observation)
  for(s in 1:nsamples_two){
    logit(rho[s]) <- inprod(
      c_species_city[2:ncov_det,species_idx[s],city_idx[s]],
      rho_covs[s,2:ncov_det]
    ) + ssc_rho[combo_idx[s]]
    y[s] ~ dbinom(rho[s] * z[s], J[s])
  }
}