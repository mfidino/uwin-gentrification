model_code <- nimbleCode({
# The multi-species occupancy model
# to add
# among-city regression (latent-state)
  for(j in 1:nspecies){
    for(k in 1:ncity){
      logit(omega_mu[j,k]) <- inprod(
        a_species[1:2,j],
        omegacov[j,k,1:2]
      )
      x[j,k] ~ dbern(
        omega_mu[j,k]
      )
    }
  }
# latent state model
# First seasons of data. fs = first_season
# Note: ssc_psi is the model intercept, b_species_city[1,,] is the
#  average intercept for a species in a city ACROSS seasons of sampling.
#  We are not using it here, but store it to make use of it for model
#  predictions.
  for(fs in 1:nfirsts){
    logit(psi[fs]) <- inprod(
      b_species_city[1:npsi, species_vec[fs], city_vec[fs]],
      psicov[fs, 1:npsi]
    ) + ssc_psi[combo_vec[fs]]
    z[fs] ~ dbern(
      psi[fs] * x[species_vec[fs], city_vec[fs]]
    )
  }
# Remaining seasons of data. tr = "the rest" (of the seasons).
# This is essentially the same as above except it has the
# first-order autologistic term.
  for(tr in (nfirsts + 1):ndata){
    logit(psi[tr]) <- inprod(
      b_species_city[1:npsi, species_vec[tr], city_vec[tr]],
      psicov[tr, 1:npsi]
    ) + 
      ssc_psi[combo_vec[tr]] +
      theta_psi[species_vec[tr], city_vec[tr]] * z[last_sample_vec[tr]]
    z[tr] ~ dbern(
      psi[tr] * x[species_vec[tr], city_vec[tr]]
    )
  }
  # data model. dp = data point
  for(dp in 1:ndata){
    logit(rho[dp]) <- inprod(
      c_species_city[1:nrho, species_vec[dp], city_vec[dp]],
      rhocov[dp, 1:nrho]
    ) + ssc_rho[combo_vec[dp]]
    y[dp] ~ dbin(
      rho[dp], 
      J[dp]
    )
  }
  ############################
  # Priors
  ############################
  #
  # Species presence in sampling area priors
  for(i in 1:nomega){
    a_community[i] ~ dt(0, 2.5, 1)
    tau_community_omega[i] ~ dgamma(1, 1)
  }
  for(i in 1:nomega){
    for(j in 1:nspecies){
      a_species[i,j] ~ dnorm(a_community[i], tau_community_omega[i])
    }
  }
  # Occupancy priors
  # Random intercepts & slopes within city
  for(l in 1:npsi){
    b_community[l] ~ dt(0, 2.5, 1)
    tau_community_psi[l] ~ dgamma(1, 1)
    # hyper priors for species random effects
    tau_shape_psi[l] ~ dunif(0.01, 10) 
    tau_rate_psi[l] ~ dunif(0.01, 10)
    for(j in 1:nspecies){
      b_species[l,j] ~ dnorm(b_community[l], tau_community_psi[l])
      tau_species_psi[l,j] ~ dgamma(tau_shape_psi[l], tau_rate_psi[l])
      for(k in 1:ncity){
        b_species_city[l,j,k] ~ dnorm(b_species[l,j], tau_species_psi[l,j])
      }
    }
  }
  theta_community ~ dt(0, 2.5, 1)
  tau_community_theta ~ dgamma(1, 1)
  # hyper priors for species random effects
  tau_shape_theta ~ dunif(0.01, 10) 
  tau_rate_theta ~ dunif(0.01, 10)
  for(j in 1:nspecies){
    theta_species[j] ~ dnorm(theta_community, tau_community_theta)
    tau_species_theta[j] ~ dgamma(tau_shape_theta, tau_rate_theta)
    for(k in 1:ncity){
      theta_psi[j,k] ~ dnorm(theta_species[j], tau_species_theta[j])
    }
  }
  # Detection priors
  for(l in 1:nrho){
    c_community[l] ~ dt(0, 2.5, 1)
    tau_community_rho[l] ~ dgamma(1, 1)
    tau_shape_rho[l] ~ dunif(0.01, 10) # hyper priors for species random effects.
    tau_rate_rho[l] ~ dunif(0.01, 10)
    for(j in 1:nspecies){
      c_species[l,j] ~ dnorm(c_community[l], tau_community_rho[l])
      tau_species_rho[l,j] ~ dgamma(tau_shape_rho[l], tau_rate_rho[l])
      for(k in 1:ncity){
        c_species_city[l,j,k] ~ dnorm(c_species[l,j], tau_species_rho[l,j])
      }
    }
  }
  # Setting up within-city seasonal variation now that we have
  #  a structure set up for each species in a city.
  city_shape_psi ~ dunif(0.001, 10)
  city_shape_rho ~ dunif(0.001, 10)
  city_rate_psi ~ dunif(0.001, 10)
  city_rate_rho ~ dunif(0.001, 10)
  for(k in 1:ncity){
    city_tau_psi[k] ~ dgamma(city_shape_psi, city_rate_psi)
    city_tau_rho[k] ~ dgamma(city_shape_rho, city_rate_rho)
  }
  # Refererence the correct season_city_tau for the seasonal random effect.
  for(m in 1:nseason_params){
    # species, season, and city seasonal variation. (ssc)
    ssc_psi[m] ~ dnorm(
      0,
      city_tau_psi[combo_city_vec[m]]
    )
    ssc_rho[m] ~ dnorm(
      0,
      city_tau_rho[combo_city_vec[m]]
    )
  }
})

cli::cli_alert_info("Model code created as object 'model_code'")

