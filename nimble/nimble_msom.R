msom_mcmc <- function(seed, data_list, constant_list){
library(nimble)

msom_code <- nimbleCode({
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
      rho[dp] * z[dp], 
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
      a_species_re[i,j] ~ dnorm(0, 1)
      a_species[i,j] <- a_community[i] + 
        tau_community_omega[i] * a_species_re[i,j]
      
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
      b_species_re[l,j] ~ dnorm(0,1)
      b_species[l,j] <- b_community[l] + tau_community_psi[l] * 
        b_species_re[l,j]
      tau_species_psi[l,j] ~ dgamma(tau_shape_psi[l], tau_rate_psi[l])
      for(k in 1:ncity){
        b_species_city_re[l,j,k] ~ dnorm(0,1)
        b_species_city[l,j,k] <-  b_community[l] +
          tau_community_psi[l] * b_species_re[l,j] +
          tau_species_psi[l,j] * b_species_city_re[l,j,k] 
      }
    }
  }
  theta_community ~ dt(0, 2.5, 1)
  tau_community_theta ~ dgamma(1, 1)
  # hyper priors for species random effects
  tau_shape_theta ~ dunif(0.01, 10) 
  tau_rate_theta ~ dunif(0.01, 10)
  for(j in 1:nspecies){
    theta_species_re[j] ~ dnorm(0,1)
    theta_species[j] <- theta_community +
      tau_community_theta * theta_species_re[j]
    tau_species_theta[j] ~ dgamma(tau_shape_theta, tau_rate_theta)
    for(k in 1:ncity){
      theta_psi_re[j,k] ~ dnorm(0,1)
      theta_psi[j,k] <- theta_community +
        tau_community_theta * theta_species_re[j] +
        tau_species_theta[j] * theta_psi_re[j,k]
    }
  }
  # Detection priors
  for(l in 1:nrho){
    c_community[l] ~ dt(0, 2.5, 1)
    tau_community_rho[l] ~ dgamma(1, 1)
    tau_shape_rho[l] ~ dunif(0.01, 10) # hyper priors for species random effects.
    tau_rate_rho[l] ~ dunif(0.01, 10)
    for(j in 1:nspecies){
      c_species_re[l,j] ~ dnorm(0,1)
      c_species[l,j] <- c_community[l] +
        tau_community_rho[l] * c_species_re[l,j]
      tau_species_rho[l,j] ~ dgamma(tau_shape_rho[l], tau_rate_rho[l])
      for(k in 1:ncity){
        c_species_city_re[l,j,k] ~ dnorm(0,1)
        c_species_city[l,j,k] <- c_community[l] +
          tau_community_rho[l] * c_species_re[l,j] +
          tau_species_rho[l,j] * c_species_city_re[l,j,k]
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
    city_psi_re[m] ~ dnorm(0,1)
    # species, season, and city seasonal variation. (ssc)
    ssc_psi[m] <- city_tau_psi[combo_city_vec[m]] * city_psi_re[m]
    city_rho_re[m] ~ dnorm(0,1)
    ssc_rho[m] <- city_tau_rho[combo_city_vec[m]] * city_rho_re[m]
  }
})


my_inits <- function(){
    list( 
      z = rep(
        1,
        constant_list$ndata
      ),
      x = matrix(
        1,
        ncol = constant_list$ncity,
        nrow = constant_list$nspecies
      ),
      a_community = rnorm(
        constant_list$nomega
      ),
      tau_community_omega = rgamma(
        constant_list$nomega,
        1,
        1
      ),
      a_species_re = matrix(
        rnorm(
          prod(
            unlist(
              constant_list[c('nomega', 'nspecies')]
            )
          )
        ),
        nrow = constant_list$nomega,
        ncol = constant_list$nspecies
      ),
      b_community = rnorm(
        constant_list$npsi
      ),
      tau_community_psi = rgamma(
        constant_list$npsi,
        1,
        1
      ),
      tau_shape_psi = runif(
        constant_list$npsi,
        0.5,
        2
      ),
      tau_rate_psi = runif(
        constant_list$npsi,
        0.5,
        2
      ),
      b_species_re = matrix(
        rnorm(
          prod(
            unlist(
              constant_list[c('npsi', 'npspecies')]
            )
          )
        ),
        nrow = constant_list$npsi,
        ncol = constant_list$nspecies
      ),
      tau_species_psi = matrix(
        rgamma(
          prod(
            unlist(
              constant_list[c('npsi', 'nspecies')]
            )
          ),
          1,
          1
        ),
        nrow = constant_list$npsi,
        ncol = constant_list$nspecies
      ),
      b_species_city_re = array(
        rnorm(
          prod(
            unlist(
              constant_list[c('npsi', 'nspecies', 'ncity')]
            )
          )
        ),
        dim = unlist(
          constant_list[c('npsi', 'nspecies', 'ncity')]
        )
      ),
      theta_community = rnorm(
        1
      ),
      tau_community_theta = rgamma(
        1,
        1,
        1
      ),
      tau_shape_theta = runif(
        1,
        0.5,
        2
      ),
      tau_rate_theta = runif(
        1,
        0.5,
        2
      ),
      theta_species_re = rnorm(
        constant_list$nspecies
      ),
      tau_species_theta = rgamma(
        constant_list$nspecies,
        1,
        1
      ),
      theta_psi_re = matrix(
        rnorm(
          prod(
            unlist(
              constant_list[c('nspecies', 'ncity')]
            )
          )
        ),
        nrow = constant_list$nspecies,
        ncol = constant_list$ncity
      ),
      c_community = rnorm(
        constant_list$nrho
      ),
      tau_community_rho = rgamma(
        constant_list$nrho,
        1,
        1
      ),
      tau_shape_rho = runif(
        constant_list$nrho,
        0.5,
        2
      ),
      tau_rate_rho = runif(
        constant_list$nrho,
        0.5,
        2
      ),
      c_species_re = matrix(
        rnorm(
          prod(
            unlist(
              constant_list[c('nrho', 'nspecies')]
            )
          )
        ),
        nrow = constant_list$nrho,
        ncol = constant_list$nspecies
      ),
      tau_species_rho = matrix(
        rgamma(
          prod(
            unlist(
              constant_list[c('nrho', 'nspecies')]
            )
          ),
          1,
          1
        ),
        nrow = constant_list$nrho,
        ncol = constant_list$nspecies
      ),
      c_species_city_re = array(
        rnorm(
          prod(
            unlist(
              constant_list[c('nrho', 'nspecies', 'ncity')]
            )
          )
        ),
        dim = unlist(
          constant_list[c('nrho', 'nspecies', 'ncity')]
        )
      ),
      city_shape_psi = runif(
        1,
        0.5,
        2
      ),
      city_shape_rho = runif(
        1,
        0.5,
        2
      ),
      city_rate_psi = runif(
        1,
        0.5,
        2
      ),
      city_rate_rho = runif(
        1,
        0.5,
        2
      ),
      city_tau_psi = rgamma(
        constant_list$ncity,
        1,
        1
      ),
      city_tau_rho = rgamma(
        constant_list$ncity,
        1,
        1
      ),
      city_psi_re = rnorm(
        constant_list$nseason_params
      ),
      city_rho_re = rnorm(
        constant_list$nseason_params
      )
      #ssc_psi = rnorm(
      #  constant_list$nseason_params
      #),
      #ssc_rho = rnorm(
      #  constant_list$nseason_params
      #)
    )
}

msom <- nimble::nimbleModel(
  code = msom_code,
  name = "msom",
  constants = constant_list,
  data = data_list,
  dimensions = list(
    z = constant_list$ndata,
    ssc_rho = constant_list$nseason_params,
    ssc_psi = constant_list$nseason_params),
  inits = my_inits()
)

#tmp_msom <- compileNimble(msom)

msom_configure <- nimble::configureMCMC(
  msom,
  monitors = to_monitor,
  monitors2 = "z",
  thin = 1,
  thin2 = 10,
  onlySlice = TRUE
)

msom_configure$removeSampler(
  c("a_community", "a_species_re")
)


msom_build <- nimble::buildMCMC(
  msom_configure
)

msom_compile <- nimble::compileNimble(
  msom, 
  msom_build,
  resetFunctions = TRUE
)


results <- nimble::runMCMC(
  msom_compile$msom_build,
  niter = 50,
  nburnin = 25,
  inits = my_inits()
)

plot(results$samples[,11])

hm <- MCMCvis::MCMCsummary(results$samples[8000:10000,])

round(hm[,3:5],2)
}

cli::cli_alert_info("Model code created as object 'model_code'")

