model{
  for(i in 1:ndata){
    # total unique species at each site pair, used for 
    #  variance term in model.
    # total_rich[i] ~ dlnorm(
    #  log_mu[i],
    #  log_tau[i]
    # )
    # The 'intercept' of the model, basically the
    #  gentrification components.
    intercept[i] <- exp(
      inprod(
        gamma[city_vec[i], ],
        gent_dm[i,]
      ) + 
        site_re[siteA_vec[i]] +
        site_re[siteB_vec[i]]
    )
    # adding in the spline design matrix to account for
    #  geographic distance and impervious cover.
    linpred[i] <- intercept[i] + 
      inprod(
        beta_exp[city_vec[i], ],
        spline_matrix[i, ] 
      )
    # calculate variance for normal distribution using the
    #  binomial variance as per Ferrier 2007.
    # step 1, convert linpred to probability
    mu[i] <- 1 - exp(-linpred[i])
    # step 2, Calculate the binomial variance
    binomial_var[i] <- (mu[i] * (1 - mu[i])) / total_rich[i]
    # Add variability to the linear predictor based on
    #  the expected richness between two sites and the
    #  estimated similarity among the sites.
    linpred_with_binomial_var[i] ~ dnorm(
      linpred[i],
      (1 / binomial_var[i])
    )T(0,)
    # convert to probability
    mu_with_binomial_var[i] <- 1 - exp(-linpred_with_binomial_var[i])
    # We have uncertainty in the beta diversity estimate, so we
    #  need to add in that uncertainty here.
    dissim[i] ~ dnorm(
      mu_with_binomial_var[i],
      (1 / var_dissim[i])
    )T(0,)
  }
  # site random effects
  site_mu ~ dnorm(0, 0.01)
  site_tau ~ dgamma(1,1)
  site_sd <- 1 / sqrt(site_tau)
  for(s in 1:nsite){
    site_re[s] ~ dnorm(site_mu, site_tau)
  }
  # city-specific parameters for gentrification terms
  for(l in 1:npar_gent){
    gamma_mu[l] ~ dnorm(0, 0.01)
    gamma_tau[l] ~ dgamma(1,1)
    gamma_sd[l] <- 1 / sqrt(gamma_tau[l])
    for(k in 1:ncity){
      gamma[k,l] ~ dnorm(gamma_mu[l], gamma_tau[l])
    }
  }
  # city-specific parameters for splines
  for(j in 1:npar_spline){
    beta_mu[j] ~ dnorm(0, 0.01)
    beta_tau[j] ~ dgamma(1,1)
    beta_sd[j] <- 1 / sqrt(beta_tau[j])
    for(k in 1:ncity){
      beta[k,j] ~ dnorm(beta_mu[j], beta_tau[j])
      beta_exp[k,j] <- exp(beta[k,j])
    }
  }
}