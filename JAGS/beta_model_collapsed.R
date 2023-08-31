model{
  for(i in 1:ndata){
    # spline dm has the intercept, 3 splines for geographic distance,
    #   3 splines for housing density, and one binary term
    #   for sites being gentrifying
    linpred[i] <-
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

  # city-specific parameters for splines
  for(j in 1:npar_spline){
    beta_mu[j] ~ dnorm(-3, 0.01)
    beta_tau[j] ~ dgamma(1,1)
    beta_sd[j] <- 1 / sqrt(beta_tau[j])
    for(k in 1:ncity){
      beta[k,j] ~ dnorm(beta_mu[j], beta_tau[j])
      beta_exp[k,j] <- exp(beta[k,j])
    }
  }
}