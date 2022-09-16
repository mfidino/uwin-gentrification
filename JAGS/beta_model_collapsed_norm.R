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
    mu_with_binomial_var[i] ~ dnorm(
      mu[i],
      (1 / binomial_var[i])
    )T(0,)
    # convert to probability
    #mu_with_binomial_var[i] <- 1 - exp(-linpred_with_binomial_var[i])
    # We have uncertainty in the beta diversity estimate, so we
    #  need to add in that uncertainty here.
    dissim[i] ~ dnorm(
      mu_with_binomial_var[i],
      (1 / var_dissim[i])
    )T(0,)
  }

  # city-specific parameters for splines
  for(j in 1:npar_spline){
    beta_mu[j] ~ dnorm(0, 0.1)T(0,)
    tau[j] ~ dgamma(1, 1)
    # use moment of methods to get mean and sd
    #  for each spline.
    #beta_mu[j] <- beta_a[j] / beta_b[j]
    beta_sd[j] <- 1 / sqrt(tau[j])
    for(k in 1:ncity){
      # need to add a small value to the scale parameter or else
      #  we get stuck at infinite density.
      beta_exp[k,j] ~ dnorm(beta_mu[j], tau[j])T(0,)
     # beta_exp[k,j] ~ dgamma(beta_a[j] + 1.0E-3, beta_b[j])
    }
  }
}