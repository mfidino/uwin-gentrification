model{
  for(i in 1:ndata){
    # total unique species at each site pair, used for 
    #  variance term in model.
    #total_rich[i] ~ dlnorm(
    #  log_mu[i],
    #  pow(log_sd[i], -2)
    #)
    # linear predictor, all coefficients are
    #  strictly positive.
    linpred[i] <-
      inprod(
        beta_exp,
        design_matrix_beta[i,] 
      )

    # calculate variance for normal distribution using the
    #  binomial variance as per Ferrier 2007.
    # step 1 make probability
    mu[i] <- 1 - exp(-linpred[i])
    # step 2, variance divided by expected species
    #  richness.
    my_var[i] <- (mu[i] * (1 - mu[i])) / total_rich[i]
    # Add variability to the linear predictor based on
    #  the expected richness between two sites and the
    #  estimated similarity among the sites.
    linpred2[i] ~ dnorm(
      linpred[i],
      (1 / my_var[i])
    )T(0,)
    # convert to probability
    mu2[i] <- 1 - exp(-linpred2[i])
    # The data we have, plus how uncertain we are
    #  in this estimate.
    simdata[i] ~ dnorm(
      mu2[i],
      pow(simsd[i],-2)
    )T(0,)
  }
  for(j in 1:npar){
    beta[j] ~ dnorm(0, 0.01)
    beta_exp[j] <- exp(beta[j])
  }
  
}