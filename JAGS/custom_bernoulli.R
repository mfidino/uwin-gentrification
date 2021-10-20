model{
  for(i in 1:n){
    logit(pi[i]) <- inprod(x[i,], mu)
    pi2[i] <- pi[i] * step(y2[i]-1)
    beta_lik[i] <- (logfact(y2[i]) - (logfact(y1[i]) + logfact(y2[i] - y1[i]))) + 
      (y1[i] * log(pi2[i])) + ((y2[i] - y1[i]) * log(1 - pi2[i]))
    beta_ones[i] ~ dbern( 
      (exp(beta_lik[i])/CONSTANT)
    )
  }
  mu[1] ~ dlogis(0,1)
  mu[2] ~ dlogis(0,1)
}