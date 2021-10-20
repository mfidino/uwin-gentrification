model{
  for(i in 1:n){
    beta_lik[i] <- (logfact(y2[i]) - (logfact(y1[i]) + logfact(y2[i] - y1[i]))) + 
      (y1[i] * log(pi)) + ((y2[i] - y1[i]) * log(1 - pi))
    beta_ones[i] ~ dbern( 
      exp(beta_lik[i])/CONSTANT
    )
  }
  pi ~ dunif(0,1)
}