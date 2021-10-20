model{
  for(i in 1:n){
    logit(pi[i]) <- b0 + inprod(
      beta,
      dmat[i,]
    )
    y[i,1] ~ dbin(
      pi[i],
      y[i,2]
    )
  }
  for(j in 1:npar){
    beta_log[j] ~ dnorm(0, 0.01)
    beta[j] <-  exp(beta_log[j])
  }
  b0 ~ dlogis(0,1)
}