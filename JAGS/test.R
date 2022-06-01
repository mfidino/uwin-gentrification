model{
  y ~ dlnorm(log_mu, pow(log_sd, -2))
  y2 ~ dnorm(mu, pow(sd, -2))
}