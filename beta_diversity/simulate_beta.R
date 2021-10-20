#######################################
#
# simulate and fit beta diversity stuff
#
#
#
#######################################

# Step 1. assume that we have perfect detection


set.seed(122)

library(mvtnorm)

nsite <- 50
nspecies <- 30

# generate community averages
bmu <- c(-0.7, 0.5)
bsd <- sqrt(c(1,2))

# Get species specific estimates
bspecies <- rmvnorm(
  nspecies,
  bmu,
  diag(bsd)
)

# get environmental covariate
x <-cbind(1,sort(rnorm(nsite)))

# and get  occupancy
psi <- plogis(x %*% t(bspecies))

# estimate species presence

# get matrix of same dimensions as psi
z <- psi

# fill in the z matrix
z[] <- rbinom(length(z), 1, psi)

ed <- data.frame(
  site = paste0("s", 1:nsite),
  lat = runif(nsite, -47, -44),
  long = runif(nsite, 111, 117),
  x = x[,2]
)

# make spline matrix
mydm <- spline_dm(ed, "site", "long", "lat")

y <- zmat_long(
  z,
  mydm
)
# spin up site matrix for intial modeling
# calculate Y

# matrix math to get the dissimilar species along the upper triangle
# This gets number of dissimilar species along the upper triangle
#aa <- (1 - z) %*% t(z)



# The upper triangle is the union
#bb <-  ncol(z) - ((1 - z) %*%(1- t(z)))

data_list <- list(
  y = y,
  dmat = mydm$dmat,
  n = nrow(y),
  npar = ncol(mydm$dmat)
)

library(runjags)


m2 <- run.jags(
  model = "./JAGS/gdm_nmm_log.R",
  monitor = c("b0", "beta"),
  data = data_list,
  n.chains = 5,
  modules = "glm",
  method = "parallel",
  adapt = 1000,
  burnin = 5000,
  sample = 10000,
  thin = 2
)

my_mcmc <- do.call("rbind", m2$mcmc)

# do spline over the one gradient
emmies <- vector("list", length = nrow(mydm$metadata))
betas <- my_mcmc[,grep("beta\\[", colnames(my_mcmc))]
locies <- matrix(1:ncol(betas), ncol = 3, nrow = nrow(mydm$metadata),
                 byrow = TRUE)
for(i in 1:nrow(locies)){
  tmp <-betas[,locies[i,]]
  emmies[[i]] <- spline_pred(
    as.numeric(mydm$metadata[i,c("min","median","max")]),
    tmp
  )
}

plot(
  emmies[[2]]$y[,2] ~  emmies[[2]]$x, type = 'l',
  xlab = "Environmental covariate",
  ylab = "f(Environmental covariate)",
  bty = "l",
  lwd = 3,
  ylim = c(0,2),
  las = 1,
  main = "the rate of beta diversity change along this gradient"
  )
lines(emmies[[2]]$y[,1] ~  emmies[[2]]$x, lty = 2)
lines(emmies[[2]]$y[,3] ~  emmies[[2]]$x, lty = 2)

yy <- data_list$y[,1] / data_list$y[,2]










