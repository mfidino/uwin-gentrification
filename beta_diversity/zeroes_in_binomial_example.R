# An example of how the zeroes in the beta diversity analysis
#  do not have an influence on the model parameters.

library(runjags)

# Step 1. Simulate some binomial data without zeroes
set.seed(500)
# Number of "species" at each site pair.
y2 <- c(rep(10,400))

# environmental covariate
x <- cbind(1, rnorm(length(y2)))
# Probability of success
pi <- plogis(x %*%  c(0.5, -1))

# Number of dissimilar "species" among site paris.
y1 <- rbinom(length(y2), y2, pi)

# Data for the model
data_list <- list(
  y1 = y1,
  y2 = y2,
  x = x,
  CONSTANT = 10000,
  n = length(y2),
  beta_ones = rep(1, length(y2))
)

# Fit the model
m1 <- run.jags(
  model = "./JAGS/custom_bernoulli.R",
  monitor = c("mu"),
  data = data_list
)

# Summarise the model
msum <- summary(m1)

# Simulate some more data, but with zeroes this time.
set.seed(500)
y2 <- c(rep(10,400), rep(0,400))

# Add on some extra "sites"
x <- cbind(1, c(x[,2], rnorm(length(y2) - nrow(x))))
# Probability of success
pi <- plogis(x %*%  c(0.5, -1))
# Add zeroes for dissimilar species as well.
y1 <- c(y1, rep(0,length(y2) - length(y1)))

# data for model
data_list <- list(
  y1 = y1,
  y2 = y2,
  x = x,
  CONSTANT = 10000,
  n = length(y2),
  beta_ones = rep(1, length(y2))
)

# fit the model
m2 <- run.jags(
  model = "./JAGS/custom_bernoulli.R",
  monitor = c("mu"),
  data = data_list
)

# Summaries
msum2 <- summary(m2)
round(msum2,2)

# bind them together
tmp <- rbind(msum[,1:3], msum2[,1:3])


row.names(tmp) <- c(
  "intercept-model w/o zeroes", "slope-model w/o zeroes",
  "intercept-model w/ zeroes", "slope-model w/ zeroes"
)

tmp
# Outputs are essentially identical
#                              Lower95     Median    Upper95
#intercept-model w/o zeroes  0.4199253  0.4877677  0.5595989
#slope-model w/o zeroes     -1.0954611 -1.0136315 -0.9322009
#intercept-model w/ zeroes   0.4180432  0.4875917  0.5599157
#slope-model w/ zeroes      -1.0996608 -1.0132789 -0.9308773





