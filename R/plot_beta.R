

library(sf)
sf::sf_use_s2(FALSE)

library(runjags)
library(dplyr)
library(cli)
library(MCMCvis)

source("./beta_diversity/gdm_functions.R")
source("./R/alpha_beta_functions.R")

# read in the knots
my_knots <- read.csv(
  "./mcmc_output/beta_output/knots.csv"
)

# read in the mcmc
my_mcmc <- readRDS(
  "./mcmc_output/beta_output/beta_results_collapsed.RDS"
)

my_mcmc <- do.call(
  "rbind",
  my_mcmc$mcmc
)

# generate predictions along a gradient for
#  each site
ncity <- length(
  unique(
    my_knots$City
  )
)
mc <- split_mcmc(my_mcmc)

cities <- unique(my_knots$City)

city_pred <- vector("list", length = ncity)
for(city in 1:ncity){
  knots <- as.numeric(my_knots[
    my_knots$City == cities[city] &
    my_knots$covariate == "mean_19",
    c("min","median","max")])
  tmp_mcmc <- mc$beta_exp[,city,]
  city_pred[[city]] <- spline_pred(
    knots = knots,
    mcmc = tmp_mcmc[,5:7],
    intercept = tmp_mcmc[,1]
  )
  
}



plot(1~1, type = "n", ylim = c(-0,0.6), xlim = c(0,0.8),
     xlab = "Impervious cover",
     ylab = "f(Impervious cover)",
     bty = "l",
     las = 1
     )
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "gray70", border = "black")
tmpcol <- rainbow(ncity)
for(city in 1:ncity){
  lines(
    city_pred[[city]]$y[,2] ~ city_pred[[city]]$x,
    col = scales::alpha(tmpcol[city], 0.8),
    lwd = 5
  )

}

legend(
  "topleft",
  unique(sp_dat$City)[c(4,3,5,1,2,6)],
  col = tmpcol[c(4,3,5,1,2,6)], lwd = 5,
  bty = "n"
  )
