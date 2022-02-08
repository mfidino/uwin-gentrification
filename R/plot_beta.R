

library(sf)
library(runjags)
library(dplyr)
library(cli)
library(MCMCvis)

# read in the knots
my_knots <- read.csv(
  "./mcmc_output/beta_output/knots.csv"
)

# read in the mcmc
my_mcmc <- read.csv(
  "./mcmc_output/beta_output/beta_mcmc.csv"
)

# generate predictions along a gradient for
#  each site
ncity <- length(
  unique(
    my_knots$City_id
  )
)
city_pred <- vector("list", length = ncity)
for(city in 1:ncity){
  knots <- as.numeric(my_knots[city,c("min","median","max")])
  tmp_mcmc <- grab(
    my_mcmc,
    paste0(
      "beta_exp\\.",city,"\\."
    )
  )
  city_pred[[city]] <- spline_pred(
    knots = knots,
    mcmc = tmp_mcmc,
    intercept = grab(my_mcmc,paste0("b0\\.",city,"\\."))
  )
  
}


plot(1~1, type = "n", ylim = c(-0,1), xlim = c(0,0.8),
     xlab = "geographic distance",
     ylab = "f(geographic distance)",
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
