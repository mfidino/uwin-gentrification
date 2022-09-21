library(sf)
sf::sf_use_s2(FALSE)
library(runjags)
library(dplyr)
library(cli)
library(MCMCvis)

source("./beta_diversity/gdm_functions.R")
source("./R/alpha_beta_functions.R")

# read in alpha diversity analysis


alpha_mc <- readRDS(
  "./mcmc_output/alpha_output/alpha_mcmc.RDS"
)

alpha_mc <- do.call(
  "rbind",
  alpha_mc$mcmc
)

alpha_mc <- split_mcmc(
  alpha_mc
)

# and beta diversity analysis

beta_mc <- readRDS(
  "./mcmc_output/beta_output/beta_results_collapsed_norm.RDS"
)

beta_mc <- do.call(
  "rbind",
  beta_mc$mcmc
)

beta_mc <- split_mcmc(
  beta_mc
)

# get parameter summaries
alpha_parms <- c(
  "int",
  "gent",
  "imp",
  "impxgent"
)

alpha_city <- apply(
  alpha_mc$alpha,
  c(2,3),
  qwrap
)

alpha_mu <- apply(
  alpha_mc$alpha_mu,
  2,
  qwrap
)

alpha_sd <- apply(
  alpha_mc$alpha_sd,
  2,
  qwrap
)

# same with betas
beta_parms <- c(
  "int",
  "geo",
  "imp",
  "gent"
)

overall_effects <- data.frame(
  city = rep(1:data_list$ncity, each = 4),
  parm = rep(c("int", "geo","imp","gent"), data_list$ncity),
  lo = NA,
  est = NA,
  hi= NA
)
for(city in 1:data_list$ncity){
  geo_tmp <- rowSums(
    beta_mc$beta_exp[,city,2:4]
  )
  imp_tmp <- rowSums(
    beta_mc$beta_exp[,city,5:7]
  )
  int_tmp <- qwrap(beta_mc$beta_exp[,city,1])
  geo_tmp <- qwrap(geo_tmp)
  imp_tmp <- qwrap(imp_tmp)
  gen_tmp <- qwrap(beta_mc$beta_exp[,city,8])
  to_send <- rbind(int_tmp, geo_tmp, imp_tmp, gen_tmp)
  overall_effects[
    overall_effects$city == city,
    c("lo","est","hi")
  ] <- to_send
  
}

overall_effects <- split(
  overall_effects,
  factor(overall_effects$city)
)

# get global average as well
overall_mean <- data.frame(
  parm = c("int", "geo","imp","gent")
)
geo_tmp <- rowSums(
  beta_mc$beta_mu[,2:4]
)
imp_tmp <- rowSums(
  beta_mc$beta_mu[,5:7]
)
int_tmp <- qwrap(beta_mc$beta_mu[,1])
geo_tmp <- qwrap(geo_tmp)
imp_tmp <- qwrap(imp_tmp)
gen_tmp <- qwrap(beta_mc$beta_mu[,8])
to_send <- rbind(int_tmp, geo_tmp, imp_tmp, gen_tmp)

overall_mean <- cbind(overall_mean,to_send)

# make a 2 x 1 plot
windows(4,6)
m <- matrix(
    1:2,
    ncol = 1,
    nrow = 2
  )
layout(m)
par(mar = c(2,1,1,1), oma = c(4,6,0,0), xpd = NA)

set.seed(2222)
bbplot::blank( xlim = c(-1.5, 3), ylim = c(0,4.5), xaxs = "i", yaxs = "i")
bbplot::axis_blank(1, tck = -0.03)
bbplot::axis_text(side = 1, line = 0.35)

bbplot::axis_text(
  c(
    #"Intercept    ",
    "Gentrification",
    "Impervious  ",
    "Impervious x\nGentrification"
  ),
  2,
  line = -0,
  at = rev(1:4),
  las = 1
)

for(i in 1:4){
  # location for plotting
  j <- 5 - i

  # my_jitter <- rev(
  #   seq(
  #     -0.45,
  #     0.45,
  #     length.out = data_list$ncity
  #     )
  # )
  my_jitter <- runif(data_list$ncity, -0.4,0.4)
  tmp_vals <- alpha_city[,,i]
  tmp_vals <- tmp_vals[,order(tmp_vals[2,])]
  for(k in 1:data_list$ncity){
    # lines(
    #   x = tmp_vals[c(1,3), k],
    #   y = rep(j+my_jitter[k], 2),
    #   col = scales::alpha("black", 0.3)
    # )
  }
  points(
    x = tmp_vals[2,],
    y = j + my_jitter,
    bg = scales::alpha("#00AADE", 0.7),
    pch = 21
  )
  
  rect(
    xleft = alpha_mu[1,i],
    xright = alpha_mu[3,i],
    ybottom = j - 0.45,
    ytop = j + 0.45,
    col = scales::alpha("gray40",0.5),
    border = NA
  )
  lines(
    x = rep(alpha_mu[2,i],2),
    y = c(j-0.45, j + 0.45),
    lwd = 4,
    lty = 1,
    lend = 1,
    col = "gray30"
  )


}

bbplot::axis_text(
  "Alpha diversity",
  side = 2,
  line = 5.5,
  cex = 1.2
)




windows(4,6)
m <- matrix(
  1:2,
  ncol = 1,
  nrow = 2
)
layout(m)
par(mar = c(2,1,1,1), oma = c(2,6,0,0), xpd = NA)
{
set.seed(2222)
bbplot::blank( xlim = c(-1, 1), ylim = c(0,3.5), xaxs = "i", yaxs = "i")
bbplot::axis_blank(1, tck = -0.03,
                   at = seq(-1, 1, 0.5))
bbplot::axis_text(side = 1, line = 0.35,
                  text = seq(-1,1, 0.5),
                  at = seq(-1, 1, 0.5))
lines(
  x = c(0,0),
  y = c(0,3.5),
  lty = 2
)
bbplot::axis_text(
  c(
    #"Intercept    ",
    "Gentrification",
    "Impervious  ",
    "Impervious x\nGentrification"
  ),
  2,
  line = -0,
  at = rev(1:3),
  las = 1
)

for(i in 1:3){
  # location for plotting
  j <- 4 - i
  
  # my_jitter <- rev(
  #   seq(
  #     -0.45,
  #     0.45,
  #     length.out = data_list$ncity
  #     )
  # )
  my_jitter <- runif(data_list$ncity, -0.4,0.4)
  tmp_vals <- alpha_city[,,i+1]
  tmp_vals <- tmp_vals[,order(tmp_vals[2,])]
  for(k in 1:data_list$ncity){
    # lines(
    #   x = tmp_vals[c(1,3), k],
    #   y = rep(j+my_jitter[k], 2),
    #   col = scales::alpha("black", 0.3)
    # )
  }
  points(
    x = tmp_vals[2,],
    y = j + my_jitter,
    bg = scales::alpha("#00AADE", 0.7),
    pch = 21
  )
  
  rect(
    xleft = alpha_mu[1,i+1],
    xright = alpha_mu[3,i+1],
    ybottom = j - 0.45,
    ytop = j + 0.45,
    col = scales::alpha("gray40",0.5),
    border = NA
  )
  lines(
    x = rep(alpha_mu[2,i+1],2),
    y = c(j-0.45, j + 0.45),
    lwd = 4,
    lty = 1,
    lend = 1,
    col = "gray10"
  )
  
  
}

bbplot::axis_text(
  "Alpha diversity",
  side = 2,
  line = 5.5,
  cex = 1.2
)

# and the same for beta diversity

set.seed(2222)
bbplot::blank( xlim = c(0, 0.4), ylim = c(0,3.5), xaxs = "i", yaxs = "i")
bbplot::axis_blank(1, tck = -0.03)
bbplot::axis_text(side = 1, line = 0.35)

bbplot::axis_text(
  c(
    #"Intercept    ",
    "Geographic\ndistance",
    "Impervious",
    "Gentrification"
  ),
  2,
  line = -0,
  at = rev(1:3),
  las = 1
)

for(i in 1:3){
  # location for plotting
  j <- 4 - i
  
  # my_jitter <- rev(
  #   seq(
  #     -0.45,
  #     0.45,
  #     length.out = data_list$ncity
  #     )
  # )
  my_jitter <- runif(data_list$ncity, -0.4,0.4)
  # for(k in 1:data_list$ncity){
  #   # lines(
  #   #   x = tmp_vals[c(1,3), k],
  #   #   y = rep(j+my_jitter[k], 2),
  #   #   col = scales::alpha("black", 0.3)
  #   # )
  # }
  tmp_vals <- sapply(
    overall_effects,
    function(x){
      x[i+1,4]
    }
  )
  points(
    x = tmp_vals,
    y = j + my_jitter,
    bg = scales::alpha("#00AADE", 0.7),
    pch = 21
  )
  
  rect(
    xleft = overall_mean[i+1,2],
    xright = overall_mean[i+1,4],
    ybottom = j - 0.45,
    ytop = j + 0.45,
    col = scales::alpha("gray40",0.5),
    border = NA
  )
  lines(
    x = rep(overall_mean[i+1,3],2),
    y = c(j-0.45, j + 0.45),
    lwd = 4,
    lty = 1,
    lend = 1,
    col = "gray10"
  )
  
  
}

bbplot::axis_text(
  "Beta diversity",
  side = 2,
  line = 5.5,
  cex = 1.2
)
bbplot::axis_text(
  "Effect size",
  side = 1,
  line = 2,
  cex = 1.2
)
rect()

}

p <- seq(0.001,0.999, length.out = 20)

cl <- function(x) -log(1-x)

pp <- cl(p)

icl <- function(x) 1 - exp(-x)

pp2 <- icl(pp)
