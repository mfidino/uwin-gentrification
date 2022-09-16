

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
  "./mcmc_output/beta_output/beta_results_collapsed_norm.RDS"
)

data_list <- readRDS(
  "./mcmc_output/beta_output/data_list.RDS"
)


dmat <- read.csv(
  "./mcmc_output/beta_output/dmat_site_ids_collapsed.csv"
)

cities <- sapply(
  strsplit(
    dmat$siteA,
    "-"
  ), 
  "[[",
  1
)
city_map <- data.frame(
  city = unique(cities)
)


city_map$pname <- c(
  "Athens, GA",
  "Bay Area, CA",
  "Boston, MA",
  "Chicago, IL",
  "Denver, CO",
  "Houston, TX",
  "Indianapolis, IN",
  "Iowa City, IA",
  "Jackson, MS",
  "Little Rock, AR",
  "Madison, WI",
  "Metro LA, CA",
  "National Capital",
  "Phoenix, AZ",
  "Portland, OR",
  "Rochester, NY",
  "Sanford, FL",
  "Salt Lake City, UT",
  "Seattle, WA",
  "St. Louis, MO",
  "Tacoma, WA",
  "Urbana, IL",
  "Wilmington, DE"
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
cc <- c("intercept", "geo1", "geo2","geo3", "imp1", "imp2", "imp3","gent")
pdf("./plots/mc_hist_exp.pdf")
for(i in 1:23){
  for(j in 1:8){
   hist(
     mc$beta_exp[,i,j],
     main = paste0(
       city_map$pname[i], ": ",
       cc[j]
     )
   ) 
  }
}
dev.off()


cities <- unique(my_knots$City)

city_pred <- vector("list", length = ncity)


# compare nearby sites vs those gentrifying
mu_1 <-100 * (1 - exp(-mc$beta_mu[,1]))
mu_2 <-100 *  (1 - exp(-rowSums(mc$beta_mu[,c(1,8)])))
mu_1 <- quantile(mu_2 - mu_1, probs = c(0.025,0.5,0.975))


city_mu1 <- matrix(
  ncol = 3,
  nrow = data_list$ncity
)

for(i in 1:data_list$ncity){
  tmp_mu_1 <-100 *  (1 - exp(-mc$beta_exp[,i,1]))
  tmp_mu_2 <-100 * ( 1 - exp(-rowSums(mc$beta_exp[,i,c(1,8)])))
  city_mu1[i,] <- quantile(tmp_mu_2 - tmp_mu_1, probs = c(0.025,0.5,0.975))
}

windows(4,6)
my_order <- order(city_mu1[,2])
par(mar = c(5,1,0.5,0.2), xpd = NA)
bbplot::blank(xlim = c(0,25), ylim = c(0.5,1.5))
bbplot::axis_blank(1)
bbplot::axis_text(side = 1, line = 0.75)
bbplot::axis_text("Percent difference in beta diversity at sites that \nare gentrifying vs. those that are not",
                  side = 1, line = 3.5)


rect(
  xleft = mu_1[1],
  xright = mu_1[3],
  ybottom = 0.5,
  ytop = 1.5,
  col = "gray70",
  border = NA
)
lines(
  x = rep(mu_1[2],2),
  y = c(0.5,1.5),
  lty = 2
)


yloc <- seq(0.5,1.5, length.out = data_list$ncity)
for(i in 1:data_list$ncity){
  tt <- city_mu1[my_order[i],]
  lines(
    x = tt[-2],
    y = rep(yloc[i],2)
  )
  
  text(
    x = city_mu1[my_order[i],3] + 0.1,
    y = yloc[i],
    labels = city_map$pname[my_order[i]],
    pos = 4,
    cex = 0.8
  )
}
points(
  x = city_mu1[my_order,2],
  y = yloc,
  pch = 19
)

# get biggest effect across each covariate per city
overall_effects <- data.frame(
  city = rep(city_map$city,each = 3),
  parm = rep(c("geo","imp","gent"), data_list$ncity),
  lo = NA,
  est = NA,
  hi= NA
)
for(city in 1:ncity){
  geo_tmp <- rowSums(
    mc$beta_exp[,city,2:4]
  )
  imp_tmp <- rowSums(
    mc$beta_exp[,city,5:7]
  )
  geo_tmp <- quantile(geo_tmp, probs = c(0.025,0.5,0.975))
  imp_tmp <- quantile(imp_tmp, probs = c(0.025,0.5,0.975))
  gen_tmp <- quantile(mc$beta_exp[,city,8], probs = c(0.025,0.5,0.975))
  to_send <- rbind(geo_tmp, imp_tmp, gen_tmp)
  overall_effects[
    overall_effects$city == city_map$city[city],
    c("lo","est","hi")
  ] <- to_send
    
}

# get global average as well
overall_mean <- data.frame(
  parm = c("geo","imp","gent")
)
geo_tmp <- rowSums(
    mc$beta_mu[,2:4]
)
imp_tmp <- rowSums(
    mc$beta_mu[,5:7]
)
geo_tmp <- quantile(geo_tmp, probs = c(0.025,0.5,0.975))
imp_tmp <- quantile(imp_tmp, probs = c(0.025,0.5,0.975))
gen_tmp <- quantile(mc$beta_mu[,8], probs = c(0.025,0.5,0.975))
to_send <- rbind(geo_tmp, imp_tmp, gen_tmp)

overall_mean <- cbind(overall_mean,to_send)
colnames(overall_mean) <- colnames(overall_effects)[-1]

overall_effects <- split(
  overall_effects,
  factor(overall_effects$parm)
)
set.seed(888)
lshift <- jitter(rep(0, data_list$ncity), factor = 10 )
{
bbplot::blank(xlim = c(0,0.5), ylim = c(0.5,3.5))
bbplot::axis_blank(1)
bbplot::axis_text(side = 1, line = 0.75)

bbplot::axis_text(
  "Sum of spline coefficients",
  side = 1,
  line = 3
)


bbplot::axis_text(
  c(
    "Gentrification",
    "Geographic\ndistance",
    "Impervious\ncover"
  ),
  at = c(1:3),
  line = -0.75,
  side = 2,
  las = 1
)

for(i in 1:3){
  
  for(city in 1:data_list$ncity){
    lines(
      x = unlist(overall_effects[[i]][city,c("lo","hi")]),
      y = rep(i + lshift[city],2),
      col = scales::alpha("black", 0.5)
    )
  }
  
  points(
    x = overall_effects[[i]]$est,
    y = i + lshift,
    pch = 19,
    col = scales::alpha("black", 0.3)
  )
  
}


overall_mean <- overall_mean[order(overall_mean$parm),]
for(i in 1:3){
  lines(
    x = unlist(overall_mean[i,c("lo","hi")]),
    y = rep(i ,2),
    lwd = 5
  )
}
  points(
    x = overall_mean$est,
    y = c(1:3 ),
    pch = 21,
    bg = "gray70",
    cex = 2
  )
}

legend(
  "bottomright",
  legend = c(
    "Among-city average estimate and 95% CI",
    "City-specific estimate and 95% CI"
  ),
  lty = c(1,1),
  lwd = c(5,1),
  pch = c(21,21),
  pt.cex = c(2,1),
  pt.bg = c("gray70", scales::alpha("black",0.3)),
  col = c("black", scales::alpha("black",0.5)),
  bty = "n",
  pt.lwd = 1
  
)

for(city in 1:ncity){
  knots <- as.numeric(my_knots[
    my_knots$City == cities[city] &
    my_knots$covariate == "mean_19",
    c("min","median","max")])
  tmp_mcmc <- mc$beta_exp[,city,]
  city_pred[[city]] <- spline_pred(
    knots = knots,
    mcmc = tmp_mcmc[,5:7]#,
    #intercept = tmp_mcmc[,1]
  )
  
}

city_ed <- vector("list", length = ncity)

for(city in 1:ncity){
  tmp_mcmc <- mc$beta_exp[,city,]
  tmp_spline <- data_list$spline_matrix[
    data_list$city_vec == city,
  ]
  tmp_dissim <- data_list$dissim[
    data_list$city_vec == city
  ]
  my_data <- tmp_dissim
  my_splines <- tmp_spline
  mcmc_mat <- tmp_mcmc
  
  city_ed[[city]] <- ed_pred(
    mcmc = tmp_mcmc[,5:7]#,
    #intercept = tmp_mcmc[,1]
  )
  
}

plot(1~1, type = "n", ylim = c(-0,0.6), xlim = c(0,0.8),
     xlab = "Impervious cover",
     ylab = "f(Impervious cover)",
     bty = "l",
     las = 1
     )

bbplot::ribbon(
  city_pred[[4]]$x,
  city_pred[[4]]$y[,-2],
  col = "purple",
  alpha = 0.5
)
lines(
  city_pred[[4]]$x,
  city_pred[[4]]$y[,2],
  lwd = 1
)

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
