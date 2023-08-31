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
  "./mcmc_output/beta_output/beta_results_collapsed_norm_vegan.RDS"
)

data_list <- readRDS(
  "./mcmc_output/beta_output/data_list.RDS"
)


dmat <- read.csv(
  "./mcmc_output/beta_output/dmat_site_ids_collapsed.csv"
)

covs <- read.csv(
  "./data/cleaned_data/covariates/site_covariates.csv"
)

# match site names to dmat
covs$Site <- paste0(
  covs$City,"-",covs$Site
)

# add covariates to dmat to calculate impervious cover
#  distance
dcovs <- dplyr::left_join(
  dmat,
  covs,
  by = c("siteA" = "Site")
)

dcovs <- dplyr::left_join(
  dcovs,
  covs,
  by = c("siteB" = "Site")
)

dcovs$imp_diff <- abs(
  dcovs$mean_19.x - dcovs$mean_19.y
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

# pretty names for each city
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
  "Washington D.C.",
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

# order of parameters
cc <- c("intercept", "geo1", "geo2","geo3", "imp1", "imp2", "imp3","gent")


cities <- unique(my_knots$City)

city_pred <- vector("list", length = ncity)


# compare nearby sites vs those gentrifying
mu_1 <-(1 - exp(-mc$beta_mu[,1]))
mu_2 <-(1 - exp(-rowSums(mc$beta_mu[,c(1,8)])))
quantile(mu_1)

mu_diff <- quantile(mu_2 - mu_1, probs = c(0.025, 0.05, 0.5,0.95, 0.975))
round(mu_diff,2)

# get the same bit of info for some west coast cities
la_1 <-(1 - exp(-mc$beta_exp[,12,1]))
la_2 <-(1 - exp(-rowSums(mc$beta_exp[,12,c(1,8)])))
quantile(mu_1)

la_diff <- quantile(la_2 - la_1, probs = c(0.025, 0.05, 0.5,0.95, 0.975))
round(la_diff,2)
round(
  quantile(
    la_2,
    probs = c(0.025, 0.05, 0.5,0.95, 0.975)
  ),2
)
round(
  quantile(
    la_1,
    probs = c(0.025, 0.05, 0.5,0.95, 0.975)
  ),2
)
round(
  quantile(
    la_2 / la_1,
    probs = c(0.025, 0.05, 0.5,0.95, 0.975)
  ),2
)
round(
  quantile(
    (la_2 / la_1) / (mu_2 / mu_1),
    probs = c(0.025, 0.05, 0.5,0.95, 0.975)
  ),2
)

city_mu1 <- matrix(
  ncol = 3,
  nrow = data_list$ncity
)

for(i in 1:data_list$ncity){
  tmp_mu_1 <-(1 - exp(-mc$beta_exp[,i,1]))
  tmp_mu_2 <-( 1 - exp(-rowSums(mc$beta_exp[,i,c(1,8)])))
  city_mu1[i,] <- quantile(tmp_mu_2 - tmp_mu_1, probs = c(0.025,0.5,0.975))
}


tiff(
  "./plots/beta_avg_gentrification.tiff",
  width = 4,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)
my_order <- order(city_mu1[,2])
par(mar = c(5,1,0.5,0.2), xpd = NA)
bbplot::blank(xlim = c(0,.25), ylim = c(0.5,1.5))
bbplot::axis_blank(1)
bbplot::axis_text(side = 1, line = 0.75)
bbplot::axis_text("Difference in beta diversity at sites that \nare gentrifying vs. those that are not",
                  side = 1, line = 3.5)


rect(
  xleft = mu_1[1],
  xright = mu_1[3],
  ybottom = 0.46,
  ytop = 1.54,
  col = "gray80",
  border = NA
)
lines(
  x = rep(mu_1[2],2),
  y = c(0.46,1.54),
  lty = 2
)


yloc <- seq(0.5,1.5, length.out = data_list$ncity)
for(i in 1:data_list$ncity){
  tt <- city_mu1[my_order[i],]
  lines(
    x = tt[-2],
    y = rep(yloc[i],2),
    lwd = 2
  )
  tmp_city <- city_map$pname[my_order[i]]
  text(
    x = city_mu1[my_order[i],3] + 0.0025,
    y = yloc[i],
    labels = substitute(bold(tmp_city), env = list(tmp_city = tmp_city)),
    pos = 4,
    cex = 0.8
  )
}
points(
  x = city_mu1[my_order,2],
  y = yloc,
  pch = 21,
  bg = "#00AADE",
  cex = 1.2
)

legend(
  "bottomright",
  legend = c(
    "E(sampled city)",
    "E(among city)",
    "95% among city CI",
    "95% sampled city CI"
  ),
  lty = c(NA, 2, NA, 1),
  lwd = c(NA, 2, NA, 2),
  pch = c(21, NA, 15, NA),
  pt.cex = c(1.3, 1.3, 2, NA),
  col = c(
    "black", "black",
    "gray80",
    "black"
  ),
  pt.bg = c(
    "#00AADE",
    NA, 
    "gray80",
    NA
  ),
  cex = 0.9,
  # box.col = "white",
  bty = "n"
)

dev.off()

city_pred <- vector(
  "list",
  length = ncity
)
for(city in 1:ncity){
  knots <- my_knots[
    my_knots$City == cities[city],
    c("min","median","max")]
  tmp_mcmc <- mc$beta_exp[,city,]

  city_pred[[city]] <- spline_pred_gradient(
    knots = knots,
    mcmc = tmp_mcmc
  )
  
}


my_cities <- c(
  'mela' = 12,
  'naca' = 13,
  'phaz' = 14
)

city_examples <- city_pred[my_cities]
windows(8,3)
{
  tiff(
    "./plots/beta_examples.tiff",
    height = 3,
    width = 8,
    units = "in",
    res = 1200,
    compression = "lzw"
  )
  par(mar = c(1,1,1,1), oma = c(4,4,0,0))
  m <- matrix(1:3, ncol = 3)
  ribbon_cols <- pals::brewer.seqseq2(9)[c(6,8)]
  layout(m)
  
  sub_names <- c(
    "A) Metro LA, CA",
    "B) Washington D.C.",
    "C) Phoenix, AZ"
  )
  for(i in 1:3){
    bbplot::blank(xlim = c(0,0.8), ylim = c(0.1,0.8), bty = "l")
    u <- par("usr")
    bbplot::axis_blank(1)
    bbplot::axis_blank(2)
    if(i == 1){
      bbplot::axis_text(side = 2, line = 0.7, las = 1)
    }
    text(
      x = u[1] + 0.04,
      y = u[4] - 0.04, pos = 4,
      labels = sub_names[i],
      cex = 1.5
    )
    bbplot::axis_text(side = 1, line = 0.7)
    if(i == 2){
      bbplot::axis_text(
        "Impervious cover (proportion)",
        1,
        outer = TRUE,
        at = NA,
        line = 2,
        cex = 1.5
      )
      bbplot::axis_text(
        "Beta diversity",
        2,
        outer = TRUE,
        at = NA,
        line = 2,
        cex = 1.5
      )
    }
    non_gent <- cbind(
      city_examples[[i]]$x[1:200],
      city_examples[[i]]$y[1:200,]
    )
    non_gent <- non_gent[non_gent[,1]<=0.8,]
    gent <- cbind(
      city_examples[[i]]$x[-c(1:200)],
      city_examples[[i]]$y[-c(1:200),]
    )
    gent <- gent[gent[,1] <=0.8,]
    bbplot::ribbon(
      x = gent[,1],
      y = gent[,c(2,4)],
      col = ribbon_cols[1],
      alpha = 0.5
    )
    bbplot::ribbon(
      x = non_gent[,1],
      y = non_gent[,c(2,4)],
      col = ribbon_cols[2],
      alpha = 0.5
    )
    lines(
      x = gent[,1],
      y = gent[,3],
      col = ribbon_cols[1],
      lwd = 3,
      lend = 2
    )
    lines(
      x = non_gent[,1],
      y = non_gent[,3],
      col = ribbon_cols[2],
      lwd = 3,
      lend = 2
    )
    
    
  }
  
  my_legend(
    x = 0.2,
    y = 0.35,
    legend = c("Gentrifying", "Non-gentrifying"),
    text.col = "white",
    fill = c(
      scales::alpha(ribbon_cols[1], 0.5),
      scales::alpha(ribbon_cols[2], 0.5)
    ),
    cex = 1.5,
    box.cex = c(1.25,1.25),
    y.intersp = 1.25,
    border = c(
      scales::alpha(ribbon_cols[1], 0.5),
      scales::alpha(ribbon_cols[2], 0.5)
    ),
    bty = "n"
  )
  my_legend(
    x = 0.2,
    y = 0.35,
    legend = c("Gentrified", "Not gentrified"),
    lty = c(1,1),
    lwd = 3,
    col = ribbon_cols,
    cex = 1.5,
    y.intersp = 1.25,
    bty = "n",
    seg.len = 1.25
  )
  dev.off()
}



city_pred <- vector("list", length = ncity)
for(city in 1:ncity){
  knots <- as.numeric(my_knots[
    my_knots$City == cities[city] &
      my_knots$covariate == "mean_19",
    c("min","median","max")])
  tmp_mcmc <- mc$beta_exp[,city,]
  
  city_pred[[city]] <- spline_pred(
    knots = knots,
    mcmc = tmp_mcmc[,c(5:7)]
  )
  
}

windows(12,8)
# try out a bivariate legend
svg("./plots/beta_impervious.svg",
    width = 12, height = 8)
{
  
  m <- matrix(
    c(1:5,24,6:23),
    ncol = 6,
    nrow = 4,
    byrow = TRUE
  )
  layout(m)
  par(mar = c(1,1,1,1), oma = c(6,6,0,0))
  
  
  # bbplot::blank(ylim = c(0,1.2), xlim = c(0,1), xaxs = "i", yaxs = "i")
  # rect(0.4,0.4,0.7,0.7, col = my_cols[3])
  # rect(0.7,0.7,1,1, col = my_cols[2])
  # rect(0.4,0.7,0.7,1, col = my_cols[4])
  # rect(0.7,0.4,1,0.7, col = my_cols[1])
  # text(x = 0.55, y = 0.3,  "N", adj = 0.5, cex = 2.25)
  # text(x = 0.85, y = 0.3,  "Y", adj = 0.5, cex = 2.25)
  # text(x = 0.3, y = 0.55, "N", adj = 0.5, cex = 2.25)
  # text(x = 0.3, y = 0.85, "Y", adj = 0.5, cex = 2.25)
  # text("Impervious", x = 0.67, y = 0.12, ad = 0.5, cex = 1.67)
  # text("Gentrification", x = 0.1, y = 0.62, ad = 0.5, cex = 1.67, srt = 90)
  
  biggest_diff <- rep(NA, 23)
  for(i in 1:23){
    biggest_diff[i] <- city_pred[[i]]$y[200,2] - city_pred[[i]]$y[1,2]
  }
  my_order <- order(biggest_diff, decreasing = TRUE)
  for(i in 1:24){
    
    if(i == 24){
      bbplot::blank()
      par(xpd = NA)
      legend(
        "center",
        legend = c(
          "Median estimate",
          "95% CI"
        ),
        lty = c(1,NA),
        pch = c(NA,22),
        pt.cex = c(NA, 4),
        pt.bg = c(
          NA,
          "#6495ED"
        ),
        lwd = c(3,NA),
        bty = "n",
        y.intersp = 1,
        cex = 1.2
      )
      bbplot::axis_text(
        "Impervious cover (proportion)",
        1,
        outer = TRUE,
        line = 3.25,
        at = NA,
        cex = 1.5
      )
      bbplot::axis_text(
        "f(Impervious cover)",
        2,
        outer = TRUE,
        line = 3.25,
        at = NA,
        cex = 1.5
      )
      next
    }
    my_city <- city_map$city[my_order[i]]
    
    
    bbplot::blank(ylim = c(0,1), xlim = c(0, 1), bty = "l",
                  main = city_map$pname[my_order[i]])
    bbplot::axis_blank(1)
    bbplot::axis_blank(2)
    if(i %in% m[,1]){
      bbplot::axis_text(side = 2, las = 1, line = 0.5)
    }
    if(i %in% m[4,]){
      bbplot::axis_text(side = 1, line = 0.75)
    }
    
    bbplot::ribbon(
      x = city_pred[[my_order[i]]]$x,
      y = city_pred[[my_order[i]]]$y[,-2],
      alpha = 0.8, col = "#6495ED"
    )
    
    
    lines(
      x = city_pred[[my_order[i]]]$x,
      y = city_pred[[my_order[i]]]$y[,2], lwd = 2
    )
    #lines(x = xx, y = city_gent[[i]][,2], lwd = 1, col = "white")
    
  }
  
  
  
  
}
dev.off()


plot(city_pred[[13]]$y[,2] ~ city_pred[[13]]$x, type = "l")

my_dissim <- data_list$dissim[dcovs$City.x == "naca"]
my_x <- dcovs$imp_diff[dcovs$City.x == "naca"]/100


test <- cut(
  my_x,
  breaks = 10
)
test2 <- split(
  my_dissim,
  test
)
test2 <- sapply(
  test2,
  median
)
my_points <- levels(test)
my_points <- gsub(
  "\\(|\\]",
  "",
  my_points
)
my_points <- strsplit(
  my_points,
  ","
)
my_points <- sapply(
  my_points,
  function(x) mean(as.numeric(x))
)

plot(test2 ~ my_points, type = "p",  bty = "l",
     ylab = "Dissimilarity",
     xlab = "Absolute difference in impervious cover between two sites",
     ylim = c(0,0.6),
     las = 1)


bbplot::ribbon(
  x = to_return$x,
  y = to_return$y[,-2],
  col = "red",
  alpha = 0.5
)

lines(
  x = to_return$x,
  y = to_return$y[,2],
  col = "red",
  lwd = 2

)




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
    my_splines = tmp_spline,
    mcmc_mat = tmp_mcmc,
    my_data = tmp_dissm
    #intercept = tmp_mcmc[,1]
  )
  
}
