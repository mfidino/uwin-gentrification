library(pander)
library(dplyr)
library(bbplot)
library(viridis)
library(raster)
library(sf)
library(runjags)
library(prettymapr)
source("./R/census_functions.R")
source("./R/gdm_functions.R")
source("./R/alpha_beta_functions.R")
sf::sf_use_s2(FALSE)


gent_list  <- readRDS(
  "./data/census_data/gent_sites.rds"
)
gent_list_poly <- readRDS(
  "./data/census_data/gent_census_tracts.rds"
)

imp_census <- read.csv(
  "./data/cleaned_data/covariates/census_tract_imperv_cover.csv"
)
imp_census <- split(
  imp_census,
  factor(imp_census$city)
)
pretty_names <- c(
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
  "Portland, Oregon",
  "Rochester, NY",
  "Sanford, FL",
  "Salt Lake City, UT",
  "Seattle, WA",
  "Saint Louis, MO",
  "Tacoma, WA",
  "Urbana, IL",
  "Wilmington, DE"
)
ribbon_cols <- pals::brewer.seqseq2(9)[c(6,8)]

# read in the mcmc stuff too and do predictions
my_knots <- read.csv(
  "./mcmc_output/beta_output/knots.csv"
)

# read in the mcmc
my_mcmc <- readRDS(
  "./mcmc_output/beta_output/beta_results_collapsed_norm_vegan.RDS"
)
my_mcmc <- do.call(
  "rbind",
  my_mcmc$mcmc
)
mc <- split_mcmc(my_mcmc)

city_examples <- vector(
  "list",
  length = length(gent_list_poly)
)
cities <- unique(my_knots$City)
for(city in 1:length(city_examples)){
  knots <- my_knots[
    my_knots$City == cities[city],
    c("min","median","max")]
  tmp_mcmc <- mc$beta_exp[,city,]
  
  city_examples[[city]] <- spline_pred_gradient(
    knots = knots,
    mcmc = tmp_mcmc
  )
  
}

# do the same for alpha diversity
sp_rich <- read.csv(
  "./results/alpha_for_stage_two_collapsed.csv"
)


# covariates
my_site <- read.csv(
  "./data/cleaned_data/covariates/site_covariates.csv"
)

# divide by 100
my_site$imp <- my_site$mean_19 / 100

# join to data
sp_rich <-  dplyr::inner_join(
  sp_rich,
  my_site[,c("City", "Site", "imp")],
  by = c("City", "Site")
)
# for plotting
city_map <- data.frame(
  city = unique(sp_rich$City)
)

my_site <- split(
  my_site,
  factor(my_site$City)
)
# pretty name
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



my_mc <- readRDS(
  "./mcmc_output/alpha_output/alpha_mcmc.RDS"
)

my_mc <- do.call(
  "rbind",
  my_mc$mcmc
)

my_mc <- split_mcmc(
  my_mc
)
xx <- seq(0, 0.8, length.out = 300)

alpha_examples <- vector(
  "list",
  length = length(gent_list_poly)
)
for(i in 1:length(gent_list_poly)){
  pp <-  my_mc$alpha[,i,] %*% t(cbind(1,0,xx,0))
  pp <- exp(pp)
  
  pp2 <-  my_mc$alpha[,i,] %*% t(cbind(1,1,xx,xx))
  pp2 <- exp(pp2)
  
  alpha_examples[[i]]$gent <- apply(
    pp2,
    2,
    qwrap
  )
  alpha_examples[[i]]$non_gent <- apply(
    pp,
    2,
    qwrap
  )
}
roundUpNice <- function(x, nice=c(1,2,3,4,5,6,7,8,9,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

gamma_samples <- read.csv(
  "./results/city_gamma_samples.csv"
)

windows(12, 10)


pb <- txtProgressBar(max = length(gent_list_poly))
for(i in 1:length(gent_list_poly)){
  setTxtProgressBar(pb, i)
  pdf(
    paste0(
      "./supp_plots/supp_1/",cities[i], ".pdf"
    ),
    height = 10,
    width = 12
  )
  # Figure out all the vulnerable census tracts across a city
  m <- matrix(
    c(1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,
      2,2,2,2,2,2,2,2,
      3,3,3,3,4,4,4,4,
      3,3,3,3,4,4,4,4,
      5,5,5,5,6,6,6,6,
      5,5,5,5,6,6,6,6), ncol = 8, nrow = 8, byrow = TRUE)
  layout(m)
  par(oma = c(0,0,5,0))
  # get impervious cover range from this city
  my_range <- range(city_examples[[i]]$x) * 100
  my_range <- roundUpNice(my_range[2]) / 100
  plot(
    gent_list_poly[[i]]["gentrifying"],
    main = "",
    reset = FALSE,
    pal = scales::alpha(rev(ribbon_cols),0.5),
    key.pos = NULL,
    mar = c(4,4,4,4)
  )
  bbplot::axis_text(
    text = paste0(
      pretty_names[i]
    ),
    outer = TRUE,
    at = NA,
    cex = 2
  )
  
  my_bgs <- ifelse(
    gent_list[[i]]$gentrifying,
    "black",
    "white"
  )
  points(
    st_coordinates(
      gent_list[[i]]
    ),
    pch =21,
    cex = 2,
    bg = my_bgs,
  )
  prettymapr::addnortharrow(pos = "topleft", padin = c(1.1,0.19), scale = 0.7)
  prettymapr::addscalebar(plotepsg = 4326, style = "ticks",
                          padin =c(0.25, 0.00), lwd = 2, label.cex = 1)
  par(mar = c(0,0,0,0), xpd = NA)
  blank(xlim = c(1,10), ylim = c(1,10), xaxs = "i")
  
  points(
    x = 2.5,
    y = 7.5,
    pch = 15,
    col = scales::alpha(ribbon_cols,0.5)[1],
    cex = 5
  )
  text(x = 2.65, y = 7.5, labels ="gentrified", pos = 4, cex = 1.8)
  points(
    x = 5.5,
    y = 7.5,
    pch = 15,
    col = scales::alpha(ribbon_cols,0.5)[2],
    cex = 5
  )
  text(x = 5.65, y = 7.5, labels ="not gentrified", pos = 4, cex = 1.8)
  
  points(
    x = 2.5,
    y = 3.5,
    pch = 21,
    bg = "black",
    cex = 5
  )
  text(x = 2.65, y = 3.5, labels ="site < 500m of gentrified census tract", pos = 4, cex = 1.8)
  points(
    x = 5.5,
    y = 3.5,
    pch = 21,
    bg = "white",
    cex = 5
  )
  text(x = 5.65, y = 3.5, labels ="site > 500m of gentrified census tract", pos = 4, cex = 1.8)
  
  par(mar = c(6,6,0,1), xpd = TRUE)
  # alpha diversity plot
  alpha_ylim <- max(
    c(
      alpha_examples[[i]]$gent[3,],
      alpha_examples[[i]]$non_gent[3,]
    )
  )
  alpha_ylim <- ceiling(alpha_ylim)
  bbplot::blank(
    xlim = c(0,my_range),
    ylim = c(0,alpha_ylim),
    bty = "l"
  )
  
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  
  bbplot::axis_text(side = 2, line = 0.7, las = 1)
  
  bbplot::axis_text(side = 1, line = 0.7)
  bbplot::axis_text(
    "Impervious cover (proportion)",
    1,
    line = 3.5,
    cex = 1.5
  )
  bbplot::axis_text(
    "Species richness",
    2,
    line = 3.5,
    cex = 1.5
  )
  to_keep <- which(xx <= max(city_examples[[i]]$x))
  bbplot::ribbon(
    x = xx[to_keep],
    y = t(alpha_examples[[i]]$gent)[to_keep,-2],
    col = ribbon_cols[1],
    alpha = 0.5
  )
  bbplot::ribbon(
    x = xx[to_keep],
    y = t(alpha_examples[[i]]$non_gent)[to_keep,-2],
    col = ribbon_cols[2],
    alpha = 0.5
  )
  lines(
    x = xx[to_keep],
    y = alpha_examples[[i]]$gent[2,to_keep],
    col = ribbon_cols[1],
    lwd = 3
  )
  lines(
    x = xx[to_keep],
    y = alpha_examples[[i]]$non_gent[2,to_keep],
    col = ribbon_cols[2],
    lwd = 3
  )
  
  beta_ylim <- max(
    city_examples[[i]]$y[,3]
  )
  beta_ylim <-roundUpNice(beta_ylim * 100)/100
  bbplot::blank(
    xlim = c(0,my_range),
    ylim = c(0,beta_ylim),
    bty = "l"
  )
  
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, line = 0.7, las = 1)
  bbplot::axis_text(side = 1, line = 0.7)
  bbplot::axis_text(
    "Impervious cover (proportion)",
    1,
    line = 3.5,
    cex = 1.5
  )
  bbplot::axis_text(
    "Beta diversity",
    2,
    at = NA,
    line = 3.5,
    cex = 1.5
  )
  
  
  non_gent <- cbind(
    city_examples[[i]]$x[1:200],
    city_examples[[i]]$y[1:200,]
  )
  gent <- cbind(
    city_examples[[i]]$x[-c(1:200)],
    city_examples[[i]]$y[-c(1:200),]
  )
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
    lwd = 3
  )
  lines(
    x = non_gent[,1],
    y = non_gent[,3],
    col = ribbon_cols[2],
    lwd = 3
  )
  #    my_legend(
  #      "topleft",
  #   legend = c("Gentrifying", "Non-gentrifying"),
  #   text.col = "white",
  #   fill = c(
  #     scales::alpha(ribbon_cols[1], 0.5),
  #     scales::alpha(ribbon_cols[2], 0.5)
  #   ),
  #   cex = 1.5,
  #   box.cex = c(1.25,1.25),
  #   y.intersp = 1.25,
  #   border = c(
  #     scales::alpha(ribbon_cols[1], 0.5),
  #     scales::alpha(ribbon_cols[2], 0.5)
  #   ),
  #   bty = "n"
  # )
  # my_legend(
  #   "topleft",
  #   legend = c("Gentrifying", "Non-gentrifying"),
  #   lty = c(1,1),
  #   lwd = 3,
  #   col = ribbon_cols,
  #   cex = 1.5,
  #   y.intersp = 1.25,
  #   bty = "n",
  #   seg.len = 1.25
  # )
  hist_data <- my_site[[i]]
  hist_data$mean_19 <- hist_data$mean_19 / 100
  
  yo1 <- hist(
    hist_data$mean_19[
      hist_data$gent == 1
    ],
    plot = FALSE
  )
  yo2 <- hist(
    hist_data$mean_19[
      hist_data$gent == 0
    ],
    plot = FALSE
  )
  y_range <- c(0, max(yo1$counts, yo2$counts))
  x_range <- c(0, max(yo1$breaks, yo2$breaks))
  u <- par("usr")
  bbplot::blank(
    ylim = y_range,
    xlim = x_range,
    bty = "l"
  )
  plot(yo2, add = TRUE, col = scales::alpha(ribbon_cols[2], 0.5))
  plot(yo1, add = TRUE, col = scales::alpha(ribbon_cols[1], 0.5))
  
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, line = 0.7, las = 1)
  bbplot::axis_text(side = 1, line = 0.7)
  bbplot::axis_text(
    "Impervious cover (proportion)",
    1,
    line = 3.5,
    cex = 1.5
  )
  bbplot::axis_text(
    "Site frequency",
    2,
    at = NA,
    line = 3.5,
    cex = 1.5
  )
  rich_samples <- gamma_samples[,i]
  rich_samples <- prop.table(table(rich_samples))
  bbplot::blank(
    xlim = c(0,21),
    ylim = c(0,1),
    bty ="l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, line = 0.7, las = 1)
  bbplot::axis_text(side = 1, line = 0.7)
  bbplot::axis_text(
    "Overall richness",
    1,
    line = 3.5,
    cex = 1.5
  )
  bbplot::axis_text(
    "Probability",
    2,
    at = NA,
    line = 3.5,
    cex = 1.5
  )
  rich_plot <- data.frame(
    x = 0:21,
    y = 0
  )
  rich_plot$y[rich_plot$x %in% as.numeric(names(rich_samples))] <- rich_samples
  lines(
    x = rich_plot$x,
    y = rich_plot$y,
    lwd = 3,
    col = "gray40"
  )
  points(
    x = rich_plot$x,
    y = rich_plot$y,
    pch = 19,
    col  = "white",
    cex = 2
  )
  points(
    x = rich_plot$x,
    y = rich_plot$y,
    pch = 19,
    col  = "gray40",
    cex = 1.6
  )
  points(
    x = rich_plot$x,
    y = rich_plot$y,
    pch = 19,
    col = "black",
    cex = 1.2
  )
  
  dev.off()
}
