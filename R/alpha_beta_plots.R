library(sf)
sf::sf_use_s2(FALSE)
library(runjags)
library(dplyr)
library(cli)
library(MCMCvis)
library(pals)

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
  "./mcmc_output/beta_output/beta_results_collapsed_norm_vegan.RDS"
)

data_list <- readRDS(
  "./mcmc_output/beta_output/data_list.RDS"
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

tiff(
  "./plots/alpha_beta_results.tiff",
  height = 6,
  width = 4,
  units = "in",
  res = 1200,
  compression = "lzw"
)
m <- matrix(
  1:2,
  ncol = 1,
  nrow = 2
)
layout(m)
par(mar = c(2,1,1,1), oma = c(2,6,0,0), xpd = NA)
{
set.seed(2222)
bbplot::blank( xlim = c(-1, 1), ylim = c(0.25,3.25), xaxs = "i", yaxs = "i")
bbplot::axis_blank(1, tck = -0.03,
                   at = seq(-1, 1, 0.5))
bbplot::axis_text(side = 1, line = 0.35,
                  text = seq(-1,1, 0.5),
                  at = seq(-1, 1, 0.5))
lines(
  x = c(0,0),
  y = c(0.25,3.25),
  lty = 2
)
bbplot::axis_text(
  c(
    "Gentrification",
    "Impervious  ",
    "Impervious x\nGentrification"
  ),
  2,
  line = 0.1,
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
bbplot::blank( xlim = c(0, 0.4), ylim = c(0.25,3.25), xaxs = "i", yaxs = "i")
bbplot::axis_blank(1, tck = -0.03)
bbplot::axis_text(side = 1, line = 0.35)

bbplot::axis_text(
  c(
    "Gentrification",
    "Impervious",
    "Geographic\ndistance  "
  ),
  2,
  line = 0.1,
  at = rev(1:3),
  las = 1
)
po <- rev(1:3)
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
  my_jitter <- runif(data_list$ncity, -0.3,0.3)
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
      x[po[i]+1,4]
    }
  )
  points(
    x = tmp_vals,
    y = j + my_jitter,
    bg = scales::alpha("#00AADE", 0.7),
    pch = 21
  )
  
  rect(
    xleft = overall_mean[po[i]+1,2],
    xright = overall_mean[po[i]+1,4],
    ybottom = j - 0.45,
    ytop = j + 0.45,
    col = scales::alpha("gray40",0.5),
    border = NA
  )
  lines(
    x = rep(overall_mean[po[i]+1,3],2),
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

}

dev.off()
# p <- seq(0.001,0.999, length.out = 20)
# 
# cl <- function(x) -log(1-x)
# 
# pp <- cl(p)
# 
# icl <- function(x) 1 - exp(-x)
# 
# pp2 <- icl(pp)

# map of all the study cities and try to show their overall alpha and beta effects?

counties <- sf::read_sf(
  "D:/GIS/counties/cb_2020_us_county_500k"
)

# Only keep contiguous US
counties <- counties[-which(
  counties$STATE_NAME %in% c("Hawaii", "Alaska", "Puerto Rico",
                             "Guam", "Commonwealth of the Northern Mariana Islands",
                             "American Samoa", "United States Virgin Islands")),]
# split by state name
counties <- split(
  counties,
  factor(counties$STATE_NAME)
)

# join all counties for each state
counties <- lapply(
  counties,
  sf::st_union
)
# change back to sf objects
counties <- lapply(
  counties,
  sf::st_as_sf
)

# combine back to one data set
counties <- dplyr::bind_rows(
  counties
)

# read in site coordinates
coords <- read.csv(
  "./data/gentri_all_coordinates.csv"
)

coords <- coords %>% 
  dplyr::group_by(City) %>% 
  dplyr::summarise(
    Long = mean(Long),
    Lat = mean(Lat)
  )

coords <- sf::st_as_sf(
  coords,
  coords = c("Long","Lat"),
  crs = 4326
)
coords <- sf::st_coordinates(
  coords
)

# get E(increase) in species richness, alpha diversity
e_alpha <- exp(alpha_city[2,,2])

# and calculate difference in beta diversity in gentrifying vs non-gentrifying
city_mu1 <- matrix(
  ncol = 3,
  nrow = data_list$ncity
)


for(i in 1:data_list$ncity){
  tmp_mu_1 <-(1 - exp(-beta_mc$beta_exp[,i,1]))
  tmp_mu_2 <-( 1 - exp(-rowSums(beta_mc$beta_exp[,i,c(1,8)])))
  city_mu1[i,] <- quantile(tmp_mu_2 - tmp_mu_1, probs = c(0.025,0.5,0.975))
}

# combine both for plotting purposes
bb <- data.frame(
  alpha = e_alpha,
  beta = city_mu1[,2]
)
# adding city names to switch seattle and tacoma for plotting
bb$city_names <- c(
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

# break each response into groups of 3
bb <- bb %>% 
  dplyr::mutate(
    alpha_cat = cut(alpha, breaks = 3),
    beta_cat = cut(beta, breaks = 3)
  )

# number of categories
ncat <- length(
  levels(
    bb$alpha_cat
  )
)

# names of the different scales (to match
#  colors)
scale_names <- paste0(
  rep(1:ncat, each = ncat), " - ", rep(1:ncat, ncat)
)

# the colors
scale_vals <- pals::brewer.seqseq2()

# make group names to match with scale_names
bb$group <- paste0(
  as.numeric(bb$alpha_cat)," - ", as.numeric(bb$beta_cat)
)
# give them the scale_names as levels
bb$group <- factor(
  bb$group,
  levels = scale_names
)

# convert to numeric to match to colors
bb$color <- scale_vals[
  as.numeric(bb$group)
]


# function to get aspect ratio of a plot,
#  using to make a square legend
get.asp <-
  function(){
    pin <- par('pin')
    usr <- par('usr')
    asp <- (pin[2]/(usr[4]-usr[3])) / (pin[1]/(usr[2]-usr[1]))
    return(asp)
  }

tiff(
  "./plots/gentrification_map.tiff",
  height = 6,
  width = 8,
  units = "in",
  res = 1200,
  compression = "lzw"
)
par(mar = c(2,1,0,0), xpd = NA)
# plot out counties real quick
bb$x <- coords[,1]
bb$y <- coords[,2]
# flip tacoma and seattle so seattle is on top
tmp <- bb[
  grep("Seattle", bb$city_names),
]
bb[
  grep("Seattle", bb$city_names),
] <- bb[grep("Tacoma", bb$city_names),]
bb[
  grep("Tacoma", bb$city_names)[2],
] <- tmp

{
plot(counties, lwd = 2, reset = FALSE)
u <- par("usr")
# get aspect ratio (ratio of width to height)

my_asp <- get.asp()
my_x <- seq(
  -120,
  -113,
  length.out = ncat
)
my_diff <- my_x[length(my_x)] - my_x[1]
my_y <- seq(22, 22 * my_asp, length.out = ncat)

image(
  my_x,
  my_y,
  z =matrix(1:9, nrow=3),
  add = TRUE,
  col = pals::brewer.seqseq2()
)
text(
  x = mean(my_x),
  y = max(my_y) + 3.75,
  labels = "Relative strength of gentrification on...",
  pos = 1
)
text(
  x = min(my_x) - 4.5,
  y = mean(my_y),
  labels = "Alpha diversity",
  srt = 90
)
arrows(
  x0 = min(my_x)-2.8,
  y0 = min(my_y) - 1,
  y1 = max(my_y) + 1,
  length = 0.1,
  angle = 30,
  lwd = 1.5
)

text(
  x = mean(my_x),
  y = min(my_y) - 3.5,
  labels = "Beta diversity"
)
arrows(
  x0 = min(my_x)-1.6,
  x1 = max(my_x) + 1.6,
  y0 = min(my_y) - 2.5,
  length = 0.1,
  angle = 30,
  lwd = 1.5
)
points(
  bb$x,
  bb$y,
  pch = 21,
  cex = 2,
  bg = bb$color
)

}
dev.off()


