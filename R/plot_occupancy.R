library(runjags)
library(dplyr)
library(bbplot)

# Prep the data for the model
species_map <- read.csv("./data/species_in_analysis.csv")

source("./R/alpha_beta_functions.R")

cat("loading in run.jags file...\n")
# Load in the occupancy model results
mout <- readRDS(
  "./results/occupancy_model_fit_simpler2.RDS"
)

cat("binding posterior simulations...\n")
# Compile the posterior
mcmc <- do.call(
  "rbind",
  mout$mcmc
)

# take a random sample to iterate through
set.seed(11556644)
mcsamp <- mcmc[sample(1:nrow(mcmc), 10000),]

rm(mout, mcmc)
gc()

mc <- split_mcmc(mcsamp)

gent_mu <- t(apply(
  mc$b_species[,2,],
  2,
  quantile,
  probs = c(0.05,0.5,0.95)
))

gent_imp <- t(apply(
  mc$b_species[,4,],
  2,
  quantile,
  probs = c(0.05,0.5,0.95)
))

row.names(gent_mu) <- species_map$species
row.names(gent_imp) <- species_map$species
gent_mu <- gent_mu[-4,]
gent_imp <- gent_imp[-4,]
gent_imp <- gent_imp[order(gent_imp[,2], decreasing = TRUE),]
gent_mu <- gent_mu[row.names(gent_imp),]


pretty_sp <- rev(c(
  "nine-banded armadillo",
  "gray squirrel sp.",
  "coyote",
  "red squirrel",
  "red fox",
  "eastern chipmunk",
  "badger",
  "striped skunk",
  "mule deer",
  "mink",
  "flying squirrel sp.",
  "beaver",
  "raccoon",
  "cottontail sp.",
  "elk",
  "white-tailed deer",
  "woodchuck",
  "Virginia opossum",
  "fox squirrel",
  "bobcat",
  "gray fox"
))

windows(6.5,8)

tiff(
  "./plots/occupancy_results.tiff",
  height = 8,
  width = 6.5,
  units = "in",
  compression = "lzw",
  res = 600
)

par(mar = c(4,9,0.5,6))
{
bbplot::blank(
  ylim = c(0.5,21.5),
  xlim = c(-4,4),
  xaxs = "i",
  yaxs = "i"
)
bbplot::axis_blank(1)
bbplot::axis_text(side = 1, line = 0.6)
u <- par("usr")
for(i in 1:21){
  if(i %% 2 != 0){
    next
  }
  rect(
    u[1],
    i - 0.5,
    u[2],
    i + 0.5,
    col = "gray80",
    border = "gray80"
    
  )
}

bbplot::axis_text(
  pretty_sp,
  at = 1:21,
  side = 2,
  line = 0.3,
  las = 2
)

test <- pals::coolwarm(4)
gent_colors <- c("#130A51", "#9d8ec3", "#e0827b", "#6C0015")
lines(
  x = c(0,0),
  y = u[3:4],
  lty = 2
)
for(i in 1:nrow(gent_mu)){
  tmp <- gent_mu[i,]
  if(length(unique(sign(tmp)))==1){
    tmp_col <- gent_colors[1]
  }else{
    tmp_col <- gent_colors[2]
  }
  rect(
    tmp[1],
    i + 0.15,
    tmp[3],
    i + 0.25,
    col = tmp_col,
    border = tmp_col
  )
  points(
    x = tmp[2],
    y = i + 0.2,
    pch = 16,
    col = tmp_col,
    cex = 1.5
  )
}

for(i in 1:nrow(gent_imp)){
  tmp <- gent_imp[i,]
  if(length(unique(sign(tmp)))==1){
    tmp_col <- gent_colors[4]
  }else{
    tmp_col <- gent_colors[3]
  }
  rect(
    tmp[1],
    i - 0.15,
    tmp[3],
    i - 0.25,
    col = tmp_col,
    border = tmp_col
  )
  points(
    x = tmp[2],
    y = i - 0.2,
    pch = 16,
    col = tmp_col,
    cex = 1.5
  )
}
}
par(xpd = NA)

legend(
  x = 4.1,
  y = mean(u[3:4]) + 7+3,
  c("Yes","No"),
  fill = gent_colors[1:2],
  title = "90% CI\nexcludes zero?\n\nGentrification",
  bty = "n"
)

legend(
  x = 4.25,
  y = mean(u[3:4]) + 2.7+3,
  c("Yes","No"),
  fill = rev(gent_colors[3:4]),
  title = "Gentrification\nX impervious",
  bty = "n"
)
lines(
  x = c(4.15,6.7),
  y = rep(mean(u[3:4]) + 5.8+3,2),
  lwd = 3
)
lines(
  x = c(4.3,6.6),
  y = rep(mean(u[3:4]) + 4.7+3,2)
)
lines(
  x = c(4.3,6.6),
  y = rep(mean(u[3:4])+1.6+3 ,2)
)
bbplot::axis_text(
  text = "Effect size (logit scale)",
  side = 1,
  line = 2,
  cex = 1.25
)
dev.off()
