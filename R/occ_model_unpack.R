mout <- readRDS(
  "./results/occupancy_model_fit_simpler_summary2.RDS"
)

source("./R/alpha_beta_functions.R")

mout <- readRDS(
  "./results/occupancy_model_fit_simpler2.RDS"
)

mc <- do.call("rbind", mout$mcmc)
mc <- split_mcmc(mc)
rm(mout)

# get gent
g1 <- t(apply(
  mc$b_species[,2,],
  2,
  function(x){
    round(
      quantile(
        x,
        probs = c(0.025,0.05,0.5,0.95,0.975)
      ),
      2
    )
  }
)
)

# get gent inxs
g2 <- t(apply(
  mc$b_species[,4,],
  2,
  function(x){
    round(
      quantile(
        x,
        probs = c(0.025,0.05,0.5,0.95,0.975)
      ),
      2
    )
  }
)
)



b_within <-
  apply(
    mc$b_within,
    2,
    function(x){
      round(
        quantile(
          x,
          probs = c(0.025,0.05,0.5,0.95,0.975)
        ),
        2
      )
    }
  )




sum(g1[,2]>0)


at_95 <-  unique(
  c(
    which(g1[,1]>0),
    which(g2[,5]<0)
  )
)
at_90 <-  unique(
  c(
    which(g1[,2]>0),
    which(g2[,4]<0)
  )
)
# order goes intercept, gent, imperv, gent * imperv
sp <- read.csv(
  "./data/species_in_analysis.csv"
)

sp[at_95,]
sp[at_90[-which(at_90%in%at_95)],]



# species specific regression coefficients
sc <- mout[
  grep(
    "b_species_city",
    row.names(mout)),1:3
]
# autologistic term
theta <- mout[
  grep(
    "theta\\[",
    row.names(mout)),2
]
# species present in city
xx <- mout[
  grep(
    "x\\[",
    row.names(mout)),2
]



cities <- c(
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


xx <- matrix(
  as.logical(xx),
  ncol = length(cities),
  nrow = nrow(sp)
)
source("../stats_help/skunk/R/mcmc_functions.R")
# get gent results for metro la


tmp_mat <- split_mcmc(t(sc))
theta  <- matrix(
  theta,
  ncol = length(cities),
  nrow = nrow(sp)
)

# read in beta diversity stuff
my_mcmc <- readRDS(
  "./mcmc_output/beta_output/beta_results_collapsed_norm_vegan.RDS"
)

my_mcmc <- do.call(
  "rbind",
  my_mcmc$mcmc
)
mc <- split_mcmc(my_mcmc)

beta_gent <- apply(
  mc$beta_exp[,,8],
  2,
  median
)
gent_order <- order(beta_gent, decreasing = TRUE)

gent_occ <- tmp_mat

# order everything appropriately
windows(30,15)
par(mfrow = c(4,6))
for(i in 1:23){

  bbplot::blank(xlim = c(-1,1), ylim = c(0,23))
  bbplot::axis_blank(1)
  bbplot::axis_text(side = 1, line = 0.5)
  bbplot::axis_text(cities[gent_order[i]], side = 3)
  my_lines <- tmp_mat$b_species_city[,2,,gent_order[i]]
  my_lines <- t(my_lines)
  my_lines <- my_lines[which(xx[,gent_order[i]]),]
  for(j in 1:nrow(my_lines)){
    lines(
      x = my_lines[j,-2],
      y = rep(j,2)
    )
  }
  points(
    x = my_lines[,2],
    y = 1:nrow(my_lines),
    pch = 19
  )
}
bbplot::blank(main = "gent effect")
windows(30,15)
par(mfrow = c(4,6))
for(i in 1:23){
  
  bbplot::blank(xlim = c(-6,2.75), ylim = c(-6,2.75))
  bbplot::axis_blank(1)
  bbplot::axis_text(side = 1, line = 0.5)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 2, line = 0.5, las = 1)
  bbplot::axis_text(cities[gent_order[i]], side = 3)
  my_lines <- tmp_mat$b_species_city[,1,,gent_order[i]]
  my_lines <- t(my_lines)
  my_lines <- my_lines[which(xx[,gent_order[i]]),]
  my_lines2 <- tmp_mat$b_species_city[,4,,gent_order[i]]
  my_lines2 <- t(my_lines2)
  my_lines2 <- my_lines2[which(xx[,gent_order[i]]),]
  points(
    y = my_lines[,2],
    x = my_lines2[,2],
    pch = 19
  )
}
bbplot::blank(main = "gent imperv interaction")
windows(30,15)
par(mfrow = c(4,6))

# MELA, raccoon occupancy increased by 10%,
# striped skunk occupancy went down 5%

for(i in 1:23){
  my_lines <- tmp_mat$b_species_city[2,c(1,3),,gent_order[i]]
  my_lines[2,] <- -0.5 * my_lines[2,]
  my_lines <- colSums(my_lines)
  my_psi <- plogis(my_lines[which(xx[,gent_order[i]])])
  my_psi2<- plogis(
    my_lines[which(xx[,gent_order[i]])] +
    theta[which(xx[,gent_order[i]]), gent_order[i]]
  )
  my_psi <- my_psi / (my_psi + (1 - my_psi2))
  my_gent_psi <- tmp_mat$b_species_city[2,1:4,,gent_order[i]]
  my_gent_psi <- t(my_gent_psi) %*% rbind(1, 1, -0.5, -0.5)
  my_gent_psi <- plogis(my_gent_psi[which(xx[,gent_order[i]])])
  my_gent_psi2<- tmp_mat$b_species_city[2,1:4,,gent_order[i]]
  my_gent_psi2 <- t(my_gent_psi2) %*% rbind(1,1,-0.5, -0.5)  
  my_gent_psi2 <- my_gent_psi2 + theta[,gent_order[i]]
  my_gent_psi2 <- plogis(my_gent_psi2[which(xx[,gent_order[i]])])
  my_gent_psi <- my_gent_psi / (my_gent_psi + (1 - my_gent_psi2))
  gmp <- my_gent_psi - my_psi
  bbplot::blank(xlim = c(-0.3,0.3), ylim = c(0,23))
  bbplot::axis_blank(1)
  bbplot::axis_text(side = 1, line = 0.5)
  bbplot::axis_text(cities[gent_order[i]], side = 3)
  points(x = gmp, y = 1:length(gmp))
  
}

plot(city_mu ~ beta_gent)

within_covs <- read.csv(
  "./data/cleaned_data/covariates/site_covariates.csv",
  stringsAsFactors = FALSE
)


aa <- within_covs %>% 
  dplyr::group_by(City) %>% 
  dplyr::summarise(
    mean_gent = round(mean(gentrifying),2),
    mean_vuln = round(mean(vulnerable),2)
  ) %>% 
  data.frame()

aa <- aa[order(aa$mean_gent, decreasing = TRUE),]
hm <- within_covs[within_covs$gentrifying,] %>% 
  dplyr::group_by(City) %>% 
  dplyr::summarise(
    imp = mean(mean_19)
  ) %>% 
  data.frame()
hm$imp <- hm$imp / 100
