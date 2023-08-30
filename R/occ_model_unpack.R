library(dplyr)
library(data.table)

# read in the model output

mc <- data.table::fread(
  "./results/occupancy_mcmc.csv",
  data.table = FALSE
)
mc <- as.matrix(mc)

source("./R/alpha_beta_functions.R")


mc <- split_mcmc(mc)

species <- read.csv(
  "./data/species_in_analysis.csv"
)
cities <-  c(
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

occ_tabs <- vector("list", length = length(mc) - 1)
names(occ_tabs) <- names(mc)[1:length(occ_tabs)]
for(i in 1:length(occ_tabs)){
  tmp <- mc[[i]]
  td <- dim(tmp)
  if(length(td) == 2){
    occ_tabs[[i]] <- apply(
      tmp,
      2,
      quantile,
      probs = c(0.025,0.05,0.5,0.95,0.975)
    )
    occ_tabs[[i]] <- t(
      round(
        occ_tabs[[i]],
        2
      )
    )
  }
  if(length(td) == 3){
    occ_tabs[[i]] <- apply(
      tmp,
      c(2,3),
      quantile,
      probs = c(0.025,0.05,0.5,0.95,0.975)
    )
    occ_tabs[[i]] <- round(
      occ_tabs[[i]],
      2
    )
  }
  if(length(td) == 4){
    occ_tabs[[i]] <- apply(
      tmp,
      c(2,3,4),
      quantile,
      probs = c(0.025,0.05,0.5,0.95,0.975)
    )
    occ_tabs[[i]] <- round(
      occ_tabs[[i]],
      2
    )
  }
  
}

saveRDS(occ_tabs, "./supplemental/occupancy_summary.RDS")


