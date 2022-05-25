#####################################################
#
# Simulate latent state from occupancy model results
#
#
#
####################################################


library(runjags)
library(dplyr)

# Prep the data for the model
source("./R/prep_data_occupancy.R")

cat("loading in run.jags file...\n")
# Load in the occupancy model results
mout <- readRDS(
  "./results/occupancy_model_fit.RDS"
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
# function to split up mcsamp into a list with correctly shaped arrays
split_mcmc <- function(x){
  # get parameter names
  pars <- colnames(x)
  # unique parameters
  unq_pars <- unique(
    gsub(
      "\\[.*\\]",
      "",
      pars
    )
  )
  # make list object to store arrays in
  result_list <- vector(
    "list",
    length = length(unq_pars)
  )
  names(result_list) <- unq_pars
  # fill in the arrays
  for(i in 1:length(result_list)){
    # get just the parameters
    tmp <- pars[grep(
      paste0(
        "^",unq_pars[i], "\\["
      ),
      pars
    )]
    if(length(tmp) == 0){
      tmp <- pars[grep(
          paste0("^",unq_pars[i],"$"),
        pars
      )]
    }
    # and then the array dimensions
    arr_dim <- gsub(
      paste0(
        unq_pars[i],"\\[|\\]"
      ),
      "",
      tmp
    )
    arr_dim <- strsplit(
      arr_dim,
      ","
    )
    ndim <- length(arr_dim[[1]])
    npar <- length(arr_dim)
    # make a matrix
    arr_ind <- suppressWarnings(
      matrix(
       as.numeric(
         unlist(arr_dim)
        ),
        ncol = ndim,
        nrow = npar,
        byrow = TRUE
      )
    )
    if(nrow(arr_ind) == 1 & ncol(arr_ind) == 1){
      arr_ind[1,1] <- 1
    }
    # get max index for each
    max_ind <- apply(arr_ind, 2, max)
    # and then fill in the array
    result_list[[i]] <- array(
      x[,tmp],
      dim = c(nrow(x), max_ind)
    )
    
  }
  return(result_list)
}
# make pieces of this sample because we cannot iterate through
#  the whole thing
mcsamp_list <- vector(
  "list",
  length = 10
)

my_groups <- rep(1:10, each = 1000)
ngroup <- 10
for(i in 1:10){
  mcsamp_list[[i]] <- split_mcmc(
    mcsamp[which(my_groups == i),]
  )
}
# and then just get all of it too for some other
#  calculations
mcsamp <- split_mcmc(
  mcsamp
)

# simulate the probability of:
#  1) Species presence in a city
#  2) Species occupancy at a site
#  3) Species probability of detection

# determine if species was detected in a city. This will
#  modify our z matrix posterior calculation.
sp_observed <- apply(
  mcsamp$x,
  c(2,3),
  sum
)
sp_observed[sp_observed < nrow(mcsamp$x)] <- 0
sp_observed[sp_observed > 0] <- 1

# the z vector for one mcmc sample, plus the other objects
#  that are of similar length. The object tmp_covs has all
#  the site info and is the same length. 
site_info <- tmp_covs

# figure out how many unique site, city, seasons we have
unq_site_samps <- site_info[,1:3]
unq_site_samps <- unq_site_samps[!duplicated(unq_site_samps),]

unq_site_samps$Season <- order_seasons(unq_site_samps$Season)

# order by city, season, and then site
unq_site_samps <- unq_site_samps[
  order(unq_site_samps$City, unq_site_samps$Season, unq_site_samps$Site),
]

# and now make a matrix to store the sp_rich results
sp_rich_mcmc <- matrix(
  NA,
  ncol = nrow(unq_site_samps),
  nrow = nrow(mcsamp$a_among)
)

# and now we need to split the mcmc into pieces because we cannot store
#  all of the results at once. Splitting into pieces of 1K.

for(gr in 1:ngroup){
  cat(
    paste0("group ", gr, " of ", ngroup,"...\n")
  )
  mcsamp_piece <- mcsamp_list[[gr]]
  # where to store stuff in sp_rich_mcmc
  mc_loc <- which(my_groups == gr)
  
  source("./R/sample_z.R")
}


# calculate mean and sd
sp_rich_mean <- apply(sp_rich_mcmc, 2, mean)
sp_rich_sd <- apply(sp_rich_mcmc, 2, sd)

sp_rich$mu <- sp_rich_mean
sp_rich$sd <- sp_rich_sd
sp_rich <- sp_rich[,-which(colnames(sp_rich) == "rich")]
write.csv(
  sp_rich,
  "./results/alpha_for_stage_two.csv",
  row.names = FALSE
)




