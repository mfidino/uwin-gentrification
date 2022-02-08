###########################
#
# Alpha diversity analysis
#
# written by M. Fidino
#
###########################

library(sf)
library(runjags)
library(dplyr)
library(cli)
library(MCMCvis)

# we need to get the objects created to fit the MSOM,
# so let's get those data prepped again.

analysis <- "alpha_beta"

# the path to the data
data_path <- "../uwin-dataset/data_2021-07-27/cleaned_data/full_capture_history.csv"

# ALL OF THESE WILL BE REMOVED AT SOME POINT. 
# The city abbreviations
cities <- c(
  "ahga", "autx", "boma","buny","chil","deco","inin","ioio"#,
  #"jams","lrar","lbca","mawi","naca","phaz","phaz2","poor",
  #"rony","scut","sewa","slmo","tawa","toon","uril","wide"
)
# The spceies
species <- c(
  "cottontail_sp", "coyote", "fox_squirrel",# "gray_fox",
  "gray_squirrel_sp"#, "mule_deer", "raccoon", "red_fox",
  #"striped_skunk", "virginia_opossum", "white_tailed_deer", "woodchuck"
)

# the years of data to collect (using for now)
years <- c(
  "18|19|20|21"
)

# Test the autologistic indexing term is specified right. Works on
#  current dataset, but want to be able to test it with new data.
test_autologistic <- FALSE

# Prep data for the model
source("./R/format_data_for_msom.R")


# move on to alpha beta analysis

# read in the model
mod <- readRDS("./mcmc_output/trial_msom.RDS")

# convert to MCMC
mcmc <- do.call("rbind", mod$mcmc)

# get just the latent state
mcmc <- grab(mcmc, "^z")
# get just the 
nmcmc <- nrow(mcmc)

# beta diversity analysis
# This requires much more setup relative to the alpha diversity
#  analysis because we need to set up pairwise comparisons
#  among sites
all_clear <- rep(FALSE, nmcmc)
for(zi in 1:nmcmc){
  # for beta diversity we need to go wide with the z matrix,
  # but our current format is long.
  sp_dat$z <- mcmc[zi,]
  # going wide
  z_mat <- matrix(
    sp_dat$z,
    ncol = jags_list$nspecies,
    nrow = jags_list$ndata / jags_list$nspecies,
    byrow = TRUE
  )
  # add species names to z matrix
  colnames(z_mat) <- unique(sp_dat$Species)
  
  # need all the other data to go along with the z_matrix
  # subsetting down to the first species achieves this
  sp_dat_beta <- sp_dat[sp_dat$Species_id == 1,]
  sp_dat_beta <- sp_dat_beta[,-which(
    colnames(sp_dat_beta) %in% c(
      "City", "Species", "Season",
      "Crs","Y","J","Species_id", "Combo_id",
      "last_samplevec", "z"
    )
  )]
  
  # make the design matrix, this is unique
  #  to each sampling period and city,
  #  so we split the dataset that way
  sp_dat_split <- split(
    sp_dat_beta,
    factor(sp_dat_beta$Season_city_id)
  )
  
  # We will now go back from wide z_mat to long z_mat
  #  and also create the design matrix we need for the
  #  analysis. Because I've split everything into a 
  #  list of data.frames we will store these changes
  #  the same way (and then bind everything together
  #  afterwards).
  z_long <- vector("list", length = length(sp_dat_split))
  dm_long <- vector("list", length = length(sp_dat_split))
  
  for(i in 1:length(sp_dat_split)){
    # grab each data.frame on it's own
    tmp <- sp_dat_split[[i]]
    
    # as well as the associated part of the z matrix
    tmpz <- z_mat[
      which(sp_dat_beta$Season_city_id == i),
      ]
    
    # temporary design matrix 
    ####################################################
    ### UPDATE THIS WHEN YOU HAVE THE CORRECT COVARIATES
    # We'll want to group mean center before doing the 
    # spline matrix and whatnot.
    tmp_dm <- tmp[,-grep("_id", colnames(tmp))]
    ####################################################
    
    # make the spline matrix for each city season
    dm_long[[i]] <- make_spline_matrix(tmp_dm, "Site", "Long","Lat")
    # provide the season_city_id, which is just the iterator
    dm_long[[i]]$dmat_site_ids$Season_city_id <- i
    dm_long[[i]]$metadata$City_id <- unique(tmp$City_id)
    # make the z matrix long format
    z_long[[i]] <- as.data.frame(
      zmat_long(tmpz,dm_long[[i]]$dmat_site_ids)
    )
    # and again add the season_city_id
    z_long[[i]]$Season_city_id <- i
  }
  # get all of the metadata we need
  if(zi == 1){
  tmplist <- vector("list", length = length(dm_long))
  for(i in 1:length(tmplist)){
    tmplist[[i]] <- dm_long[[i]]$metadata
  }
  
  mdat <- do.call(
    "rbind",
    tmplist
  )  
  mdat <- mdat[!duplicated(mdat),]
  
  write.csv(
    mdat,
    "./mcmc_output/beta_output/knots.csv",
    row.names = FALSE
    )
  }
  # combine into one large data.frame
  z_long <- do.call(
    "rbind",
    z_long
  )
  # generate the beta diversity site a and b vectors
  tmp_beta <- do.call(
    "rbind",
    lapply(dm_long, function(x) x$dmat_site_ids)
  )
  # Next step is to get the correct site index across ALL
  #  locations. This was a bigger pain than I thought
  #  it would be. 
  asb <- data.frame(
    site = c(tmp_beta$siteA, tmp_beta$siteB),
    id = c(tmp_beta$siteA_id, tmp_beta$siteB_id)
  )
  # remove duplicates
  asb <- dplyr::distinct(
    asb
  )

  # tack on the City_id for each site
  asb <- dplyr::inner_join(
    asb,
    all_sites[,c("City_id","Site")],
    by = c("site" = "Site")
  )
  # sort by city, then by id
  asb <- asb[order(asb$City, asb$id),]
  # once again, split everything by city
  asb <- split(asb,  factor(asb$City))
  # The first city is fine, we need to add the last largest id value
  #  to the following data.frame, and continue that until we are done
  for(i in 2:length(asb)){
    last_val <- max(asb[[i-1]]$id)
    asb[[i]]$id <- asb[[i]]$id + last_val
  }
  # And recombine
  asb <- do.call(
    "rbind",
    asb
  )
  # Now we can rejoin back onto the tmp_beta data.frame to 
  #  create the updated site_id's that can be used across
  #  city_season id's. Need to do this individually for 
  #  site A and site B
  tmp_a <- dplyr::inner_join(
    tmp_beta[,c("siteA","Season_city_id")],
    asb[,c("site","id","City_id")],
    by = c("siteA" = "site")
  )
  tmp_b <- dplyr::inner_join(
    tmp_beta[,c("siteB","Season_city_id")],
    asb[,c("site","id", "City_id")],
    by = c("siteB" = "site")
  )
  
  
  # Make a new data.frame that has the updated ids
  dmat_site_ids <- data.frame(
    siteA = tmp_a$siteA,
    siteA_id = tmp_a$id,
    siteB = tmp_b$siteB,
    siteB_id = tmp_b$id,
    Season_city_id = tmp_a$Season_city_id,
    City_id = tmp_a$City_id
  )
  
  
  # create the data list
  data_list <- list(
    # response variable for binomial regression
    betaz = cbind(z_long$V1, z_long$V2),
    # design matrix of slope terms
    design_matrix_beta = do.call(
      "rbind",
      lapply(
        dm_long, 
        "[[",
        1
      )
    ),
    # random intercept for city
    city_vec_beta = dmat_site_ids$City_id,
    # site A in comparison ID
    siteA_id = dmat_site_ids$siteA_id,
    # site B in comparison ID
    siteB_id = dmat_site_ids$siteB_id,
    # number of data points
    ndata_beta = nrow(dmat_site_ids),
    # number of slope terms
    npar_beta = ncol(dm_long[[1]]$dmat),
    # number of cities
    ncity = max(sp_dat$City_id),
    # number of sites
    nsite = length(unique(sp_dat$Site))
  )
  # function to generate initial values
  inits <- function(chain){
    gen_list <- function(chain = chain){
      list( 
        beta_mu = rnorm(data_list$npar_beta),
        beta_tau = rgamma(data_list$npar_beta,1,1),
        b0_mu = rnorm(1),
        b0_tau = rgamma(1,1,1),
        b0 = rnorm(data_list$ncity),
        beta_site_tau = rgamma(1,1,1),
        beta_site_re = rnorm(data_list$nsite),
        beta_log = matrix(
          rnorm(
            data_list$npar_beta * data_list$ncity
          ),
          nrow = data_list$ncity,
          ncol = data_list$npar_beta
        ),
        .RNG.name = switch(chain,
                           "1" = "base::Wichmann-Hill",
                           "2" = "base::Marsaglia-Multicarry",
                           "3" = "base::Super-Duper",
                           "4" = "base::Mersenne-Twister",
                           "5" = "base::Wichmann-Hill",
                           "6" = "base::Marsaglia-Multicarry",
                           "7" = "base::Super-Duper",
                           "8" = "base::Mersenne-Twister"),
        .RNG.seed = sample(1:1e+06, 1)
      )
    }
    return(switch(chain,           
                  "1" = gen_list(chain),
                  "2" = gen_list(chain),
                  "3" = gen_list(chain),
                  "4" = gen_list(chain),
                  "5" = gen_list(chain),
                  "6" = gen_list(chain),
                  "7" = gen_list(chain),
                  "8" = gen_list(chain)
    )
    )
  }
  # fit model
  m1 <- runjags::run.jags(
    "./JAGS/impute_beta.R",
    monitor = c(
      "b0", "b0_mu", "b0_sd", "beta_exp",
      "beta_mu","beta_tau", "beta_site_re",
      "beta_site_sd"
    ),
    data = data_list,
    inits = inits,
    burnin = 5500,
    sample = 750,
    n.chains = 3,
    thin = 2,
    adapt = 1000,
    modules = "glm",
    method= "parallel"
  )
  msum <- summary(m1)
  
  # save model summary
  write.csv(
    msum,
    paste0(
      "./mcmc_output/beta_output/model_summarys/msum_",
      stringr::str_pad(
        zi,
        width = 4,
        pad = "0"
      ),
      ".csv"
    ),
    row.names = FALSE
  )
  # Check if model converged
  if(all(msum[,colnames(msum)=="psrf"] <1.1)){
    # if it did give the all_clear!
    all_clear[zi] <- TRUE
    # and save the model output
    write.csv(
      do.call("rbind", m1$mcmc),
      paste0(
        "./mcmc_output/beta_output/mcmc/posterior_",
        stringr::str_pad(
          zi,
          width = 4,
          pad = "0"
        ),
        ".csv"
      ),
      row.names = FALSE
    )
  }
  
}

if(all(all_clear)){
  print("Model converged, move forward!")
} else {
  stop("Some models did not converge. Investigate!")
}
# compile all the results
all_mcmc <- lapply(
  list.files(
    "./mcmc_output/beta_output/mcmc/",
    full.names = TRUE
  ),
  read.csv
)

# combine it
all_mcmc <- do.call("rbind", all_mcmc)

write.csv(
  all_mcmc,
  "./mcmc_output/beta_output/beta_mcmc.csv",
  row.names = FALSE
)
