###########################
#
# Alpha diversity analysis
#
# written by M. Fidino
#
###########################

library(sf)
library(run.jags)
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


sp_dat$Site_id <- as.numeric(factor(sp_dat$Site))
# Alpha diversity analysis
all_clear <- rep(FALSE, nmcmc)
for(zi in 1:nmcmc){
  print(zi)

  # add on the estimate to the species data
  sp_dat$z <- mcmc[zi,]
  
  # get alpha diversity
  alpha <- sp_dat %>% 
    dplyr::group_by(Season_city_id, Site) %>% 
    dplyr::summarise(
      rich = sum(z),
      Season_id = unique(Season_id),
      City_id = unique(City_id),
      Site_id = unique(Site_id),
      Long = unique(Long),
      Lat = unique(Lat),
      .groups = "drop_last"
    ) %>% 
    as.data.frame()

    # adjust to whatever we are using (need to update this)
    design_matrix_alpha <- matrix(
      1,
      ncol = 1,
      nrow = nrow(alpha)
    )
    # commenting off for simpler model
    #tmp <- data.frame(
    #  City = sp_dat$City,
    #  Site = sp_dat$Site
    #)
    #tmp <- dplyr::distinct(tmp)
    #tmp <- dplyr::inner_join(
    #  tmp,
    #  city_map,
    #  by = "City"
    #)
    data_list <- list(
      alphaz = alpha$rich,
      npar_alpha = ncol(design_matrix_alpha),
      ndata_alpha = nrow(alpha),
      ncity = length(unique(alpha$City_id)),
      nsite = length(unique(alpha$Site_id)),
      design_matrix_alpha = design_matrix_alpha,
      site_vec_alpha = alpha$Site_id,
      city_vec_alpha = alpha$City_id#,
     # site_vec = tmp$City_id
    )
    # initial values function
    inits <- function(chain){
      gen_list <- function(chain = chain){
        list( 
          alpha_site_re = rnorm(data_list$nsite),
          #alpha_site_shape = rgamma(1,1,1),
          #alpha_site_rate = rgamma(1,1,1),
          alpha_mu = rnorm(data_list$npar_alpha),
          alpha_tau = rgamma(data_list$npar_alpha,1,1),
          alpha = matrix(
            rnorm(
              data_list$npar_alpha * data_list$ncity
            ),
            nrow = data_list$ncity,
            ncol = data_list$npar_alpha
          ),
          alpha_site_tau = rgamma(1,1,1),
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
    "./JAGS/impute_alpha.R",
    monitor = c(
      "alpha", "alpha_mu", "alpha_sd", "alpha_site_sd"
    ),
    data = data_list,
    inits = inits,
    burnin = 4000,
    sample = 1000,
    n.chains = 3,
    thin = 2,
    adapt = 1000,
    modules = "glm",
    method= "rjparallel"
  )
  # summarise model
  msum <- round(summary(m1),2)

  # save model summary
  write.csv(
    msum,
    paste0(
      "./mcmc_output/alpha_output/model_summarys/msum_",
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
        "./mcmc_output/alpha_output/mcmc/posterior_",
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
    "./mcmc_output/alpha_output/mcmc/",
    full.names = TRUE
  ),
  read.csv
)

# combine it
all_mcmc <- do.call("rbind", all_mcmc)

write.csv(
  all_mcmc,
  "./mcmc_output/alpha_output/alpha_mcmc.csv",
  row.names = FALSE
)



# beta diversity analysis
for(zi in 1:nmcmc){
  # to do beta diversity we need to go wide with the z matrix.
  
  z_mat <- matrix(
    sp_dat$z,
    ncol = jags_list$nspecies,
    nrow = jags_list$ndata / jags_list$nspecies,
    byrow = TRUE
  )
  colnames(z_mat) <- unique(sp_dat$Species)
  sp_dat_beta <- sp_dat[sp_dat$Species_id == 1,]
  sp_dat_beta <- sp_dat_beta[,-which(
    colnames(sp_dat_beta) %in% c(
      "City", "Species", "Season",
      "Crs","Y","J","Species_id", "Combo_id",
      "last_samplevec", "z"
    )
  )]
  
  # make the design matrix
  sp_dat_split <- split(
    sp_dat_beta,
    factor(sp_dat_beta$Season_city_id)
  )
  
  z_long <- vector("list", length = length(sp_dat_split))
  dm_long <- vector("list", length = length(sp_dat_split))
  
  for(i in 1:length(sp_dat_split)){
    tmp <- sp_dat_split[[i]]
    
    tmpz <- z_mat[
      which(sp_dat_beta$Season_city_id == i),
      ]
    tmp_dm <- tmp[,-grep("_id", colnames(tmp))]
    
    dm_long[[i]] <- make_spline_matrix(tmp_dm, "Site", "Long","Lat")
    dm_long[[i]]$dmat_site_ids$Season_city_id <- i
    z_long[[i]] <- as.data.frame(zmat_long(tmpz,dm_long[[i]]$dmat_site_ids))
    z_long[[i]]$Season_city_id <- i
  }
  z_long <- dplyr::bind_rows(z_long)
  # generate the beta diversity site and b vectors
  tmp_beta <- dplyr::bind_rows(
    lapply(dm_long, function(x) x$dmat_site_ids)
  )
  # get all the sites and id's in one vector
  asb <- data.frame(
    site = c(tmp_beta$siteA, tmp_beta$siteB),
    id = c(tmp_beta$siteA_id, tmp_beta$siteB_id)
  )
  asb <- dplyr::distinct(asb)
  asb <- dplyr::inner_join(
    asb,
    all_sites[,c("City","Site")],
    by = c("site" = "Site")
  )
  asb <- asb[order(asb$City, asb$id),]
  asb <- split(asb,  factor(asb$City))
  for(i in 2:length(asb)){
    last_val <- max(asb[[i-1]]$id)
    asb[[i]]$id <- asb[[i]]$id + last_val
  }
  asb <- dplyr::bind_rows(asb)
  test_a <- dplyr::inner_join(
    tmp_beta[,c("siteA","Season_city_id")],
    asb[,c("site","id")],
    by = c("siteA" = "site")
  )
  test_b <- dplyr::inner_join(
    tmp_beta[,c("siteB","Season_city_id")],
    asb[,c("site","id")],
    by = c("siteB" = "site")
  )
  
}