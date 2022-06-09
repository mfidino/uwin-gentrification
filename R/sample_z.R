


nmcmc <- nrow(mcsamp_piece$a_among)
nspecies <- data_list$nspecies
ncity <- data_list$ncity
nsamples_one <- data_list$nsamples_one
nsamples_two <- data_list$nsamples_two


z <- psi <- rho <- numerator <- denominator <- z_prob <- matrix(
  NA,
  ncol = length(data_list$y),
  nrow = nmcmc
)

psitheta <- matrix(
  NA,
  ncol = data_list$nsamples_one,
  nrow = nmcmc
)

# among-city regression (latent-state)
city_prob <- array(NA, dim = c(nmcmc, nspecies, ncity))
for(j in 1:nspecies){
  for(k in 1:ncity){
    if(sp_observed[j,k] == 1){
      city_prob[,j,k] <- 1
    } else {
      city_prob[,j,k] <- plogis(mcsamp_piece$b_among[,,j] %*%data_list$among_covs[j,k,])
    }
  }
}
# within-city regression (latent-state) for 'first' surveys
cat("\nsimulating z first 'season'\n(has 2 progress bars)...\n")
cat("Progress bar 1 of 2\n")
pb <- txtProgressBar(max= nsamples_one)
for(s in 1:nsamples_one){
  setTxtProgressBar(pb,s)
  # if observed then we say species is present
  if(data_list$y[s] > 0){
    z[,s] <- 1
    next
  }
  # Other calculate probability of occupancy 1
  psi[,s] <-  mcsamp_piece$b_species_city[
    , , data_list$species_idx[s], data_list$city_idx[s]] %*%
    data_list$psi_covs[s,]
  psitheta[,s] <- psi[,s] + 
    mcsamp_piece$theta[, data_list$species_idx[s],
                         data_list$city_idx[s]]
  # and daily probability of detection
  rho[,s] <- mcsamp_piece$c_species_city[
    , , data_list$species_idx[s], data_list$city_idx[s]] %*%
    data_list$rho_covs[s,]
}
# convert to probability
psi[,1:nsamples_one] <- plogis(psi[,1:nsamples_one])
psitheta <- 1 - plogis(psitheta)
# equilibrium occupancy at t = 1
psi[,1:nsamples_one] <- psi[,1:nsamples_one] / (
  psi[,1:nsamples_one] + psitheta
)
rho[,1:nsamples_one] <- plogis(rho[,1:nsamples_one])
psi[,which(data_list$y>0)] <- 1
rho[,which(data_list$y>0)] <- 1
# back into the loop to index the species and city
cat("\nProgress bar 2 of 2\n")
for(s in 1:nsamples_one){
  setTxtProgressBar(pb,s)
  if(data_list$y[s]>0){
    z_prob[,s] <- 1
    next
  }
  # Prob you were present and not detected
  numerator[,s] <- city_prob[,
                             data_list$species_idx[s],
                             data_list$city_id[s]] *
    psi[,s] *
    (1 - rho[,s])^data_list$J[s]
  # Either you are not there or you were there
  #  and not detected.
  denominator[,s] <- (
    city_prob[,
              data_list$species_idx[s],
              data_list$city_id[s]] * 
      (1 - psi[,s])
  ) + numerator[,s]
  z_prob[,s] <- numerator[,s] / denominator[,s]
}
# sample z
z[,1:nsamples_one] <- rbinom(
  nsamples_one * nmcmc,
  1,
  z_prob[,1:nsamples_one]
)
z[,which(data_list$y>0)[1:nsamples_one]] <- 1

# and then with the auto-regressive term added in
# within-city regression (latent-state) with auto-logistic term.
# at this point, however, we need to start calculating the probability
# for each sample
cat("\nsampling remaining seasons\n(has one progress bar)...\n")
pb <- txtProgressBar(min =nsamples_one+1, max = nsamples_two)
for(s in (nsamples_one+1):nsamples_two){
  setTxtProgressBar(pb,s)
  # if observed then we say species is present
  if(data_list$y[s] > 0){
    z[,s] <- 1
    psi[,s] <- 1
    rho[,s] <- 1
  }else{
  # Other calculate probability of occupancy
  psi[,s] <-  mcsamp_piece$b_species_city[
    , , data_list$species_idx[s], data_list$city_idx[s]] %*%
    data_list$psi_covs[s,] + 
    mcsamp_piece$theta[,
                 data_list$species_idx[s],
                 data_list$city_idx[s]
    ] * z[,data_list$last_sample_vec[s]]
  psi[,s] <- plogis(psi[,s])
  # and daily probability of detection
  rho[,s] <- mcsamp_piece$c_species_city[
    , , data_list$species_idx[s], data_list$city_idx[s]] %*%
    data_list$rho_covs[s,]
  rho[,s] <- plogis(rho[,s])
  
  # Prob you were present and not detected
  numerator[,s] <- city_prob[,
                             data_list$species_idx[s],
                             data_list$city_id[s]] *
    psi[,s] *
    (1 - rho[,s])^data_list$J[s]
  # Either you are not there or you were there
  #  and not detected.
  denominator[,s] <- (
    city_prob[,
              data_list$species_idx[s],
              data_list$city_id[s]] * 
      (1 - psi[,s])
  ) + numerator[,s]
  z_prob[,s] <- numerator[,s] / denominator[,s]
  z[,s] <- rbinom(nmcmc,1, z_prob[,s])
  }
}
if(analysis == "alpha"){
  # calculate species richness at each site. The tmp_covs object
  #   stores the site, city, and season of data collection. It
  #   is the same length as z and ordered the same way.
  #site_info$z <- z
  cat("\ncalculating species richness...\n")
  pb <- txtProgressBar(max = nmcmc)
  for(i in 1:nmcmc){
    setTxtProgressBar(pb, i)
    site_info$z <- z[i,]  
    
    
    sp_rich <- site_info %>% 
      dplyr::group_by(Site, City, Season) %>% 
      dplyr::summarise(
        rich = sum(z),
        gentrifying = all(gentrifying),
        mean_19 = unique(mean_19),
        .groups = "drop_last"
      ) %>% 
      data.frame()
    
    sp_rich$Season <- factor(
      sp_rich$Season,
      levels = order_seasons(as.character(sp_rich$Season))
    )
    
    
    sp_rich <- sp_rich[order(sp_rich$City, sp_rich$Season, sp_rich$Site),]
    sp_rich_mcmc[mc_loc[i],] <- sp_rich$rich
  }
  rm(denominator, numerator, z, z_prob, psi, rho)
  gc()
}
if(analysis == "beta"){
  rm(psi, rho, denominator, numerator, z_prob)
  gc()
  # make unique site, as there may be some sites that
  #  share names among cities.
  sp_dat$citysite <- paste0(sp_dat$City,"-",sp_dat$Site)
  cat("\ncalculating pairwise dissimilarity...\n")
  pb <- txtProgressBar(max = nrow(z))
  for(zi in 1:nrow(z)){
    setTxtProgressBar(pb, zi)
    # for beta diversity we need to go wide with the z matrix,
    # but our current format is long.
    sp_dat$z <- z[zi,]
    
    # going wide
    # need to go wide, but we have to split by city 
    #  and season first.
    sp_dat_list <- split(
     sp_dat,
     factor(
       paste0(
         sp_dat$City,"-",sp_dat$Season
       ),
       levels = unique(
         paste0(
           sp_dat$City,"-",sp_dat$Season
         )
       )
     )
    )
    # store z and the rest of the beta info for 
    #  modeling
    zlist <- sp_dat_beta <- 
      vector("list", length = length(sp_dat_list))
  
    # now make a z matrix for each sampling period.
    for(j in 1:length(sp_dat_list)){
      # reorder sp_dat_list
      sp_dat_list[[j]] <- sp_dat_list[[j]][
        order(
          sp_dat_list[[j]]$Site, 
          sp_dat_list[[j]]$Species
        ),
      ]
       zlist[[j]] <- matrix(
         sp_dat_list[[j]]$z,
         ncol = data_list$nspecies,
         nrow = length(
           unique(
             sp_dat_list[[j]]$Site
           )
         ),
         byrow = TRUE
       )
       # Add relevant column/ rownames
       colnames(zlist[[j]]) <- unique(
         sp_dat_list[[j]]$Species
       )
       rownames(zlist[[j]]) <- unique(
         sp_dat_list[[j]]$Site
       )
       # need all the other data to go along with the z_matrix
       #   subsetting down to the first species achieves this
       #   from sp_dat_list
       sp_dat_beta[[j]] <- sp_dat_list[[j]][sp_dat_list[[j]]$Species_id == 1,]
       sp_dat_beta[[j]] <- sp_dat_beta[[j]][,-which(
         colnames(sp_dat_beta[[j]]) %in% c(
           "City", "Species", "Season",
           "Crs","Y","J","Species_id", "Combo_id",
           "last_samplevec", "z", "city", "Site"
         )
       )]
    }
    names(sp_dat_beta) <- names(sp_dat_list)
    names(zlist) <- names(sp_dat_list)
    
    # We will now go back from wide z_mat to long z_mat
    #  and also create the design matrix we need for the
    #  analysis. Because I've split everything into a 
    #  list of data.frames we will store these changes
    #  the same way (and then bind everything together
    #  afterwards).
    z_long <- vector("list", length = length(sp_dat_beta))
    dm_long <- vector("list", length = length(sp_dat_beta))
    
    for(i in 1:length(sp_dat_beta)){
      # grab each data.frame on it's own
      tmpsds <- sp_dat_beta[[i]]
      
      # as well as the associated part of the z matrix
      tmpz <- zlist[[i]]
      
      # temporary design matrix 
      ####################################################
      ### UPDATE THIS WHEN YOU HAVE THE CORRECT COVARIATES
      # We'll want to group mean center before doing the 
      # spline matrix and whatnot.
      tmp_dm <- tmpsds[,-grep("_id|starter|rowID|last_sample_vec|X", colnames(tmpsds))]
      ####################################################
      
      # make the spline matrix for each city season
      dm_long[[i]] <- make_spline_matrix(tmp_dm, "citysite", "Long","Lat")
      # provide the season_city_id, which is just the iterator
      #dm_long[[i]]$dmat_site_ids$Season_city_id <- i
      #dm_long[[i]]$metadata$City_id <- unique(tmp$City_id)
      # make the z matrix long format
      z_long[[i]] <- as.data.frame(
        zmat_long(tmpz,dm_long[[i]]$dmat_site_ids)
      )
      # and again add the season_city_id
     # z_long[[i]]$Season_city_id <- i
    }
    # get all of the metadata we need
    if(zi == 1 & gr == 1){
      tmplist <- vector("list", length = length(dm_long))
      for(i in 1:length(tmplist)){
        tmplist[[i]] <- dm_long[[i]]$metadata
      }
      
      mdat <- do.call(
        "rbind",
        tmplist
      )  
      mdat$city_season <- names(sp_dat_list)
      #mdat <- mdat[!duplicated(mdat),]
      
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
    beta_results[[mc_loc[zi]]] <- z_long
    # generate the beta diversity site a and b vectors

    if(zi == 1 & gr == 1){
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
      # plus site duplicates
      asb <- asb[!duplicated(asb$site),]
      asb <- asb[order(asb$site, asb$id),]
      # tack on the City_id for each site
      asb <- dplyr::inner_join(
        asb,
        dplyr::distinct(sp_dat[,c("City_id","citysite")]),
        by = c("site" = "citysite")
      )
      # sort by city, then by id
      asb <- asb[order(asb$City, asb$id),]
      # once again, split everything by city
      asb <- split(asb,  factor(asb$City))
      # The first city is fine, we need to add the last largest id value
      #  to the following data.frame, and continue that until we are done.
      #  On top of this, since we dropped duplicate sites (sampling
      #  across multiple seasons), then we need to rewrite the
      #  id values before adding things to them
      for(i in 2:length(asb)){
        last_val <- max(asb[[i-1]]$id)
        asb[[i]]$id <- (1:nrow(asb[[i]])) + last_val
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
        data.frame(siteA= tmp_beta[,c("siteA")]),
        asb[,c("site","id","City_id")],
        by = c("siteA" = "site")
      )
      tmp_b <- dplyr::inner_join(
        data.frame(siteB = tmp_beta[,c("siteB")]),
        asb[,c("site","id", "City_id")],
        by = c("siteB" = "site")
      )
      # Make a new data.frame that has the updated ids
      dmat_site_ids <- data.frame(
        siteA = tmp_a$siteA,
        siteA_id = tmp_a$id,
        siteB = tmp_b$siteB,
        siteB_id = tmp_b$id,
        City_id = tmp_a$City_id
      )
    }
  }
  rm(z)
}

