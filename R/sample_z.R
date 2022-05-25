


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
  # Other calculate probability of occupancy
  psi[,s] <-  mcsamp_piece$b_species_city[
    , , data_list$species_idx[s], data_list$city_idx[s]] %*%
    data_list$psi_covs[s,]
  # and daily probability of detection
  rho[,s] <- mcsamp_piece$c_species_city[
    , , data_list$species_idx[s], data_list$city_idx[s]] %*%
    data_list$rho_covs[s,]
}
# convert to probability
psi[,1:nsamples_one] <- plogis(psi[,1:nsamples_one])
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
    next
  }
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
      (1 - psi[s])
  ) + numerator[,s]
  z_prob[,s] <- numerator[,s] / denominator[,s]
  z[,s] <- rbinom(nmcmc,1, z_prob[,s])
}
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
  
  sp_rich$Season <- order_seasons(sp_rich$Season)
  
  
  sp_rich <- sp_rich[order(sp_rich$City, sp_rich$Season, sp_rich$Site),]
  sp_rich_mcmc[mc_loc[i],] <- sp_rich$rich
}
rm(denominator, numerator, z, z_prob, psi, rho)
gc()

