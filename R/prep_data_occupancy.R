

# This script assumes you have ran the following scripts, which write
#  data into the "./data/cleaned_data" folder:
#  "./R/cleaning_scripts/data_cleaning_script.R"
#  "./R/cleaning_scripts/covariate_cleaning_script.R"

cat("Loading functions...\n")
# load functions used to clean data
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)

for(fn in functions_to_load){
  source(fn)
}

cat("Reading in detection data...\n")
# read in the data
city_data <- read.csv(
  "./data/detection_data.csv",
  stringsAsFactors = FALSE
)

# drop J01-NOR1
city_data <- city_data[-which(city_data$Site == "J01-NOR1"),]

# There are A LOT of NA's in this sampling array 
#  prop.table(table(city_data$J > 0))
#  FALSE  TRUE 
#  0.7   0.3

# Because of this we do not want to iterate over the observation's with 
#  no data as it greatly increases the amount of time it takes to fit
#  the model. We'll need to be able to index the species, city, and season.
#  To do this, it's going to be easiest to analyze these data in long-format.

# Additionally, we're going to only analyze a few species that were detected
#  a minimum of 10 times in the dataset. This will remove a number of rare
#  species, 

reduce_species <- table(
  city_data$Species,
  city_data$Y > 0
)
#reduce_species[order(reduce_species[,2]),]\

min_dets <- 25

first_species <- reduce_species[which(reduce_species[,2] >= min_dets),]

cdet <- city_data %>% group_by(Species, City) %>% 
  summarise(cdet = any(Y>0), .groups = "drop_last") %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  summarise(n = sum(cdet), .groups = "drop_last") %>% 
  data.frame()

cdet <- cdet[order(cdet$n),]

to_go <- cdet$Species[cdet$n <= 3]

first_species <- first_species[-which(
  row.names(first_species) %in% to_go),]

#write.csv(
#  data.frame(
#    species = row.names(first_species)
#  ),
#  "./data/species_in_analysis.csv",
#  row.names = FALSE
#)

if(exists("my_species")){
  if(length(my_species)>0){
    first_species <- first_species[
      which(row.names(first_species) %in% my_species),
      ]
  }
}

# get just these species for now
tmp <- city_data[
  which(city_data$Species %in% row.names(first_species)),
  ]


# these data.frames will help us connect the model output back to 
#  the actual species, seasons, or cities they are associated to. These
#  data.frames will also give us the numeric index we will need for
#  the JAGS model.
species_map <- data.frame(
  Species = sort(
    unique(
      tmp$Species
    )
  )
)

species_map$Species_id <- make_index(
  species_map$Species
)

season_map <- data.frame(
  Season = unique(
    tmp$Season
  )
)

season_map$Season_id <- make_index(
  season_map$Season,
  TRUE
)

city_map <- data.frame(
  City = sort(
    unique(
      tmp$City
    )
  )
)

city_map$City_id <- make_index(
  city_map$City
)

# make pretty names for cities (for plotting).
# pcity <- data.frame(
#   City = city_map$City,
#   Pretty = c(
#     "Atlanta,\nGeorgia",
#     "Austin,\nTexas",
#     "Chicago,\nIllinois",
#     "Denver,\nColorado",
#     "Edmonton,\nAlberta",
#     "Fort Collins,\nColorado",
#     "Iowa City,\nIowa",
#     "Indianapolis,\nIndiana",
#     "Jackson,\nMississippi",
#     "Manhattan,\nKansas",
#     "Madison,\nWisconsin",
#     "Metro LA,\nCalifornia",
#     "Phoenix,\nArizona",
#     "Rochester,\nNew York",
#     "Sanford,\nFlorida",
#     "Salt Lake\nCity, Utah",
#     "Seattle,\nWashington",
#     "Saint Louis,\nMissouri",
#     "Tacoma,\nWashington",
#     "Wilmington,\nDelaware"
#   )
# )

# now that we have these we can make a temporary city_data and join the maps.
tmp <- dplyr::left_join(
  tmp,
  species_map,
  by = "Species"
) %>% 
  dplyr::left_join(
    .,
    season_map,
    by = "Season"
  ) %>% 
  dplyr::left_join(
    .,
    city_map,
    by = "City"
  )

cat("Removing sites with no data...\n")
# now we just need to drop the rows where J == 0 (i.e., no sampling).
tmp2 <- tmp[tmp$J > 0,]

tmp2 <- split(
  tmp2,
  factor(
    tmp2$Season,
    order_seasons(tmp2$Season)
  )
)


for(i in 1:length(tmp2)){
  tmp2[[i]] <- table(
    paste0(tmp2[[i]]$City, "_", tmp2[[i]]$Site)
  )
  tmp2[[i]] <- tmp2[[i]][tmp2[[i]] != nrow(species_map)]
}

# remove these as well
tmp <- tmp[tmp$J > 0,]
for(i in 1:length(tmp2)){
  if(length(tmp2[[i]])==0) next
  to_go <- data.frame(
    City = sapply( strsplit(names(tmp2[[i]]),"_"), "[[",1),
    Site = sapply( strsplit(names(tmp2[[i]]),"_"), "[[",2),
    Season = names(tmp2)[i]
  )
  for(j in 1:nrow(to_go)){
    baddies <- which(
      tmp$Season == to_go$Season[j] &
        tmp$Site == to_go$Site[j] &
        tmp$City == to_go$City[j]
    )
    if(length(baddies)> 0){
      tmp <- tmp[-baddies,]
    } else {
      stop()
    }
  }
  
}


# We have multiple seasons of sampling and therefore need to account for 
#  this in the model. The issue here though is that each city does not
#  have data for each season. As such, we need to be able to index
#  the correct seasonal parameter in our long format detection array.

season_density <- t(
  table(
    tmp$Season_id,
    tmp$City_id
  )
)

# give row names that are the city names
row.names(season_density) <- city_map$City

# make this binary, 1 means that city has data on a given season.
season_density[season_density>0] <- 1

# this tells us which city and season has data.
season_has_data <- which(
  season_density ==1,
  arr.ind = TRUE
)
# give useful headers to 'season_has_data'
colnames(season_has_data) <- c(
  "city_id",
  "season_id"
)
# order by season then city (i..e, same way as species detection data).
season_has_data <- season_has_data[
  order(
    season_has_data[,"season_id"],
    season_has_data[,"city_id"]
  ),
  ]

# season_has_data is our first step, it's effectively a map that let's us
#  index what cities and seasons have data. We do model the probability
#  of each species in a city, but the seasons and cities that have data
#  vary. As such we are going to make a unique identifier for each
#  species, season, and city. 
#  It's important to note that tmp is sorted by species, season, then city.
#  Cities that do not have data on a given season are removed at this point
#  in the script as well. As such, we need to make a unique identifier
#  for every combination of species, season, and city.
combo <- tmp[
  ,
  c("Species_id", "Season_id", "City_id")
  ]
# remove duplicates, this will tell us how many parameters we are attempting
#  to estimate (for the seasonal stuff).
combo <- combo[
  -which(
    duplicated(combo)
  ),
  ]
# make the unique identifier
combo$Combo_id <- 1:nrow(combo)

# and then join this back to the species detection array
tmp <- left_join(tmp,
                 combo,
                 by = c("Species_id", "Season_id", "City_id")
)

# pull in range_data, which we use to estimate if a species is 
#  within a given city. This is for the top level part of the model
#  which helps determine the species pool of a given city (along with
#  the spatial coordinates of a city).
city_covs <- read.csv(
  "./data/cleaned_data/covariates/dist_city_to_species_edge.csv",
  row.names = 1
)

# reduce down to the species we can model, order it in the 
#  same way as species map (i.e., alphabetically), and then convert to 
#  a matrix.
city_covs <- city_covs[which(row.names(city_covs) %in% species_map$Species),]
city_covs <- city_covs[order(row.names(city_covs)),]
city_covs <- as.matrix(city_covs)


# Calculate average lat and long of each city
city_ll <- tmp[,c("City", "Lat", "Long")]
city_ll <- city_ll[-which(duplicated(city_ll)),]

city_ll <- city_ll %>% 
  dplyr::group_by(City) %>% 
  dplyr::summarise(
    Lat = mean(Lat),
    Long = mean(Long),
    .groups = "drop_last"
  )

# make sure city_ll is ordered by city
city_ll <- city_ll[order(city_ll$City),]

# make it a matrix
city_ll <- as.matrix(city_ll[,c("Lat", "Long")])

# among_covs will be input into JAGS. This is for all of the among-city
#  variables for the top-level of the model. The first level of the array
#  is all 1's for the intercept.
among_covs <- array(
  1,
  dim = c(dim(city_covs),2)
)
# add distance of species to city, and scale by standard deviation.
#  we are not subtracting the mean because we want the 0 to indicate
#  that a city is on the edge of a species range.
among_covs[,,2] <- city_covs / sd(city_covs)
# and then add the latitude 
# among_covs[,,3] <- matrix(
#   scale(city_ll[,"Lat"]),
#   ncol = ncol(city_covs),
#   nrow = nrow(city_covs),
#   byrow = TRUE
# )
# # and longitude
# among_covs[,,4] <- matrix(
#   scale(city_ll[,"Long"]),
#   ncol = ncol(city_covs),
#   nrow = nrow(city_covs),
#   byrow = TRUE
# )

cat("Creating auto-logistic id for model...\n")
# set up the model for the autologistic formulation. First off,
# get just the site, city, and season info.
scs <- dplyr::distinct(tmp[,c("Site","City","Season")])
scs$si <- paste0(scs$City,"-",scs$Site)

scs <- split(
  scs,
  factor(scs$City)
)
all_seas <- unique(tmp$Season)

first_sites <- vector("list", length = length(scs))
for(i in 1:length(scs)){
  my_seas <- unique(scs[[i]]$Season)
  my_seas <- all_seas[
    which(all_seas == my_seas[1]):
      which(all_seas == my_seas[length(my_seas)])
  ]
  
  first_sites[[i]] <- data.frame(
    sites = scs[[i]]$Site[scs[[i]]$Season == my_seas[1]],
    season = my_seas[1],
    city = unique(scs[[i]]$City)
  )
  for(j in 2:length(my_seas)){
    a1 <- scs[[i]]$Site[scs[[i]]$Season == my_seas[j]]
    if(length(a1) == 0) next
    a2 <- scs[[i]]$Site[scs[[i]]$Season == my_seas[j-1]]
    tmp_si <- which(
      !a1 %in% a2
    )
    if(length(tmp_si)>0){
      first_sites[[i]] <- rbind(
        first_sites[[i]],
        data.frame(
          sites = scs[[i]]$Site[scs[[i]]$Season == my_seas[j]][tmp_si],
          season = my_seas[j],
          city = unique(scs[[i]]$City)
        )
      )
    }
  }
}
# these are the "first" instances of sampling at a location
# that may need to be linked to sit
first_sites <- do.call("rbind", first_sites)
first_sites$season <- factor(
  first_sites$season,
  order_seasons(first_sites$season)
)
first_sites <- first_sites[order(first_sites$city, first_sites$season, first_sites$sites),]
row.names(first_sites) <- NULL
# reorder data so these are all at the very beginning
my_rows <- vector("list", length = nrow(first_sites))
for(i in 1:nrow(first_sites)){
  my_rows[[i]] <-  which(
    tmp$Site == first_sites$sites[i] &
      tmp$City == first_sites$city[i]  &
      tmp$Season == first_sites$season[i]
  )
}
  
new_tmp <- vector("list", length = length(my_rows))
for(i in 1:length(my_rows)){
  new_tmp[[i]] <- tmp[my_rows[[i]],]
}
new_tmp <- do.call("rbind", new_tmp)
new_tmp$starter <- TRUE
tmp$starter <- FALSE
tmp <- rbind(
  new_tmp,
  tmp[-which(tmp$X %in% new_tmp$X),]
)
row.names(tmp) <- NULL

# make a unique identifier for each sample
tmp$rowID <- 1:nrow(tmp)
  
sp_recs <- distinct(tmp[,c("Species","Site","City")])
sp_recs_list <- vector("list", length = nrow(sp_recs))
tmp$Season <- factor(tmp$Season, order_seasons(tmp$Season))
# vector that incidates where the previous sample was
tmp$last_sample_vec <- NA
sp_recs[sp_recs$Species == "woodchuck" & sp_recs$Site == "MLP1",]
  
pb <- txtProgressBar(max = nrow(sp_recs))
for(i in 1:nrow(sp_recs)){
  setTxtProgressBar(pb, i)
  small_dat <- tmp[
    tmp$Species == sp_recs$Species[i] &
    tmp$Site == sp_recs$Site[i] &
    tmp$City == sp_recs$City[i],
  ]
  small_dat <- small_dat[order(small_dat$Season),]
  if(all(diff(small_dat$Season_id) == 1)){
    tmp$last_sample_vec[small_dat$rowID[-1]] <- 
      small_dat$rowID[1:(nrow(small_dat)-1)]
  } else {
    s_groups <- which(small_dat$starter)
    # check if any of the s_groups only have 0 trailing seasons
    
    for(j in 1:length(s_groups)){
      if(j < length(s_groups)){
        tdat <- small_dat[s_groups[j]:(s_groups[j+1]-1),]
      } else {
        tdat <- small_dat[s_groups[j]:nrow(small_dat),]
      }
      if(nrow(tdat) == 1) next
      if(all(diff(tdat$Season_id) == 1)){
        tmp$last_sample_vec[tdat$rowID[-1]] <- 
          tdat$rowID[1:(nrow(tdat)-1)]
      } else {
        stop("investigate")
      }
    }
  }
}
  # just do a quick check to make sure this was set up correctly.
my_seas <- levels(tmp$Season)
for(i in rev(1:nrow(tmp))){
  if(tmp$starter[i]) next
  my_eval <- tmp$Species[tmp$last_sample_vec[i]] == tmp$Species[i] &
    tmp$Site[tmp$last_sample_vec[i]] == tmp$Site[i] &
    tmp$Season[tmp$last_sample_vec[i]] ==
    my_seas[which(my_seas == tmp$Season[i])-1]
  if(!my_eval){
    stop("last_sample_vec wrong, investigate.")
  }
}

# remove a bunch of stuff that is not needed
rm(list = c("new_tmp", "first_sites", "my_rows", "pb", "small_dat", "sp_recs",
     "sp_recs_list", "tdat", "tmp2", "to_go", "a1", "a2", "baddies",
     "my_eval", "my_seas"))

cat("\nReading in site covariates...\n")
# read in the within-city covariates
within_covs <- read.csv(
  "./data/cleaned_data/covariates/site_covariates.csv",
  stringsAsFactors = FALSE
)
  
# divide mean_19 by 100
within_covs$mean_19 <- within_covs$mean_19 / 100


# calculate the mean of each variable across cities
wc_mean <- within_covs %>% 
  dplyr::group_by(City) %>% 
  dplyr::summarise_if(is.numeric, mean, .groups = "drop_last")

# city-mean center the data
within_covs_centered <- split(
  within_covs,
  factor(within_covs$City)
)
for(i in 1:length(within_covs_centered)){
  within_covs_centered[[i]]$mean_19 <- 
    within_covs_centered[[i]]$mean_19 - mean(
      within_covs_centered[[i]]$mean_19
    )
}

within_covs_centered <- do.call("rbind", within_covs_centered)
# now we need to replicate this across our detection matrix.

tmp_covs <- tmp[,c("Site","City", "Season")]
tmp_covs$gentrifying <- FALSE
tmp_covs$mean_19 <- NA
for(i in 1:nrow(within_covs_centered)){
  my_loc <-    which( tmp_covs$Site == within_covs_centered$Site[i] &
    tmp_covs$City == within_covs_centered$City[i])
  if(length(my_loc) > 0){
    tmp_covs$gentrifying[my_loc] <- within_covs_centered$gentrifying[i]
    tmp_covs$mean_19[my_loc] <- within_covs_centered$mean_19[i]
  } else {
    stop()
  }
}

cat("Creating design matrices...\n")
# tmp covs is ordered like tmp, so we should be good to go with this
#  for psi, we are going to use all of the covariates
# -3 columns for site, season, and city. +1 column for intercept (i.e., -2)
psi_covs <- matrix(
  1,
  ncol = ncol(tmp_covs) - 2,
  nrow = nrow(tmp_covs)
)
psi_covs[,-1] <- as.matrix(
  tmp_covs[,-which(
    colnames(tmp_covs) %in% c("Site", "Season", "City")
  )]
)

# for detection (rho) we are going to use NDVI and population density
rho_covs <- matrix(
  1,
  ncol = ncol(tmp_covs) - 2,
  nrow = nrow(tmp_covs)
)
  
rho_covs[,-1] <- as.matrix(
  tmp_covs[,-which(
    colnames(tmp_covs) %in% c("Site", "Season", "City")
  )]
)
  
  
cat("Creating data_list and inits() function...\n")
data_list <- list(
  # detection data
  y = tmp$Y,
  J = tmp$J,
  # ids to fit the model in long format
  species_idx = tmp$Species_id,
  city_idx = tmp$City_id,
  combo_species_idx = combo$Species_id,
  combo_city_idx = combo$City_id,
  combo_idx = tmp$Combo_id,
  last_sample_vec = tmp$last_sample_vec,
  # covariates
  among_covs = among_covs,
  psi_covs = psi_covs,
  rho_covs = rho_covs,
  #season_covs = season_covs,
  # number of species, parameters, etc.
  nspecies = max(tmp$Species_id),
  ncity = max(tmp$City_id),
  ncov_among = dim(among_covs)[3],
  ncov_within = ncol(psi_covs),
  nsamples_one = sum(tmp$starter),
  nsamples_two = nrow(tmp),
  nseason_params = nrow(combo),
  #nseason_covs = ncol(season_covs),
  ncov_det = ncol(rho_covs)
)
  
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = rep(1, nrow(tmp)),
      x = matrix(1, ncol = data_list$ncity, nrow = data_list$nspecies),
      a_among = rnorm(data_list$ncov_among),
      tau_among = rgamma(data_list$ncov_among, 1, 1),
      b_among = matrix(rnorm(data_list$nspecies * data_list$ncov_among),
                       ncol = data_list$nspecies,
                       nrow = data_list$ncov_among),
      b_within = rnorm(data_list$ncov_within),
      tau_within = rgamma(data_list$ncov_within, 1, 1),
      b_species = matrix(
        rnorm(data_list$nspecies * data_list$ncov_within),
        ncol = data_list$nspecies,
        nrow = data_list$ncov_within),
      tau_species = matrix(
        rgamma(data_list$nspecies * data_list$ncov_within, 1, 1),
        ncol = data_list$nspecies,
        nrow = data_list$ncov_within),
      tau_shape = runif(data_list$ncov_within, 0.5,2),
      tau_rate = runif(data_list$ncov_within, 0.5,2),
      b_species_city = array(
        rnorm(data_list$ncov_within * data_list$nspecies * data_list$ncity),
        dim = c(data_list$ncov_within, data_list$nspecies, data_list$ncity)
      ),
      c_shape_psi = runif(1),
      c_rate_psi = runif(1),
      c_tau_psi = rgamma(data_list$ncity, 1, 1),
      ssc_psi = rnorm(data_list$nseason_params),
      c_shape_rho = runif(1),
      c_rate_rho = runif(1),
      c_tau_rho = rgamma(data_list$ncity, 1, 1),
      ssc_rho = rnorm(data_list$nseason_params),
      theta_mu = rnorm(1),
      tau_theta = rgamma(1,1,1),
      theta_shape = runif(1,0.5,1),
      theta_rate = runif(1, 0.5,1),
      theta_species = rnorm(data_list$nspecies),
      tau_theta_species = rgamma(data_list$nspecies, 1, 1),
      theta = matrix(
        rnorm(data_list$nspecies * data_list$ncity),
        nrow = data_list$nspecies,
        ncol = data_list$ncity
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

cat("prep_data_occupancy.R Complete!\n")
