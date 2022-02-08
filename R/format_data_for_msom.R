# load functions used to clean data
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)

for(fn in functions_to_load){
  source(fn)
}
cli::cli_h1("Pulling in data")
# read in the data
dat <- read.csv(
  data_path
)
cli::cli_alert_success("Data loaded from {data_path}")


cli::cli_h1("Subsetting and formatting data")
# subset down to what we need
dat <- dat[
  which(dat$City %in% cities & dat$Species %in% species),
]

# subset down years
dat <- dat[grep(years, dat$Season),]

# make both phoenix study areas the same
dat$City <- gsub("phaz2", "phaz", dat$City)

# get just the unique sites and coordinates
all_sites <- dplyr::distinct(dat[,c("City", "Site","Long","Lat")])


# We need to figure out which sites to include, and we want to have
#  an autologistic term from one time period to the next. Because of
#  this I can't just remove all sites with no data.

# Drop all missing data
tmp <- dat[dat$J > 0,]

# order seasons appropriately
tmp$Season <- factor(tmp$Season, order_seasons(tmp$Season))

cli::cli_alert("Determining seasons to include for each city...")

# get number of sites per season
site_density <- tmp %>% 
  dplyr::group_by(City, Season) %>% 
  dplyr::summarise(
    nsite = length(unique(Site)),
    .groups = "drop_last"
  ) %>% 
  data.frame()

# remove seasons that have less than 15 sites per season
to_go <- site_density[site_density$nsite < 15,]

site_density <- site_density[-as.numeric(row.names(to_go)),]

tmp <- tmp[
  -which(
    paste(tmp$City, tmp$Season) %in%
    paste(to_go$City, to_go$Season)
  ),
]

tmp <- split(tmp, factor(tmp$City))
for(i in 1:length(tmp)){
  unq_city <- names(tmp)[i]
  unq_species <- unique(tmp[[i]]$Species)
  unq_season <- unique(tmp[[i]]$Season)
  unq_sites <- unique(tmp[[i]]$Site)
  all_cat_combo <- expand.grid(
    City = unq_city,
    Species = unq_species,
    Season = unq_season,
    Site = unq_sites
  )
  all_cat_combo <- dplyr::inner_join(
    all_cat_combo,
    all_sites,
    by = c("City", "Site")
  )
  tmp[[i]] <- dplyr::left_join(
    all_cat_combo,
    tmp[[i]][,-which(colnames(tmp[[i]]) %in% c("Long","Lat"))],
    by = c("City", "Species", "Season", "Site")
  )
  
}
tmp <- dplyr::bind_rows(tmp)
tmp$J[is.na(tmp$J)] <- 0

cli::cli_alert("Generating nested indexing variables...")
# This is a data.frame of the seasonal 'chunks' that each city has.
# If they have complete data then that city only has one row in the
# data.frame
season_check <- subset_seasons(
  tmp
)

# Create the species, city, and season id now
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
season_map <- season_map[order(season_map$Season_id),]

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

# add city_id to all_sites
all_sites <- dplyr::inner_join(
  all_sites,
  city_map,
  by = "City"
)

# make pretty names for cities (for plotting).
# Fill in when all cities have data
#pcity <- data.frame(
#  City = city_map$City,
#  Pretty = c(
#    "Atlanta,\nGeorgia",
#    "Austin,\nTexas",
#    "Chicago,\nIllinois",
#    "Denver,\nColorado",
#    "Edmonton,\nAlberta",
#    "Fort Collins,\nColorado",
#    "Iowa City,\nIowa",
#    "Indianapolis,\nIndiana",
#    "Jackson,\nMississippi",
#    "Manhattan,\nKansas",
#    "Madison,\nWisconsin",
#    "Metro LA,\nCalifornia",
#    "Phoenix,\nArizona",
#    "Rochester,\nNew York",
#    "Sanford,\nFlorida",
#    "Salt Lake\nCity, Utah",
#    "Seattle,\nWashington",
#    "Saint Louis,\nMissouri",
#    "Tacoma,\nWashington",
#    "Wilmington,\nDelaware"
#  )
#)

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
  ) %>% 
  data.frame()

# this is all the unique combos of species, seasons, and cities. Can
#  be used to index the random effect we use.
combo <- tmp[,
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
tmp <- dplyr::left_join(
  tmp,
  combo,
  by = c("Species_id", "Season_id", "City_id")
)

# Split data into first season vs the rest, this is for
#  adding the autologistic term into the model
first_season <- vector("list", length = nrow(season_check))

the_rest <- tmp

# grab the "firsts" and then leave "the rest"
for(i in seq_len(nrow(season_check))){
  first_season[[i]] <- tmp[
    tmp$City == season_check$City[i] &
    tmp$Season == season_check$Season[i],
  ]
  # remove first season from the rest
  the_rest <- the_rest[-which(
    the_rest$City == season_check$City[i] & 
    the_rest$Season == season_check$Season[i]
  ),]
}


first_season <- dplyr::bind_rows(first_season)

# Will need to know this to get the model started
nfirsts <- nrow(first_season)

# recombine first_season and the_rest now so first_season is on top
tmp_reor <- dplyr::bind_rows(
  list(first_season, the_rest)
)


# This will hold the vector for prior sampling
last_sample <- rep(NA, nrow(tmp_reor))

for(i in 1:nrow(city_map)){
  the_city <- city_map$City[i]
  city_seasons <- unique(tmp$Season[tmp$City == the_city])
  # go to next if there is only one season
  if(length(city_seasons) == 1){
    next
  }
  for(j in 1:length(city_seasons)){
    # if the season is in season_check then we just go
    #  to the next. There is no autologistic term there.
    if(
      paste(the_city, city_seasons[j]) %in%
      paste(season_check$City, season_check$Season)
    ){
      next
    }
    # Otherwise compare time t to t-1 and get the index
      tmp_prior_season <- tmp_reor[
        tmp_reor$City == the_city &
        tmp_reor$Season == city_seasons[j-1],
      ]
      tmp_prior_season$row_vec <- as.numeric(row.names(tmp_prior_season))
      # join to just that specific season
      first_join <- left_join(
        tmp_reor[
          tmp_reor$City == the_city &
            tmp_reor$Season == city_seasons[j],
          
          ],
        tmp_prior_season[,c("Species", "Site", "row_vec")],
        by = c("Species", "Site")
      ) %>% 
        left_join(
          tmp_reor,
          .,
          by = c("City", "Species", "Season","Site")
        )
      if(any(!is.na(last_sample[!is.na(first_join$row_vec)]))){
        stop("you screwed up")
      } else {
        last_sample[
          !is.na(first_join$row_vec)
        ] <- first_join$row_vec[!is.na(first_join$row_vec)]
      }
    }
}

# check a random subset
if(test_autologistic){
  cli::cli_alert("Testing autologistic nested index term...")
  to_check <- last_sample[sample((nfirsts+1):length(last_sample),30000)]
  pb <- txtProgressBar(max = length(to_check))
  for(i in 1:length(to_check)){
    setTxtProgressBar(pb, i)
    # one sample
    os <- tmp_reor[to_check[i],]
    expected_loc <- which(
      tmp_reor$City == os$City &
      tmp_reor$Species == os$Species &
      tmp_reor$Season == season_map$Season[
        which(season_map$Season == os$Season) + 1] &
      tmp_reor$Site == os$Site
    )
    # one sample test
    ost <- tmp_reor[expected_loc,]
    if(
      ost$City == os$City &
      ost$Species == os$Species &
      season_map$Season[
        which(season_map$Season == ost$Season)-1
      ] == os$Season &
      ost$Site == os$Site
    ){
      next
    } else{
      cli::cli_abort("There is a mismatch in the autologistic term creation.")
    }
  }
cli::cli_alert_success("Autologistic indexing term in working order")
rm(os, ost, to_check, expected_loc, pb)
}

# add this index term to the data.
tmp_reor$last_samplevec <- last_sample

# remove a bunch of temporary objects
rm(
  tmp, to_go, tmp_prior_season,
  all_cat_combo, dat, first_join, first_season, season_check,
  site_density, the_rest, last_sample, city_seasons, fn,
  functions_to_load, i, j, the_city, unq_city, unq_season,
  unq_sites, unq_species
)

cli::cli_alert_success("Data subsetted and formatted")

cli::cli_h1("Creating data_list for JAGS")
# make the design matrices for the model
omegacov <- array(1, dim = c(nrow(species_map), nrow(city_map), 2))
omegacov[,,2] <- rnorm(prod(dim(omegacov)[1:2]))
psicov <- matrix(1, ncol = 2, nrow = nrow(tmp_reor))
psicov[,2] <- rnorm(nrow(tmp_reor))
rhocov <- matrix(1, ncol = 2, nrow = nrow(tmp_reor))
rhocov[,2] <- rnorm(nrow(tmp_reor))

data_list <- list(
  # detection / non-detection data
  y = tmp_reor$Y,
  J = tmp_reor$J,
  # design matrices
  omegacov = omegacov,
  psicov = psicov,
  rhocov = rhocov
)
constant_list <- list(
  # nested indexing
  species_vec = tmp_reor$Species_id,
  season_vec = tmp_reor$Season_id,
  city_vec = tmp_reor$City_id,
  combo_species_vec = combo$Species_id,
  combo_city_vec = combo$City_id,
  combo_vec = tmp_reor$Combo_id,
  last_sample_vec = tmp_reor$last_samplevec,
  #constants,
  nfirsts = nfirsts,
  ndata = nrow(tmp_reor),
  nspecies = nrow(species_map),
  ncity = nrow(city_map),
  nomega = dim(omegacov)[3],
  npsi = ncol(psicov),
  nrho = ncol(rhocov),
  nseason_params = nrow(combo)
)

jags_list <- list(
  # detection / non-detection data
  y = tmp_reor$Y,
  J = tmp_reor$J,
  # design matrices
  omegacov = omegacov,
  psicov = psicov,
  rhocov = rhocov,
  # nested indexing
  species_vec = tmp_reor$Species_id,
  season_vec = tmp_reor$Season_id,
  city_vec = tmp_reor$City_id,
  combo_species_vec = combo$Species_id,
  combo_city_vec = combo$City_id,
  combo_vec = tmp_reor$Combo_id,
  last_sample_vec = tmp_reor$last_samplevec,
  #constants,
  nfirsts = nfirsts,
  ndata = nrow(tmp_reor),
  nspecies = nrow(species_map),
  ncity = nrow(city_map),
  nomega = dim(omegacov)[3],
  npsi = ncol(psicov),
  nrho = ncol(rhocov),
  nseason_params = nrow(combo)
)

cli::cli_alert_success("object 'jags_list' created and ready for modeling")

if(analysis == "msom"){
  rm(psicov, rhocov, omegacov, nfirsts, tmp_reor, combo)
} else {
  sp_dat <- tmp_reor
  
  combo$Season_city_id <- paste0(combo$Season_id, "-", combo$City_id)
  combo$Season_city_id <- as.numeric(
    factor(
      combo$Season_city_id,
      levels = unique(combo$Season_city_id)
    )
  )
  sp_dat <- dplyr::inner_join(
    sp_dat,
    combo[,c("Combo_id", "Season_city_id")],
    by = "Combo_id"
  )
  
  rm(psicov, rhocov, omegacov, nfirsts, tmp_reor)
  
}


