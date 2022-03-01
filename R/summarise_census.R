############################
#
# Code to manipulate and summarise data pulled in via "./R/pull_census.R"
#
# Written by M. Fidino
#
############################

library(sf)
library(dplyr)
library(uwinspatialtools)

# read in some functions to help with this
source("./R/census_functions.R")

##
# the functions ----------------------------------------
##

# add_utms: used to add utms to the census data

# bbox_to_point: function used by add_utms()

# longlat_to_utm: self explanatory

# read_bind: read in census data and combine into list.

# interpolate_census_data: query data around census tracts, doing some
#  areal interpolation 

##
# the functions ----------------------------------------
##



# read in coordinates
coords <- read.csv(
  "./data/gentri_all_coordinates.csv"
)

# convert them to utms
coords$utms <- NA
for(i in 1:nrow(coords)){
coords$utms[i] <- longlat_to_utm(
  coords[i,c("Long","Lat")]
  )  
}

coords <- sf::st_as_sf(
  coords,
  coords = c("Long", "Lat"),
  crs = 4326
)

# split by utm crs
coords <- split(
  coords,
  factor(coords$utms)
)

# data types
my_folders <- c("education", "income", "race")

# poverty line
# read in income classifications
i_class <- read.csv(
  "./data/census_data/income_classifications.csv"
)


for(fldr in 2:length(my_folders)){

  cat(
    paste0(
      "\n census_data = ", my_folders[fldr],
      ". ", fldr, " of ", length(my_folders), ".\n"
    )
  )
  
  # read in race data
  census_data <- read_bind(
    list.files(
      paste0(
        "./data/census_data/",
        my_folders[fldr],
        "/"
      ),
      full.names = TRUE
    )
  )
  
  # number of years of census data
  nyear <- length(census_data)
  
  # add utms to each year
  for(year in 1:nyear){
    cat(paste0("\nyear ", year, " of ", nyear), "\n")
    # make column names lowercase
    colnames(census_data[[year]]) <- tolower(
      colnames(
        census_data[[year]]
      )
    )
    
    census_data[[year]] <- add_utms(
      census_data[[year]]
    )

    # convert variable names for income classes
    if(my_folders[fldr] == "income"){
      tmp_var <- census_data[[year]]$variable
      for(ii in 1:nrow(i_class)){
        tmp_var <- gsub(
          i_class$variable[ii],
          i_class[ii,year+1],
          tmp_var
        )
      }
      # and then sum across multiple records within each 'new' class
      if("estimate" %in% colnames(census_data[[year]])){
        census_data[[year]]$variable <- tmp_var
        tmp <- census_data[[year]]
        new_tmp <- tmp[,c("geoid", "name", "variable")]
        new_tmp <- new_tmp[!duplicated(new_tmp),]
        tmp <- data.frame(tmp)
        
        tmp <- tmp %>% 
          dplyr::group_by(geoid, variable) %>% 
          dplyr::summarise(
            estimate = sum(estimate, na.rm = TRUE),
            summary_est = unique(summary_est),
            crs = unique(crs),
            .groups = "drop_last"
          )
        new_tmp <- dplyr::inner_join(
          new_tmp,
          tmp,
          by = c("geoid", "variable")
        )
        census_data[[year]] <- new_tmp
        
      }else{
      census_data[[year]]$variable <- tmp_var
      tmp <- census_data[[year]]
      new_tmp <- tmp[,c("geoid", "name", "variable")]
      new_tmp <- new_tmp[!duplicated(new_tmp),]
      tmp <- data.frame(tmp)
      
      tmp <- tmp %>% 
        dplyr::group_by(geoid, variable) %>% 
        dplyr::summarise(
          value = sum(value, na.rm = TRUE),
          summary_value...total = unique(summary_value...total),
          crs = unique(crs),
          .groups = "drop_last"
        )
      new_tmp <- dplyr::inner_join(
        new_tmp,
        tmp,
        by = c("geoid", "variable")
      )
      census_data[[year]] <- new_tmp
      }
    }
  }
  
  # go through each year, split the df by utms, query data, and do areal
  #  interpolation.
  cat("\n splitting by crs.\n" )
  # split by utms
  for(year in 1:nyear){
    census_data[[year]] <- split(
      census_data[[year]],
      factor(census_data[[year]]$crs)
    )
  }
  
  # check if number of utms lines up with our sites
  if(!all(names(coords) %in% names(census_data$`2000`))){
    stop("utms are messed up, investigate.")
  }
  
  cat("\n Transforming projections.\n")
  # Convert each to their appropriate crs
  for(year in 1:nyear){
    for(crs in 1:length(census_data[[year]])){
      census_data[[year]][[crs]] <- sf::st_transform(
        census_data[[year]][[crs]],
        crs = as.numeric(names(census_data[[year]])[crs])
      )
    }
  }
  
  # same for sites
  for(crs in 1:length(coords)){
    coords[[crs]] <- sf::st_transform(
      coords[[crs]],
      crs = as.numeric(names(coords)[crs])
    )
  }
  
  cat("\n buffering sites.\n")
  my_buffs <- vector("list", length = length(coords))
  # make a buffer around the sites
  for(crs in 1:length(coords)){
    my_buffs[[crs]] <- sf::st_buffer(
      coords[[crs]], 
      dist = 1000
    )
  }
  names(my_buffs) <- names(coords)
  
  cat("\n querying data around sites.\n")
  # Abstract this code into a function so we can apply this to every 
  #  location and year. Just put the for loop in the function itself.
  
  to_save <- interpolate_census_data(
    census_data,
    my_buffs
  )
  write.csv(
    to_save,
    paste0(
      "./data/cleaned_data/", my_folders[fldr],".csv"
    ),
    row.names = FALSE
  )

}

