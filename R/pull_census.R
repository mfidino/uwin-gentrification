# pull census data at the tract level

# this code could have been abstracted out A LOT to reduce
#  it's length. But since this was a 'one and done' analysis
#  I've kept it as is.


library(dplyr)
library(tidycensus)
library(sf)

open_urls <- FALSE

# Note: you need to get an api from the census website 

# input api key
tidycensus::census_api_key(
  key = "{PLACE YOUR API KEY HERE}",
  overwrite = TRUE,
  install = TRUE
)

# The first that needs to be done is pulling in all US counties.
#  This is because we have to determine which counties to query in
#  the census data. 

if(open_urls){
  browseURL(
    "https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html"
  )
}


# make a list that stores the correct tables
#  to query for each dataset.

# load in a table of variables

all_years <- c(
  2000, 2010, 2015, 2019
)

if(
  file.exists("./data/census_data/var_tables.RDS")
){
  my_vars <- readRDS(
    "./data/census_data/var_tables.RDS"
  )
} else {

  my_vars <- list(
    yr00sf1 = tidycensus::load_variables(
      2000,
      "sf1"
    ),
    yr00sf3 = tidycensus::load_variables(
      2000,
      "sf3"
    ),
    yr10sf1 = tidycensus::load_variables(
      2010,
      "sf1"
    ),
    yr10acs5 = tidycensus::load_variables(
      2010,
      "acs5"
    ),
    yr15 = tidycensus::load_variables(
      2015,
      "acs5"
    ),
    yr19 = tidycensus::load_variables(
      2019,
      "acs5"
    )
  )
  
  saveRDS(
    my_vars,
    "./data/census_data/var_tables.RDS"
  )

}

# start figuring out the tables we need
# this has leading zeroes, so we need to figure
#  out the exact "name" for each table
#  we need tables P1 and P5.

race <- vector(
  "list",
  length = length(all_years)
)


names(race) <- all_years

race$`2000` <- list(
  year = 2000,
  census = "decennial"
)

# start looking for elements to add to this

# looks like it is just 1 record for P001
my_vars$yr00sf1[grep("P001", my_vars$yr00sf1$name),]


# now P008,
my_vars$yr00sf1[grep("P008", my_vars$yr00sf1$name),]

race$`2000`$vars <-  c(
  "Total" = "P001001",
  "White" = "P008003",
  "Black" = "P008004",
  "Asian" = "P008006",
  "Latino" = "P008010"
)


# 2010 is decennial, so different.


# looks like it is just 1 record for P001
my_vars$yr10sf1[grep("P001", my_vars$yr10sf1$name),]
my_vars$yr10sf1[grep("P005", my_vars$yr10sf1$name),]

race$`2010` <- list(
  year = 2010,
  census = "decennial",
  vars = c(
    "Total" = "P001001",
    "White" = "P005003",
    "Black" = "P005004",
    "Asian" = "P005006",
    "Latino" = "P005010"
  )
)


# check 2015
my_vars$yr15[grep("!!Hispanic or Latino", my_vars$yr15$label),1:2]
race$`2015` <- list(
  year = 2015,
  census = "acs",
  vars = c(
    "Total" = "B01001_001",
    "White" = "B03002_003",
    "Black" = "B03002_004",
    "Asian" = "B03002_006",
    "Latino" = "B03001_003"
  )
)


# check 2019, same again.
my_vars$yr19[grep("!!Not Hispanic or Latino", my_vars$yr19$label),1:2]
race$`2019` <- race$`2015`
race$`2019`$year <- 2019

# get race data from all of these years

for(year in 1:length(race)){

# read in the county data
  
  
  counties <- sf::read_sf(
    "D:/GIS/counties/cb_2020_us_county_500k"
  )
  # read in spatial data, and reproject to counties crs
  coords <- read.csv(
    "./data/gentri_all_coordinates.csv"
  )
  
  coords <- sf::st_as_sf(
    coords,
    coords = c("Long", "Lat"),
    crs = 4326
  )

  coords <- sf::st_transform(
    coords,
    sf::st_crs(counties)
  )

  # buffer coordinates, just by some
  #  small value (this is > 1000m buffer)
  coords_buffer <- sf::st_buffer(
    coords,
    dist = 0.04
  )

  # get just the counties that intersect these coordinates
  my_counties_idx <- unique(
    sort(
      unlist(
        sf::st_intersects(
          coords_buffer,
          counties,
        )
      )
    )
  )
  # subset down to just the counties we need
  counties <- counties[my_counties_idx,]

  counties <- data.frame(
    counties[,c("STATEFP", "COUNTYFP", "NAME", "STUSPS")]
  )[,1:4]


counties <- split(
  counties,
  factor(counties$STUSPS)
)

my_results <- vector(
  "list",
  length = length(counties)
)

for(i in 1:length(my_results)){
  # split counties by state
  longshot <- TRUE
  mc <- 1
  while(longshot){
    
    if(race[[year]]$census == "decennial"){
      tmp <- try(
        get_decennial(
          geography = "tract",
          variables = race[[year]]$vars[-1],
          year = race[[year]]$year,
          summary_var = race[[year]]$vars[1],
          state = unique(counties[[i]]$STATEFP),
          county = counties[[i]]$COUNTYFP,
          geometry = TRUE,
          show_call = TRUE
        ),
        silent = TRUE
      )
    } else {
      tmp <- try(
        get_acs(
          geography = "tract",
          variables = race[[year]]$vars[-1],
          year = race[[year]]$year,
          summary_var = race[[year]]$vars[1],
          state = unique(counties[[i]]$STATEFP),
          county = counties[[i]]$COUNTYFP,
          geometry = TRUE,
          show_call = TRUE
        ),
        silent = TRUE
      )
    }
    if(any(class(tmp) == "sf")){
      longshot <- FALSE
    } else {
      mc <- mc+1
    }
    if(mc>10){
      stop("badness.")
    }
  }
  
  my_results[[i]] <- tmp
  
}

saveRDS(
  my_results,
  paste0(
    "./data/census_data/race/race_",
    race[[year]]$year,
    ".RDS"
  )
)

}

# Do the same with income totals data

# just going to clear up the names and whatnot
#  from the variables


income <- vector(
  "list",
  length = length(all_years)
)

names(income) <- all_years

# this should clean all the income stuff
inc_clean <- function(x){
  x <- gsub("Estimate!!","", x)
  x <- gsub("!!", "_", x)
  x <- gsub(",000","", x)
  x <- gsub(",999", "", x)
  x <- gsub("\\$", "", x)
  x <- gsub(" to ", "_to_", x)
  x <- gsub("Less than ", "under_", x)
  x <- gsub("Total_", "btw_", x)
  x <- gsub("btw_200 or more", "over_200", x)
  x <- gsub("btw_under", "under", x)
  
  x <- tolower(x)
  return(x)
}

# 2000
#View(my_vars$yr00sf3)
tmp <- my_vars$yr00sf3[grep("P076",  my_vars$yr00sf3$name),]

income$`2000` <- list(
  year = 2000,
  census = "decennial",
  vars = tmp$name
)

names(income$`2000`$vars) <- inc_clean(
  tmp$label
)

# 2010
tmp <- my_vars$yr10acs5[grep("B19001", my_vars$yr10acs5$name),][1:17,]

income$`2010` <- list(
  year = 2010,
  census = "acs",
  vars = tmp$name
)

names(income$`2010`$vars) <- inc_clean(
  tmp$label
)

# 2015, same as 2010
my_vars$yr15[grep("B19001", my_vars$yr15$name),1:2]
income$`2015` <- income$`2010`
income$`2015`$year <- 2015

# 2019, same as 2010
my_vars$yr19[grep("B19001", my_vars$yr19$name),1:2]
income$`2019` <- income$`2010`
income$`2019`$year <- 2019



for(year in 1:length(income)){
  
  # read in the county data
  
  
  counties <- sf::read_sf(
    "D:/GIS/counties/cb_2020_us_county_500k"
  )
  # read in spatial data, and reproject to counties crs
  coords <- read.csv(
    "./data/gentri_all_coordinates.csv"
  )
  
  coords <- sf::st_as_sf(
    coords,
    coords = c("Long", "Lat"),
    crs = 4326
  )
  
  coords <- sf::st_transform(
    coords,
    sf::st_crs(counties)
  )
  
  # buffer coordinates, just by some
  #  small value (this is > 1000m buffer)
  coords_buffer <- sf::st_buffer(
    coords,
    dist = 0.04
  )
  
  # get just the counties that intersect these coordinates
  my_counties_idx <- unique(
    sort(
      unlist(
        sf::st_intersects(
          coords_buffer,
          counties,
        )
      )
    )
  )
  # subset down to just the counties we need
  counties <- counties[my_counties_idx,]
  
  counties <- data.frame(
    counties[,c("STATEFP", "COUNTYFP", "NAME", "STUSPS")]
  )[,1:4]
  
  
  counties <- split(
    counties,
    factor(counties$STUSPS)
  )
  
  my_results <- vector(
    "list",
    length = length(counties)
  )
  
  for(i in 1:length(my_results)){
    # split counties by state
    longshot <- TRUE
    mc <- 1
    while(longshot){
      
      if(income[[year]]$census == "decennial"){
        tmp <- try(
          get_decennial(
            geography = "tract",
            variables = income[[year]]$vars[-1],
            year = income[[year]]$year,
            summary_var = income[[year]]$vars[1],
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE,
            sumfile = "sf3"
          ),
          silent = TRUE
        )
      } else {
        tmp <- try(
          get_acs(
            geography = "tract",
            variables = income[[year]]$vars[-1],
            year = income[[year]]$year,
            summary_var = income[[year]]$vars[1],
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE
          ),
          silent = TRUE
        )
      }
      if(any(class(tmp) == "sf")){
        longshot <- FALSE
      } else {
        mc <- mc+1
      }
      if(mc>10){
        stop("badness.")
      }
    }
    
    my_results[[i]] <- tmp
    
  }
  
  saveRDS(
    my_results,
    paste0(
      "./data/census_data/income/income_",
      income[[year]]$year,
      ".RDS"
    )
  )
  
}


## Education

edu <- vector(
  "list",
  length(all_years)
)

names(edu) <- all_years

# this should clean all the edu stuff as needed
#  that being said, we'll need to apply this cleaning
#  after we query the data.


edu_clean <- function(x){
  x <- gsub("Estimate!!","", x)
  x <- gsub("!!", "_", x)
  x <- gsub("_Male|_Female","", x)
  # drop the 2nd and 3rd totals
  x
  x <- gsub(
    "Total_No schooling completed|.*grade.*",
    "No HS diploma",
    x
  )
  x <- gsub(".*graduate.*", "HS diploma", x)
  x <- gsub(".*Some college.*", "Some college", x)
  x <- gsub(".*Associate.*", "College graduate", x)
  x <- gsub(".*Bachelor.*", "College graduate", x)
  x <- gsub(".*degree", "Advanced college degree",x)
  x <- tolower(x)
  return(x)
}

# 2000
tmp <- my_vars$yr00sf3[grep("P037", my_vars$yr00sf3$name), 1:2]
# drop the sub-totals that we do not need
tmp <- tmp[-which(tmp$label %in% c("Total!!Male", "Total!!Female")),]


edu$`2000` <- list(
  year = 2000,
  census = "decennial",
  vars = tmp$name
)

names(edu$`2000`$vars) <- tmp$label

# 2010
my_vars$yr10acs5[grep("Doctorate", my_vars$yr10acs5$label),]
tmp <- my_vars$yr10acs5[grep("B15002", my_vars$yr10acs5$name),]
tmp$label <- gsub("Estimate!!", "", tmp$label)
tmp <- tmp[-which(tmp$label %in% c("Total!!Male", "Total!!Female")),]

edu$`2010` <- list(
  year = 2010,
  census = "acs",
  vars = tmp$name
)

names(edu$`2010`$vars) <- tmp$label


# 2015, same as 2010
my_vars$yr15[grep("Doctorate", my_vars$yr15$label),]
edu$`2015` <- edu$`2010`
edu$`2015`$year <- 2015

# 2019, same as 2010
my_vars$yr19[grep("Doctorate", my_vars$yr19$label),]
edu$`2019` <- edu$`2010`
edu$`2019`$year <- 2019



for(year in 1:length(edu)){
  
  # read in the county data
  
  
  counties <- sf::read_sf(
    "D:/GIS/counties/cb_2020_us_county_500k"
  )
  # read in spatial data, and reproject to counties crs
  coords <- read.csv(
    "./data/gentri_all_coordinates.csv"
  )
  
  coords <- sf::st_as_sf(
    coords,
    coords = c("Long", "Lat"),
    crs = 4326
  )
  
  coords <- sf::st_transform(
    coords,
    sf::st_crs(counties)
  )
  
  # buffer coordinates, just by some
  #  small value (this is > 1000m buffer)
  coords_buffer <- sf::st_buffer(
    coords,
    dist = 0.04
  )
  
  # get just the counties that intersect these coordinates
  my_counties_idx <- unique(
    sort(
      unlist(
        sf::st_intersects(
          coords_buffer,
          counties,
        )
      )
    )
  )
  # subset down to just the counties we need
  counties <- counties[my_counties_idx,]
  
  counties <- data.frame(
    counties[,c("STATEFP", "COUNTYFP", "NAME", "STUSPS")]
  )[,1:4]
  
  
  counties <- split(
    counties,
    factor(counties$STUSPS)
  )
  
  my_results <- vector(
    "list",
    length = length(counties)
  )
  
  for(i in 1:length(my_results)){
    print(i)
    # split counties by state
    longshot <- TRUE
    mc <- 1
    while(longshot){
      print(paste0("mc = ", mc))
      
      if(edu[[year]]$census == "decennial"){
        tmp <- try(
          get_decennial(
            geography = "tract",
            variables = edu[[year]]$vars[-1],
            year = edu[[year]]$year,
            summary_var = edu[[year]]$vars[1],
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE,
            sumfile = "sf3"
          ),
          silent = TRUE
        )
        if( any(class(tmp) == "sf")){
          # clean up variable names
          tmp$variable <- edu_clean(tmp$variable)
          # and summarise based on the groups we created
          tmp <- tmp %>% dplyr::group_by(
            GEOID, variable, summary_value...Total
          ) %>%
            dplyr::summarise(
              value = sum(value, na.rm = TRUE),
              .groups = "drop_last"
            )
        }
      } else {
        tmp <- try(
          get_acs(
            geography = "tract",
            variables = edu[[year]]$vars[-1],
            year = edu[[year]]$year,
            summary_var = edu[[year]]$vars[1],
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE
          ),
          silent = TRUE
        )
        # clean up variable names
        if( any(class(tmp) == "sf")){
          tmp$variable <- edu_clean(tmp$variable)
          # and summarise based on the groups we created
          tmp <- tmp %>% dplyr::group_by(
            GEOID, variable, summary_est
          ) %>%
            dplyr::summarise(
              value = sum(estimate, na.rm = TRUE),
              .groups = "drop_last"
            )
          }
      }
      if(any(class(tmp) == "sf")){
        longshot <- FALSE
      } else {
        mc <- mc+1
      }
      if(mc>100){
        stop("badness.")
      }
    }
    
    my_results[[i]] <- tmp
    
  }
  
  saveRDS(
    my_results,
    paste0(
      "./data/census_data/education/education_",
      edu[[year]]$year,
      ".RDS"
    )
  )
}


## Housing

hou <- vector(
  "list",
  length(all_years)
)

names(hou) <- all_years

# 2000
tmp <- my_vars$yr00sf3[grep("H001001", my_vars$yr00sf3$name), 1:2]
# drop the sub-totals that we do not need


hou$`2000` <- list(
  year = 2000,
  census = "decennial",
  vars = tmp$name
)

names(hou$`2000`$vars) <- tmp$label

# 2010

#View(my_vars$yr10sf1)

tmp <- my_vars$yr10acs5[grep("B25001_001", my_vars$yr10acs5$name),]
tmp$label <- gsub("Estimate!!", "", tmp$label)


hou$`2010` <- list(
  year = 2010,
  census = "acs",
  vars = tmp$name
)

names(hou$`2010`$vars) <- tmp$label


# 2015, same as 2010
my_vars$yr15[grep("B25001_001", my_vars$yr15$name),]
hou$`2015` <- hou$`2010`
hou$`2015`$year <- 2015

# 2019, same as 2015
my_vars$yr19[grep("B25001_001", my_vars$yr19$name),]
hou$`2019` <- hou$`2010`
hou$`2019`$year <- 2019



for(year in 1:length(hou)){
  
  # read in the county data
  
  
  counties <- sf::read_sf(
    "D:/GIS/counties/cb_2020_us_county_500k"
  )
  # read in spatial data, and reproject to counties crs
  coords <- read.csv(
    "./data/gentri_all_coordinates.csv"
  )
  
  coords <- sf::st_as_sf(
    coords,
    coords = c("Long", "Lat"),
    crs = 4326
  )
  
  coords <- sf::st_transform(
    coords,
    sf::st_crs(counties)
  )
  
  # buffer coordinates, just by some
  #  small value (this is > 1000m buffer)
  coords_buffer <- sf::st_buffer(
    coords,
    dist = 0.04
  )
  
  # get just the counties that intersect these coordinates
  my_counties_idx <- unique(
    sort(
      unlist(
        sf::st_intersects(
          coords_buffer,
          counties,
        )
      )
    )
  )
  # subset down to just the counties we need
  counties <- counties[my_counties_idx,]
  
  counties <- data.frame(
    counties[,c("STATEFP", "COUNTYFP", "NAME", "STUSPS")]
  )[,1:4]
  
  
  counties <- split(
    counties,
    factor(counties$STUSPS)
  )
  
  my_results <- vector(
    "list",
    length = length(counties)
  )
  
  for(i in 1:length(my_results)){
    print(i)
    # split counties by state
    longshot <- TRUE
    mc <- 1
    while(longshot){
      print(paste0("mc = ", mc))
      
      if(hou[[year]]$census == "decennial"){
        tmp <- try(
          get_decennial(
            geography = "tract",
            variables = hou[[year]]$vars,
            year = hou[[year]]$year,
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE,
            sumfile = "sf3"
          ),
          silent = TRUE
        )
        if( any(class(tmp) == "sf")){
          # clean up variable names
          # and summarise based on the groups we created
          tmp <- tmp %>% dplyr::group_by(
            GEOID, variable
          ) %>%
            dplyr::summarise(
              value = sum(value, na.rm = TRUE),
              .groups = "drop_last"
            )
        }
      } else {
        tmp <- try(
          get_acs(
            geography = "tract",
            variables = hou[[year]]$vars,
            year = hou[[year]]$year,
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE
          ),
          silent = TRUE
        )
        # clean up variable names
        if( any(class(tmp) == "sf")){
          # and summarise based on the groups we created
          tmp <- tmp %>% dplyr::group_by(
            GEOID, variable
          ) %>%
            dplyr::summarise(
              value = sum(estimate, na.rm = TRUE),
              .groups = "drop_last"
            )
        }
      }
      if(any(class(tmp) == "sf")){
        longshot <- FALSE
      } else {
        mc <- mc+1
      }
      if(mc>100){
        stop("badness.")
      }
    }
    
    my_results[[i]] <- tmp
    
  }
  
  saveRDS(
    my_results,
    paste0(
      "./data/census_data/housing/housing_",
      hou[[year]]$year,
      ".RDS"
    )
  )
}



## Housing prices

pri <- vector(
  "list",
  length(all_years)
)

names(pri) <- all_years

# 2000
#View(my_vars$yr00sf3)
tmp <- my_vars$yr00sf3[grep("H076001", my_vars$yr00sf3$name), 1:2]
# drop the sub-totals that we do not need


pri$`2000` <- list(
  year = 2000,
  census = "decennial",
  vars = tmp$name
)

names(pri$`2000`$vars) <- tmp$label

# 2010

#View(my_vars$yr10acs5)

tmp <- my_vars$yr10acs5[grep("B25077_001", my_vars$yr10acs5$name),]
tmp$label <- gsub("Estimate!!", "", tmp$label)
tmp$label <- gsub("\\s\\(dollars\\)", "", tmp$label)

pri$`2010` <- list(
  year = 2010,
  census = "acs",
  vars = tmp$name
)

names(pri$`2010`$vars) <- tmp$label


# 2015, same as 2010
my_vars$yr15[grep("B25077_001", my_vars$yr15$name),]
pri$`2015` <- pri$`2010`
pri$`2015`$year <- 2015

# 2019, same as 2015
my_vars$yr19[grep("B25077_001", my_vars$yr19$name),]
pri$`2019` <- pri$`2010`
pri$`2019`$year <- 2019



for(year in 1:length(pri)){
  
  # read in the county data
  
  
  counties <- sf::read_sf(
    "D:/GIS/counties/cb_2020_us_county_500k"
  )
  # read in spatial data, and reproject to counties crs
  coords <- read.csv(
    "./data/gentri_all_coordinates.csv"
  )
  
  coords <- sf::st_as_sf(
    coords,
    coords = c("Long", "Lat"),
    crs = 4326
  )
  
  coords <- sf::st_transform(
    coords,
    sf::st_crs(counties)
  )
  
  # buffer coordinates, just by some
  #  small value (this is > 1000m buffer)
  coords_buffer <- sf::st_buffer(
    coords,
    dist = 0.04
  )
  
  # get just the counties that intersect these coordinates
  my_counties_idx <- unique(
    sort(
      unlist(
        sf::st_intersects(
          coords_buffer,
          counties,
        )
      )
    )
  )
  # subset down to just the counties we need
  counties <- counties[my_counties_idx,]
  
  counties <- data.frame(
    counties[,c("STATEFP", "COUNTYFP", "NAME", "STUSPS")]
  )[,1:4]
  
  
  counties <- split(
    counties,
    factor(counties$STUSPS)
  )
  
  my_results <- vector(
    "list",
    length = length(counties)
  )
  
  for(i in 1:length(my_results)){
    print(i)
    # split counties by state
    longshot <- TRUE
    mc <- 1
    while(longshot){
      print(paste0("mc = ", mc))
      
      if(pri[[year]]$census == "decennial"){
        tmp <- try(
          get_decennial(
            geography = "tract",
            variables = pri[[year]]$vars,
            year = pri[[year]]$year,
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE,
            sumfile = "sf3"
          ),
          silent = TRUE
        )
        if( any(class(tmp) == "sf")){
          # clean up variable names
          # and summarise based on the groups we created
          tmp <- tmp %>% dplyr::group_by(
            GEOID, variable
          ) %>%
            dplyr::summarise(
              value = sum(value, na.rm = TRUE),
              .groups = "drop_last"
            )
        }
      } else {
        tmp <- try(
          get_acs(
            geography = "tract",
            variables = pri[[year]]$vars,
            year = pri[[year]]$year,
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE
          ),
          silent = TRUE
        )
        # clean up variable names
        if( any(class(tmp) == "sf")){
          # and summarise based on the groups we created
          tmp <- tmp %>% dplyr::group_by(
            GEOID, variable
          ) %>%
            dplyr::summarise(
              value = sum(estimate, na.rm = TRUE),
              .groups = "drop_last"
            )
        }
      }
      if(any(class(tmp) == "sf")){
        longshot <- FALSE
      } else {
        mc <- mc+1
      }
      if(mc>100){
        stop("badness.")
      }
    }
    
    my_results[[i]] <- tmp
    
  }
  
  saveRDS(
    my_results,
    paste0(
      "./data/census_data/housing_price/housing_price_",
      pri[[year]]$year,
      ".RDS"
    )
  )
}


## median income


# Do the same with income totals data

# just going to clear up the names and whatnot
#  from the variables


med_income <- vector(
  "list",
  length = length(all_years)
)

names(med_income) <- all_years


# 2000

tmp <- my_vars$yr00sf3[grep("HCT012001",  my_vars$yr00sf3$name),]

med_income$`2000` <- list(
  year = 2000,
  census = "decennial",
  vars = tmp$name
)


# 2010

tmp <- my_vars$yr10acs5[grep("B19113_001", my_vars$yr10acs5$name),]

med_income$`2010` <- list(
  year = 2010,
  census = "acs",
  vars = tmp$name
)


# 2015, same as 2010
my_vars$yr15[grep("B19113_001", my_vars$yr15$name),1:2]
med_income$`2015` <- med_income$`2010`
med_income$`2015`$year <- 2015

# 2019, same as 2010
my_vars$yr19[grep("B19113_001", my_vars$yr19$name),1:2]
med_income$`2019` <- med_income$`2010`
med_income$`2019`$year <- 2019



for(year in 2:length(med_income)){
  
  # read in the county data
  
  
  counties <- sf::read_sf(
    "D:/GIS/counties/cb_2020_us_county_500k"
  )
  # read in spatial data, and reproject to counties crs
  coords <- read.csv(
    "./data/gentri_all_coordinates.csv"
  )
  
  coords <- sf::st_as_sf(
    coords,
    coords = c("Long", "Lat"),
    crs = 4326
  )
  
  coords <- sf::st_transform(
    coords,
    sf::st_crs(counties)
  )
  
  # buffer coordinates, just by some
  #  small value (this is > 1000m buffer)
  coords_buffer <- sf::st_buffer(
    coords,
    dist = 0.04
  )
  
  # get just the counties that intersect these coordinates
  my_counties_idx <- unique(
    sort(
      unlist(
        sf::st_intersects(
          coords_buffer,
          counties,
        )
      )
    )
  )
  # subset down to just the counties we need
  counties <- counties[my_counties_idx,]
  
  counties <- data.frame(
    counties[,c("STATEFP", "COUNTYFP", "NAME", "STUSPS")]
  )[,1:4]
  
  
  counties <- split(
    counties,
    factor(counties$STUSPS)
  )
  
  my_results <- vector(
    "list",
    length = length(counties)
  )
  
  for(i in 1:length(my_results)){
    # split counties by state
    longshot <- TRUE
    mc <- 1
    while(longshot){
      
      if(med_income[[year]]$census == "decennial"){
        tmp <- try(
          get_decennial(
            geography = "tract",
            variables = med_income[[year]]$vars,
            year = med_income[[year]]$year,
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE,
            sumfile = "sf3"
          ),
          silent = TRUE
        )
      } else {
        tmp <- try(
          get_acs(
            geography = "tract",
            variables = med_income[[year]]$vars,
            year = med_income[[year]]$year,
            state = unique(counties[[i]]$STATEFP),
            county = counties[[i]]$COUNTYFP,
            geometry = TRUE,
            show_call = TRUE
          ),
          silent = TRUE
        )
      }
      if(any(class(tmp) == "sf")){
        longshot <- FALSE
      } else {
        mc <- mc+1
      }
      if(mc>10){
        stop("badness.")
      }
    }
    
    my_results[[i]] <- tmp
    
  }
  
  saveRDS(
    my_results,
    paste0(
      "./data/census_data/med_income/med_income_",
      med_income[[year]]$year,
      ".RDS"
    )
  )
  
}

