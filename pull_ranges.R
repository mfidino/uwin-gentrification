###################################
#
# Pull range info
# Written by M. Fidino
#
#
####################################

### load packages
library(sf)
library(dplyr)
library(mapview)
library(smoothr)

# load functions for this project
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)

for(fn in functions_to_load){
  source(fn)
}

# I am using the IUCN Spatial data to determine the range of terrestrial
#  mammals throughout north america. On 6/9/2020 these data could be
#  queried at (uncomment and run next line to open up link):
#
#browseURL("https://www.iucnredlist.org/resources/spatial-data-download")

# I am using the world continents data from the HIFLD site. On 6/9/2020
#  these could be found at (uncomment and run next line to open up link):
#
#browseURL("https://hifld-geoplatform.opendata.arcgis.com/datasets/bee7adfd918e4393995f64e155a1bbdf_0")
#
#browseURL(http://gis-michigan.opendata.arcgis.com/datasets/6031c4fb8cac48649f2e0a98999d1248_0)
#
# I am using the great lakes data from the state of Michigan. On 6/11/2020
#  these could be found at (uncomment and run next line to open up link):


# After downloading and unzipping I am assuming that the folders of spatial
#  data is here (update as necessary):
mammal_path <- "D:/TERRESTRIAL_MAMMALS"
cont_path <- "D:/GIS/continents"
lake_path <- "D:/GIS/great_lakes"

if(!file.exists(mammal_path)){
  stop("You did not specify where the IUCN spatial data is in the script 'pull_range_info.R'")
}

if(!file.exists(cont_path)){
  stop("You did not specify where the continent shapefile is in the script 'pull_range_info.R'")
}

if(!file.exists(lake_path)){
  stop("You did not specify where the great lakes shapefile is in the script 'pull_range_info.R'")
}


# Read in the continent data
countries <- sf::read_sf(
  cont_path
) %>% 
  dplyr::filter(COUNTRY %in% c("USA", "CAN")) 

# Reducing to Canada and lower 48
countries <- countries[
  -which(
    countries$NAME %in% 
      c(
        "Hawaii",
        "Puerto Rico",
        "United States Virgin Islands",
        "Navassa Island",
        "water/agua/d'eau"
      )
  ),]


countries <- countries %>% 
  group_by(NAME) %>% 
  summarise(SHAPE__Are = sum(SHAPE__Are))


# make alabama and georgia a LITTLE bigger so we dont have a line
#  between florida and these two states. Since this is just for plotting
#  increasing the size of these doesn't really matter all taht much.
countries[countries$NAME %in% c("Georgia", "Alabama"),] <-
  sf::st_buffer(countries[countries$NAME %in% c("Georgia", "Alabama"),], 400)

# new bounding box
new_bound <- c(
  xmin = -15696048,
  ymin = 2819949,
  xmax = -5862282,
  ymax = 11323599 # removing some of the Northern territory & Nunavut
) 

# crop and combine them
countries <- sf::st_crop(
  countries,
  new_bound
) %>% 
  summarise(SHAPE__Are = sum(SHAPE__Are))

# The mammal data is too large for me to pull entirely onto my laptop.
#  So we're going to need to write a SQL query to collect the specific
#  species for our analysis. First step is to get one row of data
#  from this table to determine column names for querying
one_row <- sf::read_sf(
  "D:/TERRESTRIAL_MAMMALS",
  query = "SELECT * FROM TERRESTRIAL_MAMMALS LIMIT 1"
)




# read in the species data to figure out which ones we need to query
unq_spe <- read.csv(
  "./data/species_in_analysis.csv"
)
# unique species

# Given this we need to convert our species names we generated BACK into
#  binomial nomenclature. I made a function that does this translation
#  for us (will only work for this analysis, it's hardcoded).
bspe <- common_to_binomial(
  unq_spe$species
)


# Create the SQL query.  I made a function that will take a character vector
#  and will convert it to the list for a SQL IN statement. Want to see the 
#  columns available?
#
#  run this to see columns available
#one_row <- sf::read_sf(
# "D:/TERRESTRIAL_MAMMALS",
# query = "SELECT * FROM TERRESTRIAL_MAMMALS LIMIT 1"
#)

### IMPORTANT COLUMNS:
# 'binomial': binomial nomenclature. First letter capitalized in genus.
# 'legend': We only want "Extant & Introduced (resident)" or 
#             "Extant (resident)" in this column.



mammal_qry <- paste0(
  "SELECT tm.binomial AS binomial FROM TERRESTRIAL_MAMMALS tm\n",
  "WHERE tm.binomial IN ", sql_IN(bspe$Binomial), "\n",
  "AND tm.legend IN ", sql_IN(c("Extant & Introduced (resident)",
                                "Extant (resident)"))
)

mams <- sf::read_sf("D:/TERRESTRIAL_MAMMALS",
                    query = mammal_qry
)

mams <- sf::st_transform(
  mams,
  sf::st_crs(countries)
)

# add our common species name to this
mams <- dplyr::left_join(
  mams,
  bspe,
  by = c("binomial" = "Binomial")
)

mams <- mams %>%
  dplyr::group_by(Species) %>% 
  dplyr::summarize(geometry = sf::st_union(geometry)) 

# calculate area after that.
mams$Area <- sf::st_area(
  sf::st_geometry(
    mams
  )[1:nrow(mams),]
)

mams <- sf::st_buffer(
  mams[1:nrow(mams),],
  0
)

# It's currently failing to save as a shapefile, so saving as an RDS file
#  in case my computer stops working at some point.
saveRDS(
  mams,
  "./data/cleaned_data/species_range_maps.RDS"
)

# get the unique sites from city data
sites <- read.csv("./data/gentri_all_coordinates.csv")

sites <- sf::st_as_sf(
  sites, 
  coords = c(3:4),
  crs = 4326
)

sites <- sf::st_transform(
  sites,
  sf::st_crs(countries)
)

# get mean sampled area of each city
city_location <- aggregate(
  sites,
  list(sites$City),
  mean
)

city_location <- data.frame(
  st_coordinates(sites),
  City = sites$City
) %>% 
  dplyr::group_by(City) %>% 
  dplyr::summarise(X = mean(X),
                   Y = mean(Y)) %>% 
  sf::st_as_sf(.,
               coords = 2:3,
               crs = st_crs(countries))

# determine which cities are in a species range
in_city <- sf::st_intersects(
  st_buffer(city_location, 5000),
  mams[1:nrow(mams),],
  sparse = FALSE
) %>% 
  t # transposing to species on rows

# add city names in rows
colnames(in_city) <- tolower(city_location$City)
row.names(in_city) <- mams$Species


write.csv(
  in_city,
  "./data/cleaned_data/covariates/city_in_species_range.csv"
)
# we're going to calculate the distance to the nearest edge of the range.
#  if the city is outside of the range we are going to make that number
#  negative. HOWEVER, the edge of the great lakes will get selected often
#  For Chicago, which we do not want. Will need to rectify that. This is
#  done better if you do each lake individually.

lakes <- rbind(sf::read_sf(
  paste0(lake_path, "/hydro_p_LakeErie")
),
sf::read_sf(
  paste0(lake_path, "/hydro_p_LakeHuron")
),
sf::read_sf(
  paste0(lake_path, "/hydro_p_LakeMichigan")
),
sf::read_sf(
  paste0(lake_path, "/hydro_p_LakeOntario")
),
sf::read_sf(
  paste0(lake_path, "/hydro_p_LakeStClair")
),
sf::read_sf(
  paste0(lake_path, "/hydro_p_LakeSuperior")
)
)

# fix name of st clair lake
lakes$NAMEEN[grep("Clair", lakes$NAMEEN)] <- "Lake St. Clair"
lakes <- lakes %>% 
  dplyr::group_by(NAMEEN) %>% 
  dplyr::summarize(geometry = sf::st_union(geometry)) 

lakes <- st_transform(
  lakes,
  st_crs(countries)
)

# add a small buffer
lakes <- sf::st_buffer(lakes, 2000)
# fill all of the holes in the lakes
lakes <- smoothr::fill_holes(lakes, 1e20)

mam_list <- vector("list", length = nrow(mams))

for(i in 1:nrow(mams)){
  print(i)
  # get one species
  tmp_mams <- mams[i,]
  for(j in 1:nrow(lakes)){
    # convert the lake to a line
    tmp_lake <- sf::st_cast(lakes[j,], "LINESTRING")[1,]
    # check for overlap
    does_intersect <- sf::st_intersects(
      sf::st_geometry(tmp_mams),
      sf::st_geometry(tmp_lake),
      sparse = FALSE
    )
    if(does_intersect){
      
      # get just the overlapping part.
      tmp_lake <- sf::st_intersection(tmp_lake, tmp_mams)
      
      # We then need to make the line back into a polygon. However,
      #  there are times when parts of the line have too few points
      #  (less than 4). So we are going to remove those when 
      #  they occur. It was rather rare in this dataset.
      points_on_line <- sf::st_cast(tmp_lake, "LINESTRING") %>% 
        mapview::npts(., TRUE)
      if(any(points_on_line < 4)){
        tmp_lake <- sf::st_cast(tmp_lake, "LINESTRING")
        tmp_lake <- tmp_lake[-which(points_on_line < 4),]
        tmp_lake <- st_combine(tmp_lake)
      }
      if(as.numeric(st_length(tmp_lake)) > as.numeric( sum(st_length(sf::st_cast(lakes[j,], "LINESTRING"))) * 0.9)){
        # if it's most of the lake just use the lake polygon
        tmp_lake <- st_buffer(lakes[j,], 0)
      }
      
      tmp_lake <- sf::st_cast(tmp_lake, "POLYGON")
      tmp_lake <- sf::st_make_valid(tmp_lake)
      tmp_mams <- sf::st_union(tmp_mams, tmp_lake)
      
      # fill any potential holes
      
      
    } 
  }
  if(any(class(st_geometry(tmp_mams)) == "sfc_GEOMETRYCOLLECTION")){
    tmp_mams <- st_collection_extract(tmp_mams, "POLYGON")
  }
  
  tmp_mams <- fill_holes(tmp_mams, 1e10)
  
  mam_list[[i]] <- tmp_mams[, c("Species", "geometry", "Area")]
  
}
# some of these need to be combined!
for(i in 1:length(mam_list)){
  if(nrow(mam_list[[i]])>1){
    mam_list[[i]] <- mam_list[[i]] %>% 
      dplyr::group_by(Species) %>% 
      dplyr::summarise(
        geometry = sf::st_union(geometry))
  }
  mam_list[[i]]$Area <- sf::st_area(mam_list[[i]])
}


mam_filled_range <- do.call("rbind", mam_list)

# saving this object as it takes some time to run.
saveRDS(
  mam_filled_range,
  "./data/cleaned_data/species_range_greatlakes_filled.RDS"
)
dist_to_edge <- sf::st_distance(
  city_location,
  sf::st_cast(
    mam_filled_range[1:nrow(mams),],
    "MULTILINESTRING"
  )
) %>% 
  t # transpose so species on rows

# now multiply these by negative 1 if they are outside the range
for(i in 1:nrow(dist_to_edge)){
  tmp <- dist_to_edge[i,]
  tmp2 <- as.numeric(in_city[i,])
  if(any(tmp2 == 0)){
    tmp2[tmp2==0] <- -1
  }
  
  tmp3 <- tmp * tmp2
  dist_to_edge[i,] <- tmp3
  
}
# transpose and divide by 1000 to make it km
dist_to_edge <- dist_to_edge/1000

colnames(dist_to_edge) <- city_location$City
row.names(dist_to_edge) <- mams$Species

# save dist_to_edge
write.csv(
  dist_to_edge,
  "./data/cleaned_data/covariates/dist_city_to_species_edge.csv"
)

# plot them out
sp <- mams$Species
for(i in 1:nrow(mams)){
  tmp <- mams[i,]
  tmp <- sf::st_buffer(tmp, 0)
  tmp <- sf::st_intersection(countries, tmp)
  
  jpeg(
    paste0("./plots/species_ranges/",sp[i],".jpeg"),
    height = 550,
    width = 550
  ) 
  par(mar = c(0,0,0,0))
  plot(countries, reset = FALSE , col = "white", main = sp[i])
  plot(sf::st_geometry(tmp), add = TRUE, col = "gray50")
  
  cols <- rep("black", ncol(dist_to_edge))
  cols[as.numeric(dist_to_edge[i,])<0] <- "red"
  
  plot(sf::st_geometry(city_location), 
       add = TRUE, bg = cols,
       pch = 21, cex = 1.5)
  legend("bottomleft", pch = 21, pt.bg = c("black", "red"),
         legend = c("City in species range",
                    "City outside of species range"),
         bty="n"
  )
  dev.off()
  
}




