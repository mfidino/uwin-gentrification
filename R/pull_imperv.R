# Query all impervious surface data for UWIN

library(dplyr)
library(sf)
library(FedData)
library(raster)

grab <- function(x,y) x[,grep(y,colnames(x))]


# read in city info
cities <- read.csv(
  "./data/gentri_all_coordinates.csv"
)


cities <- split(cities, factor(cities$City))

yrs <- c(2001, 2004, 2006, 2008, 2011, 2016, 2019)

city_result <- vector("list", length = length(cities))

for(i in 1:length(cities)){
  print(names(cities)[i])
  pb <- txtProgressBar(
    0,
    length(yrs)
  )
  for(yr in 1:length(yrs)){
    
    dat <- cities[[i]]
    
    dat <- sf::st_as_sf(
      dat,
      coords = c("Long", "Lat"),
      crs = 4326
    )
    
    my_polygon <- FedData::polygon_from_extent(
      raster::extent(
        as.numeric(
          sf::st_bbox(dat)[c("xmin","xmax","ymin","ymax")] +
            c(-.05,.05,-.05,.05)
        )
      ),
      proj4string = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    )
    
    my_raster <- FedData::get_nlcd(
      my_polygon,
      dataset = "impervious",
      year = yrs[yr],
      label =paste0(
        names(cities)[i],
        "_",
        yrs[yr]
      )
    )
    
    # reproject points to raster crs
    dat <- sf::st_transform(
      dat,
      crs = crs(my_raster)
    )
    
    # extract data
    prop_extract <- suppressWarnings(
      raster::extract(
        my_raster, 
        dat,
        fun = mean,
        buffer = 1000, 
        na.rm = TRUE
      )
    )
    
    #
    if(yr == 1){
      city_result[[i]] <- data.frame(
        city = cities[[i]]$City,
        site = cities[[i]]$Site,
        mean_01 = prop_extract
      )
    }else{
      city_result[[i]]$new <- prop_extract
      colnames(city_result[[i]])[
        ncol(city_result[[i]])
        ] <- paste0(
          "mean_",
          substr(
            as.character(yrs[yr]),3,4
          )
        )
    }
    setTxtProgressBar(
      pb,
      yr
    )
  }
}
# grab impervious for each year
cc <- do.call("rbind", city_result)

write.csv(cc, "./data/imperv.csv", row.names = FALSE)
