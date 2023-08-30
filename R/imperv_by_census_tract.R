library(dplyr)
library(sf)
library(FedData)
library(raster)

cities <- readRDS(
  "./data/census_data/gent_census_tracts.rds"
)

yrs <- c(2011,2019)

city_result <- vector("list", length = length(cities))

for(i in 1:length(cities)){
  print(
    paste(
      i,"of",length(cities)
    )
  )
  for(yr in 1:length(yrs)){
    
    dat <- cities[[i]]
    
    my_bbox <- sf::st_bbox(dat)[c("xmin","xmax","ymin","ymax")]
    
    my_polygon <- FedData::polygon_from_extent(
      raster::extent(
        as.numeric(
          my_bbox +
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
      ),
      force.redo = TRUE
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
        na.rm = TRUE
      )
    )
    
    #
    if(yr == 1){
      city_result[[i]] <- data.frame(
        city = cities[[i]]$City,
        gent = cities[[i]]$gentrifying,
        vuln = cities[[i]]$vulnerable,
        mean_11 = prop_extract
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
  }
}
# grab impervious for each year
library(lme4)
library(dplyr)
library(brms)

#saveRDS(city_result, "imp_census_tract.RDS")

city_result <- readRDS("imp_census_tract.RDS")

cc <- do.call("rbind", city_result)
cc <- cc[complete.cases(cc),]


tmp <- cc %>% 
  dplyr::group_by(city, gent, vuln) %>%
  dplyr::summarise(
    mu11 = mean(mean_11),
    mu19 = mean(mean_19)
  ) %>% 
  data.frame()




