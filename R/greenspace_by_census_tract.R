library(dplyr)
library(sf)
library(FedData)
library(raster)
library(uwinspatialtools)

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
      dataset = "landcover",
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
    spatial_summary <- function(x, ncats = my_raster@data@max, 
                                ...) {
      return(prop.table(tabulate(x, ncats)))
    }
    
    prop_extract <- suppressWarnings(
      raster::extract(my_raster, 
                        dat,
                      fun = spatial_summary,
                      buffer = my_buffer, 
                      na.rm = TRUE))
    # get 21 (developed open space)
    prop_extract <- prop_extract[,21]
    
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

saveRDS(city_result, "gs_census_tract.RDS")

city_result <- readRDS("gs_census_tract.RDS")

cc <- do.call("rbind", city_result)
cc <- cc[complete.cases(cc),]


cc$diff <- cc$mean_19 - cc$mean_11


cc <- cc %>% 
  dplyr::group_by(city, gent) %>% 
  dplyr::summarise(
    mean_diff = round(mean(diff),2),
    sd_diff = round(sd(diff),2),
    min = round(min(diff),2),
    max = round(max(diff),2)
) %>% 
  data.frame()

cc$imp <- cc$imp / 25



cc$city <- factor(cc$city)

m1 <- glmer(
  gent ~ imp + (1 + imp| city),
  data = cc,
  family = "binomial"
)

m2 <- brm(
  gent ~ imp + (1 + imp|city),
  data = cc,
  family = "bernoulli"
)

summary(m1)

imp <- seq(0,100, length.out = 100)

pred_mat <- data.frame(
  imp = (imp - mean(imp))/25
)

my_pred <- predict(m1, newdata = pred_mat, re.form=NA, type = "response")




cc$bin <- cut(
  cc$mean_19,
  breaks = seq(0,100, 10),
  include.lowest = TRUE
)

yo <- cc %>% 
  dplyr::group_by(bin) %>% 
  dplyr::summarise(
    gent_mu = mean(gent),
    gent_se = sd(gent) / sqrt(length(gent))
)
yo$mid <- seq(5,95,10)

plot(
  my_pred ~ imp, type = "l", ylim = c(0,1),
  xlab = "Impervious cover within census tract (%)",
  ylab = "Probability census tract is gentrifying",
  bty = "l",
  las = 1,
  lwd = 3
)
for(i in 1:nrow(yo)){
  lines(
    x = rep(yo$mid[i],2),
    y = c(
      yo$gent_mu[i] + yo$gent_se[i],
      yo$gent_mu[i] - yo$gent_se[i]
    ),
    lwd = 2
  )
}
points(y = yo$gent_mu, x = yo$mid, pch = 21, cex = 1.3,
       bg = "gray50")





