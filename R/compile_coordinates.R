# combine all the site coordinates
library(dplyr)

# some of the data is already ready to go.

dat <- read.csv("../uwin-dataset/data_2021-07-27/cleaned_data/full_capture_history.csv")


dat <- dat[,c("City", "Site", "Long", "Lat", "Crs")]
dat <- dat[!duplicated(dat),]

cities <- c(
  "ahga",
  "autx",
  "boma",
  "buny",
  "chil",
  "deco",
  "inin",
  "ioio",
  "jams",
  "lrar",
  "lbca",
  "mawi",
  "naca",
  "phaz",
  "phaz2",
  "poor",
  "rony",
  "scut",
  "sewa",
  "slmo",
  "tawa",
  "toon",
  "uril",
  "wide"
)

dat <- dat %>% group_by(Site,City) %>% 
  summarise(Long = mean(Long), Lat = mean(Lat)) %>% 
  as.data.frame()
row.names(dat) <- NULL

# drop sanfrancisco
dat <- dat[-which(dat$City == "safa"),]

# read in coords from uwin db
d2 <- read.csv("./data/raw_coordinates/gent_sites.csv")

# remove any city where I already have coors
d2$areaAbbr <- tolower(d2$areaAbbr)

d2 <- d2[-which(d2$areaAbbr %in% c(dat$City, "paca", "lbca")),]


d2 <- split(d2, factor(d2$areaAbbr))
library(sf)
for(i in 1:length(d2)){
  tmp <- d2[[i]]
  tmp$utmZone <- gsub("T|S|N| ", "", tmp$utmZone)
  tmp <- sf::st_as_sf(
    tmp,
    coords = c("utmEast", "utmNorth"),
    crs = as.numeric(paste0("326", unique(tmp$utmZone)))
  )
  tmp <- sf::st_transform(
    tmp,
    crs = 4326
  )
  
  response <- data.frame(
    Site = tmp$locationAbbr,
    City = tmp$areaAbbr,
    Long = sf::st_coordinates(tmp)[,1],
    Lat = sf::st_coordinates(tmp)[,2]
  )
  d2[[i]] <- response
}

d2 <- dplyr::bind_rows(d2)

d2 <- data.frame(d2)

d2 <- d2[-grep("SFCA-O01-GGP2", d2$Site),]
all_dat <- dplyr::bind_rows(
  list(dat, d2)
)

all_dat <- all_dat[all_dat$City %in% cities,]

which(!cities %in% all_dat$City)
all_dat$Crs <- "4326"
row.names(all_dat) <- NULL
write.csv(all_dat, "gentri_all_coordinates.csv", row.names = FALSE)
