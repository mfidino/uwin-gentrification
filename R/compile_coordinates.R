# combine all the site coordinates
library(dplyr)

# some of the data is already ready to go.

dat <- read.csv("../uwin-dataset/data_2022-04-11/cleaned_data/full_capture_history.csv")

#
dat$City <- gsub("phaz2", "phaz", dat$City)
dat$City <- gsub("oaca|safa", "baca", dat$City)

#dat$City <- gsub("phaz2", "phaz", dat$City)

dat2 <- read.csv("./data/detection_data.csv")

# reduce down to cities in dat2
dat <- dat[dat$City %in% c(unique(dat2$City)),]

# reduce down to correct seasons
dat <- dat[dat$Season %in% unique(dat2$Season),]

# and then species
dat <- dat[dat$Species %in% unique(dat2$Species),]

# change some city names and reorder
dat$Season <- factor(dat$Season, levels = unique(dat$Season))

# drop any sites with absolutely zero data
zdat <- dat %>% 
  dplyr::group_by(Site,City) %>% 
  dplyr::summarise(ssamp = sum(J>0)) %>% 
  dplyr::filter(ssamp == 0)


dat <- dat[-which(dat$Site %in% zdat$Site),]

# drop autx
dat <- dat[-which(dat$City == "autx"),]

write.csv(dat, "./data/detection_data.csv", row.names = FALSE)

# drop some season data here MASON

# get number of sites sampled per city
hm <- dat[dat$Species == dat$Species[1],] %>%
  dplyr::group_by(Season, City) %>% 
  dplyr::summarise(ssamp = sum(J>0))
hm <- hm[order(hm$City, hm$Season),]

hm <- split(
  hm,
  factor(hm$City)
)
# censor all of autx
hm$autx$ssamp <- 0

# turn all 10 or less to zero as well

for(i in 1:length(hm)){
  hm[[i]]$ssamp[hm[[i]]$ssamp <= 10] <- 0
}

write.csv(dat, "./data/detection_data.csv")

sf <- dat[dat$City == "safa",]

sf <- distinct(sf[,c("Site", "Long", "Lat")])

tmp <- SELECT(
  "select cl.locationAbbr, cl.utmEast, cl.utmNorth, cl.utmZone from CameraLocations cl
  inner join StudyAreas sa on sa.areaID = cl.areaID
  and sa.areaAbbr IN ('safa')"
)

library(sf)


tmp <- sf::st_as_sf(
  tmp,
  coords = c("utmEast", "utmNorth"),
  crs = 32610
)

my_df <- data.frame(
  site = tmp$locationAbbr,
  long = st_coordinates(t2)[,1],
  lat = st_coordinates(t2)[,2],
  crs = 4326,
  city = "SAFA"
)

t2 <- sf::st_transform(
  tmp,
  crs = 4326
)
plot(st_coordinates(t2))

aa <- season_availability(dat, "here.tiff", TRUE)

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
