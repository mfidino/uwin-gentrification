# combine all the site coordinates
library(dplyr)

# some of the data is already ready to go.

dat <- read.csv("../uwin-dataset/data_2022-04-11/cleaned_data/full_capture_history.csv")

dat <- read.csv("./data/detection_data.csv")

dat$Season <- factor(dat$Season, unique(dat$Season))

dat$City <- gsub("safa|oaca", "baca", dat$City)


dat$City <- gsub("phaz2", "phaz", dat$City)
dat <- dat[dat$City %in% dat2$City,]

dat <- dat[dat$Season %in% dat2$Season,]

# get number of sites sampled per city
sampled <- dat[dat$Species == dat$Species[1],] %>%
  dplyr::group_by(Season, City) %>% 
  dplyr::summarise(ssamp = sum(J>0))
sampled <- sampled[order(sampled$City, sampled$Season),]

sampled <- split(
  sampled,
  factor(sampled$City)
)

# turn all 10 or less to zero as well

for(i in 1:length(sampled)){
  sampled[[i]]$ssamp[sampled[[i]]$ssamp <= 10] <- 0
}

sampled <- dplyr::bind_rows(sampled)
# go through and convert all 0 for a city / season
for(i in 1:nrow(sampled)){
  if( sampled$ssamp[i] == 0 ){
    dat$Y[
      dat$City == sampled$City[i] & 
      dat$Season == sampled$Season[i]
    ] <- 0
    dat$J[
      dat$City == sampled$City[i] & 
      dat$Season == sampled$Season[i]
    ] <- 0
  }
}

# remove sites sampled zero times

si_go <- dat %>% 
  group_by(Site, City) %>% 
  summarise(nsamp = sum(J)) %>% 
  filter(nsamp == 0) %>% 
  data.frame()

dat <- dat[-which(dat$Site %in% si_go$Site),]


write.csv(dat, "./data/detection_data.csv")


coords <- dat[,c("Site", "City","Long", "Lat", "Crs")]

coords <- coords %>% group_by(Site,City) %>% 
  summarise(Long = mean(Long), Lat = mean(Lat)) %>% 
  as.data.frame()
row.names(coords) <- NULL

# do a histogram and inspect
hist(coords$Lat)
coords[coords$Lat < 30,]

ctmp <- split(
  coords,
  factor(coords$City)
)

# t1 <- sf::st_as_sf(
#   data.frame(x = 720932, y = 4270280),
#   coords = c("x","y"),
#   crs = 32615
# )
# 
# t1 <- sf::st_transform(
#   t1,
#   crs = 4326
# )
# 
# hist(coords$Long)
write.csv(coords, "./data/gentri_all_coordinates.csv", row.names = FALSE)
