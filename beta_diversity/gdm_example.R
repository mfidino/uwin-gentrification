library(gdm)

# reads in example input data
load(system.file("./data/gdm.RData", package="gdm"))
# columns 3-7 are soils variables, remainder are climate
gdmExpData[1:3,]
# get columns with xy, site ID, and species data
sppTab <- gdmExpData[, c("species", "site", "Lat", "Long")]

# Note: It's in long format

# get columns with env. data and xy-coordinates
envTab <- gdmExpData[, c(2:ncol(gdmExpData))]


# x-y species list example
gdmTab <- formatsitepair(
  sppTab, bioFormat=2, XColumn="Long", YColumn="Lat",
  sppColumn="species", siteColumn="site", predData=envTab
)
gdmTab[1:3,]

hm <- formatsitepair(
  data.frame(
    site = ed$site, z
    ),
  bioFormat = 1,
  dist = "bray",
  siteColumn = "site",
  XColumn = "long",
  YColumn = "lat",
  predData = ed
)


site <- unique(sppTab$site)
gdmDissim <- cbind(site, gdmDissim)

# Fit first model
gdm.1 <- gdm(gdmTab, geo=T)
gdm.2 <- gdm(hm, geo=TRUE)
summary(gdm.1)

length(gdm.1$predictors) # get idea of number of panels

windows(12,4)
plot(gdm.2, plot.layout=c(1,3))

mysd <- isplineExtract(gdm.1)

# x = seq of values from minimum to maximum value for each
#       covariate.

my_goal <- splinedat$y[,3]
head(my_goal)

# get my knots

knotties <- quantile(
  c(gdmTab$s1.phTotal, gdmTab$s2.phTotal),
  probs = c(0,0.5,1)
)
my_mcmc <- c(1.12736848908874,0.229764248308982,0)

longshot <- spline_pred(200, knotties, my_mcmc)

# assuming 3 splines per predictor






