library(censusapi)
library(dplyr)

# I think tidycensus may be the way to go, seems like
# we can get educational attainment as well from these,
# we just need ot know which ones to grab.


# Note: you need to get an api from the census website and then use

# Add key to .Renviron
 Sys.setenv(CENSUS_KEY="e962407457f5dada3d25db6173c82221336c12a2")
# Reload .Renviron
 readRenviron("~/.Renviron")
# Check to see that the expected key is output in your R console
 Sys.getenv("CENSUS_KEY")

##########
# controls
##########

# open url to api census metadata?
open_api_site <- TRUE

hm <- load_variables(2010, "sf1")

# step 1. locate the names of the data you need.
vars10 <- c("P005003", "P005004", "P005006", "P004003")
il <- get_decennial(geography = "tract", variables = vars10, year = 2010,
                    summary_var = "P001001", state = "IL", county = "Cook",
                    geometry = TRUE) %>%
  mutate(pct = 100 * (value / summary_value))

# Current website link:
if(open_api_site){
  browseURL("https://api.census.gov/data.html")
}

 apis <- listCensusApis()
 View(apis[apis$vintage == 2010,])
 

 
# Get 2010 data 
 
 # check geography
 dec_geo <- listCensusMetadata(
   name = "dec/cd113profile",
   vintage = 2010,
   type = "geography"
 )


dec_vars <- listCensusMetadata(
  name = "dec/cd113profile",
  vintage = 2010,
  type = "variables"
)

View(dec_vars[grep("incom", TRUE, dec_vars$label),])

# try it out
longshot <- censusapi::getCensus(
  name = "dec/cd110sprofile",
  vars = "PCO008033",
  vintage = 2010,
  region = "state:01"
)

getCensus(name = "timeseries/healthins/sahie",
          vars = c("NAME", "IPRCAT", "IPR_DESC", "PCTUI_PT"), 
          region = "state:01",
          time = 2018)

data2010 <- getCensus(
  name = "dec/sf1",
  vintage = 2010,
  vars = c("NAME", "P001001", "H010001"), 
  region = "metropolitan statistical area/micropolitan statistical area:*")
head(data2010)

