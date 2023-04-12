library(sf)
sf::sf_use_s2(FALSE)
library(uwinspatialtools)
library(dplyr)


# read in income data that was already queried
price_list <- readRDS(
  "./data/census_data/med_income/gent_classification_r2.rds"
)

for(i in 1:length(price_list)){
  price_list[[i]] <- sf::st_transform(
    price_list[[i]],
    crs = 4326
  )
}

income_list <- vector("list", length = length(price_list))
for(i in 1:length(income_list)){
  income_list[[i]] <- sf::st_join(
    price_list[[i]],
    gent_list_poly[[i]],
    largest = TRUE
  )
}

income_df <- as.data.frame(
  dplyr::bind_rows(
    income_list
  )
)

overall <- income_df[!is.na(income_df$gentrifying),] %>% 
  dplyr::group_by(gentrifying) %>% 
  dplyr::summarise(
    mu = mean(estimate, na.rm = TRUE),
    sd = sd(estimate, na.rm = TRUE)
  )


income_list <- vector("list", length = length(price_list))
for(i in 1:length(income_list)){
  income_list[[i]] <- sf::st_join(
    price_list[[i]],
    gent_list_poly[[i]],
    largest = TRUE
  )
}

income_df <- as.data.frame(
  dplyr::bind_rows(
    income_list
  )
)

overall <- income_df[!is.na(income_df$gentrifying),] %>% 
  dplyr::group_by(gentrifying) %>% 
  dplyr::summarise(
    mu = mean(estimate, na.rm = TRUE),
    sd = sd(estimate, na.rm = TRUE)
  )


race_list <- readRDS("./data/census_data/race/gent_classification_r2.rds")


for(i in 1:length(race_list)){
  race_list[[i]] <- sf::st_transform(
    race_list[[i]],
    crs = 4326
  )
}

race_results <- vector("list", length = length(price_list))
for(i in 1:length(income_list)){
  race_results[[i]] <- sf::st_join(
    race_list[[i]],
    gent_list_poly[[i]],
    largest = TRUE
  )
}
race_df <- as.data.frame(
  dplyr::bind_rows(
    race_results
  )
)



overall_race <- race_df[!is.na(race_df$gentrifying),] %>% 
  dplyr::group_by(gentrifying) %>% 
  dplyr::summarise(
    mu = mean((value/total), na.rm = TRUE),
    sd = sd((value/total), na.rm = TRUE)
  )


edu_list <- readRDS(
  "./data/census_data/education/gent_classification_r2.rds"
)


for(i in 1:length(edu_list)){
  edu_list[[i]] <- sf::st_transform(
    edu_list[[i]],
    crs = 4326
  )
}

edu_results <- vector("list", length = length(price_list))
for(i in 1:length(edu_list)){
  edu_results[[i]] <- sf::st_join(
    edu_list[[i]],
    gent_list_poly[[i]],
    largest = TRUE
  )
}

edu_df <- as.data.frame(
  dplyr::bind_rows(
    edu_results
  )
)


overall_edu <- edu_df[!is.na(edu_df$gentrifying),] %>% 
  dplyr::group_by(gentrifying) %>% 
  dplyr::summarise(
    mu = mean((value/total), na.rm = TRUE),
    sd = sd((value/total), na.rm = TRUE)
  )
