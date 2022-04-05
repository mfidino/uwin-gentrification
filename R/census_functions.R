# A function to read in and bind the data that
#  was previously pulled from the Census API.
#  When saved, it was still a list object with
#  a different sf for each county we need data
#  from.
read_bind <- function(x){
  tmp <- lapply(
    x,
    readRDS
  )
  for(i in 1:length(tmp)){
    for(j in 1:length(tmp[[i]])){
      if(dplyr::is.grouped_df(
        tmp[[i]][[j]]
        )
      ){
        tmp[[i]][[j]] <- dplyr::ungroup(tmp[[i]][[j]])
      }
    }
  }
  tmp <- lapply(
    tmp,
    dplyr::bind_rows
  )
  names(tmp) <- c(2000,2010,2015,2019)
  return(tmp)
}

# going from long/lat to utm crs
longlat_to_utm <- function(longlat) {
  if("sf" %in% class(longlat)){
    longlat <- c(
      mean(sf::st_bbox(longlat)[c("xmin", "xmax")]),
      mean(sf::st_bbox(longlat)[c("ymin", "ymax")])
    )
  }
  utm <- (floor((longlat[1] + 180) / 6) %%60) + 1
  if(longlat[2] > 0) {
    utm <- utm + 32600
  } else{
    utm <- utm + 32700
  }
  return(as.numeric(utm))
}



# go through the decennial and acs data to get bounding
#  box for each census tract. Then take average of
#  long / lat to get a point to convert to utms.

bbox_to_point <- function(x){
  tmp <- sf::st_bbox(x)
  tmp <- c(
    mean(
      tmp[c(1,3)]
    ),
    mean(
      tmp[c(2,4)]
    )
  )
  return(tmp)
}

add_utms <- function(df, id = "geoid"){
  tmp <- df[,id]
  tmp <- tmp[!duplicated(tmp),]
  
  cat("\n\ncompiling crs...\n\n")
  
  # go through a loop and get UTMS
  pb <- txtProgressBar(max = nrow(tmp))
  my_utms <- rep(NA, nrow(tmp))
  for(i in 1:length(my_utms)){
    response <- try(
      longlat_to_utm(
        bbox_to_point(
          tmp[i,]
        )
      ),
      silent = TRUE
    )
    if(is.numeric(response)){
      my_utms[i] <- response
    }
    setTxtProgressBar(pb, i)
  }
  
  # remove NAs
  to_go <- which(is.na(my_utms))
  if(length(to_go)>0){
    tmp <- tmp[-to_go,]
    tmp$crs <- my_utms[-to_go]
  }
  tmp <- data.frame(
    geoid = tmp$geoid,
    crs = tmp$crs
  )
  
  
  df <- dplyr::inner_join(
    df,
    tmp,
    by = id
  )
  return(df)
}


interpolate_census_data <- function(
  my_cdat,
  buff_coords
){
  if(!is.list(my_cdat)){
    stop('my_cdat needs to be a list')
  }
  if(!is.list(buff_coords)){
    stop('buff_coords needs to be a list')
  }
  if(!all(names(my_cdat[[1]]) == names(buff_coords))){
    stop("both data sets need to be split() by utm crs")
  }
  
  tmp_cdat <- my_cdat
  nyear <- length(my_cdat)
  # number of projections
  ncrs <- length(buff_coords)
  results <- vector("list", length = nyear)
  
  for(yr in 1:nyear){
    cat(
      paste0("\nyear ", yr ," of ", nyear,"\n")
    )
    results[[yr]] <- vector("list", length = ncrs)
    pb <- txtProgressBar(max = ncrs)
    for(cs in 1:ncrs){
      setTxtProgressBar(pb, cs)
      nsite <- nrow(buff_coords[[cs]])
      results[[yr]][[cs]] <- vector('list', length = nsite)
      # gets the tracts that intersect with each point
      #  rows = census tracts, cols = sites
      tmp1 <- sf::st_intersects(
        my_cdat[[yr]][[cs]],
        buff_coords[[cs]],
        sparse = FALSE
      )

      # get area for each bit of census info
      colnames(my_cdat[[yr]][[cs]]) <- tolower(colnames(my_cdat[[yr]][[cs]]))
      for(site in 1:nsite){
        total_area <- sf::st_area(
        my_cdat[[yr]][[cs]][tmp1[,site],]
        )
        suppressWarnings(
          tmp2 <- sf::st_intersection(
            my_cdat[[yr]][[cs]],
            buff_coords[[cs]][site,]
          )
        )
      
        ins_area <- sf::st_area(
          tmp2
        )
        census_weight <- ins_area / total_area
        
        if(all(c("estimate", "summary_est") %in% colnames(my_cdat[[yr]][[cs]]))){
          tmp_cdat[[yr]][[cs]]$estimate[tmp1[,site]] <- 
            census_weight *  my_cdat[[yr]][[cs]]$estimate[tmp1[,site]] 
          tmp_cdat[[yr]][[cs]]$summary_est[tmp1[,site]] <-
            census_weight * my_cdat[[yr]][[cs]]$summary_est[tmp1[,site]]
          
          # group and summarise based on the variable name
          
          results[[yr]][[cs]][[site]] <- tmp_cdat[[yr]][[cs]][tmp1[,site],] %>% 
            dplyr::group_by(variable) %>% 
            dplyr::summarise(
              value = sum(estimate, na.rm = TRUE),
              total = sum(summary_est, na.rm = TRUE))
          results[[yr]][[cs]][[site]]$site <- buff_coords[[cs]]$Site[site]
          results[[yr]][[cs]][[site]]$city <- buff_coords[[cs]]$City[site]
          results[[yr]][[cs]][[site]] <- data.frame(
            results[[yr]][[cs]][[site]]
          )[,c("city", "site", "variable", "value", "total")]
          
        }
        if("estimate" %in% colnames(my_cdat[[yr]][[cs]]) &
           !"summary_est" %in% colnames(my_cdat[[yr]][[cs]])){
          tmp_cdat[[yr]][[cs]]$estimate[tmp1[,site]] <- 
            census_weight *  my_cdat[[yr]][[cs]]$estimate[tmp1[,site]] 

          # group and summarise based on the variable name
          
          results[[yr]][[cs]][[site]] <- tmp_cdat[[yr]][[cs]][tmp1[,site],] %>% 
            dplyr::group_by(variable) %>% 
            dplyr::summarise(
              value = sum(estimate, na.rm = TRUE)
            )
          results[[yr]][[cs]][[site]]$site <- buff_coords[[cs]]$Site[site]
          results[[yr]][[cs]][[site]]$city <- buff_coords[[cs]]$City[site]
          results[[yr]][[cs]][[site]] <- data.frame(
            results[[yr]][[cs]][[site]]
          )[,c("city", "site", "variable", "value")]
          
        }
        
        if("summary_est" %in% colnames(my_cdat[[yr]][[cs]]) & 
           "hs diploma" %in% my_cdat[[yr]][[cs]]$variable){
          tmp_cdat[[yr]][[cs]]$value[tmp1[,site]] <- 
            census_weight *  my_cdat[[yr]][[cs]]$value[tmp1[,site]] 
          tmp_cdat[[yr]][[cs]]$summary_est[tmp1[,site]] <-
            census_weight * my_cdat[[yr]][[cs]]$summary_est[tmp1[,site]]
          
          # group and summarise based on the variable name
          
          results[[yr]][[cs]][[site]] <- tmp_cdat[[yr]][[cs]][tmp1[,site],] %>% 
            dplyr::group_by(variable) %>% 
            dplyr::summarise(
              value = sum(value, na.rm = TRUE),
              total = sum(summary_est, na.rm = TRUE))
          results[[yr]][[cs]][[site]]$site <- buff_coords[[cs]]$Site[site]
          results[[yr]][[cs]][[site]]$city <- buff_coords[[cs]]$City[site]
          results[[yr]][[cs]][[site]] <- data.frame(
            results[[yr]][[cs]][[site]]
          )[,c("city", "site", "variable", "value", "total")]
          
        }
        if("summary_value...total" %in% colnames(my_cdat[[yr]][[cs]])){
        
          tmp_cdat[[yr]][[cs]]$value[tmp1[,site]] <- 
            census_weight *  my_cdat[[yr]][[cs]]$value[tmp1[,site]] 
          tmp_cdat[[yr]][[cs]]$summary_value...total[tmp1[,site]] <-
            census_weight * my_cdat[[yr]][[cs]]$summary_value...total[tmp1[,site]]
          
          # group and summarise based on the variable name
          
          results[[yr]][[cs]][[site]] <- tmp_cdat[[yr]][[cs]][tmp1[,site],] %>% 
            dplyr::group_by(variable) %>% 
            dplyr::summarise(
              value = sum(value, na.rm = TRUE),
              total = sum(summary_value...total, na.rm = TRUE))
          results[[yr]][[cs]][[site]]$site <- buff_coords[[cs]]$Site[site]
          results[[yr]][[cs]][[site]]$city <- buff_coords[[cs]]$City[site]
          results[[yr]][[cs]][[site]] <- data.frame(
            results[[yr]][[cs]][[site]]
          )[,c("city", "site", "variable", "value", "total")]
        }
        if(length(colnames(my_cdat[[yr]][[cs]])) == 5 &
           "value" %in% colnames(my_cdat[[yr]][[cs]])){
          
          tmp_cdat[[yr]][[cs]]$value[tmp1[,site]] <- 
            census_weight *  my_cdat[[yr]][[cs]]$value[tmp1[,site]] 
          
          # group and summarise based on the variable name
          
          results[[yr]][[cs]][[site]] <- tmp_cdat[[yr]][[cs]][tmp1[,site],] %>% 
            dplyr::group_by(variable) %>% 
            dplyr::summarise(
              value = sum(value, na.rm = TRUE))
          results[[yr]][[cs]][[site]]$site <- buff_coords[[cs]]$Site[site]
          results[[yr]][[cs]][[site]]$city <- buff_coords[[cs]]$City[site]
          results[[yr]][[cs]][[site]] <- data.frame(
            results[[yr]][[cs]][[site]]
          )[,c("city", "site", "variable", "value")]
        }
      }
    }
  }
  tmp_list <- vector("list", length = nyear)
  to_return <- vector("list", length = nyear)
  # fix the results a bit
  for(yr in 1:nyear){
    tmp_list[[yr]] <- vector("list", length = ncrs)
    for(cs in 1:ncrs){
    tmp_list[[yr]][[cs]] <- do.call("rbind", results[[yr]][[cs]])
    tmp_list[[yr]][[cs]]$year <- yr
    }
    to_return[[yr]] <- do.call("rbind", tmp_list[[yr]])
  }
  return(
    do.call("rbind", to_return)
  )
  
}
