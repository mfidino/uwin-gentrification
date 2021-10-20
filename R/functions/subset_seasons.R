# Subset seasons

# This function figures out which seasonal "bundles" are present for each
#  city in order to accomodate variation in sampling among seasons and cities.
subset_seasons <- function(x){
  
  season_density <- t(
    table(
      x$Season,
      x$City
    )
  )

  # make this binary, 1 means that city has data on a given season.
  season_density[season_density>0] <- 1
  
  my_seasons <- colnames(season_density)
  
  # get first and last for each city
  season_qry <- apply(
    season_density,
    1,
    function(x) which(x == 1)
  )
  keep_all <- sapply(
    season_qry,
    function(x) all(diff(x) == 1)
  )
  # loop through and figure out which seasons we need to be surgical with
  response <- vector("list", length(keep_all))
  for(i in seq_len(length(keep_all))){
    if(keep_all[i]) {
      response[[i]] <- data.frame(
        City = names(keep_all[i]),
        Season = names(season_qry[[i]])[1],
        n = length(season_qry[[i]]),
        MultipleStart = FALSE
      )
      next
    }
    
    missing <- diff(season_qry[[i]])
    # tricky way to figure out the different season segments
    all_possible <- names(season_qry[[i]])
    new_start <- c(names(season_qry[[i]])[1],  
                   names(missing)[missing> 1])
    to_rep <- diff(which(all_possible %in% new_start))
    # check the last season
    if(new_start[length(new_start)] == all_possible[length(all_possible)]){
      to_rep <- c(to_rep, 1)
    } else {
      to_rep <- c(
        to_rep,
        (length(all_possible) -
          which(all_possible == new_start[length(new_start)])) + 1
      )
    }

    # find where the reset needs to happen
    response[[i]] <- data.frame(
      City = names(keep_all)[i],
      Season = c(names(season_qry[[i]])[1],  
                 names(missing)[missing> 1]),
      n = to_rep,
      MultipleStart = TRUE
    )
  }
  response <- dplyr::bind_rows(response)
  
  return(response)

}
