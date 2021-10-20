order_seasons <- function(x){
  unq_seasons <- unique(x)
  month <- substr(unq_seasons,1,2)
  
  tmp <- sapply(month, function(y) switch(y, "JA" = 0.1,
                                             "AP" = 0.2,
                                             "JU" = 0.3,
                                             "OC" = 0.4)
          )
  
  year <- as.numeric(paste0("20", substr(unq_seasons, 3,4)))
  
  to_order <- year + tmp
  
  unq_seasons <- unq_seasons[order(to_order)]
  
  return(unq_seasons)
  
}
