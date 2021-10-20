make_index <- function(x, season = FALSE){
  if(season){
    as.numeric(
      factor(x,
             levels = order_seasons(x)
      )
    )
  } else {
    as.numeric(
      factor(
        x
      )
    )
  }
}
