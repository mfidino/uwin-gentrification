convert_to_logmean <- function(mean, sd) {
  log(mean^2 / sqrt(sd^2 + mean^2))
}

convert_to_logsd <- function(mean, sd) {
  sqrt(log(1 + (sd^2 / mean^2)))
}

# function to split up mcsamp into a list with correctly shaped arrays
split_mcmc <- function(x){
  # get parameter names
  pars <- colnames(x)
  # unique parameters
  unq_pars <- unique(
    gsub(
      "\\[.*\\]",
      "",
      pars
    )
  )
  # make list object to store arrays in
  result_list <- vector(
    "list",
    length = length(unq_pars)
  )
  names(result_list) <- unq_pars
  # fill in the arrays
  for(i in 1:length(result_list)){
    # get just the parameters
    tmp <- pars[grep(
      paste0(
        "^",unq_pars[i], "\\["
      ),
      pars
    )]
    if(length(tmp) == 0){
      tmp <- pars[grep(
        paste0("^",unq_pars[i],"$"),
        pars
      )]
    }
    # and then the array dimensions
    arr_dim <- gsub(
      paste0(
        unq_pars[i],"\\[|\\]"
      ),
      "",
      tmp
    )
    arr_dim <- strsplit(
      arr_dim,
      ","
    )
    ndim <- length(arr_dim[[1]])
    npar <- length(arr_dim)
    # make a matrix
    arr_ind <- suppressWarnings(
      matrix(
        as.numeric(
          unlist(arr_dim)
        ),
        ncol = ndim,
        nrow = npar,
        byrow = TRUE
      )
    )
    if(nrow(arr_ind) == 1 & ncol(arr_ind) == 1){
      arr_ind[1,1] <- 1
    }
    # get max index for each
    max_ind <- apply(arr_ind, 2, max)
    # and then fill in the array
    result_list[[i]] <- array(
      x[,tmp],
      dim = c(nrow(x), max_ind)
    )
    
  }
  return(result_list)
}

qwrap <- function(x){
  quantile(x, probs = c(0.025,0.5,0.975))
}
