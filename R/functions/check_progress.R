check_progress_mcmc <- function(outfile, nchain){
  
  is_running <- TRUE
  iter = 1
  total_run <- (nchain * 57) - nchain * 2
  pb <- txtProgressBar(0, total_run)
  while(is_running){
    tmp <- readLines(outfile, warn = FALSE)
    tmp <- tmp[grep("^\\|-",tmp)]
    if(length(tmp) == 0){
      cat("Nimble model getting set up. Checking again in 30 seconds...\n")
      Sys.sleep(30)
      next
    } else {
      nran <- nchar(tmp[length(tmp)])
      setTxtProgressBar(pb, nran)
      Sys.sleep(10)
        
      if(substr(tmp[length(tmp)],nran,nran) == "|"){
        is_running <- FALSE
        cat("\n Model running complete!\n")
      }
    }
  }
}

check_progress_mcmc(outfile = "hmmm.txt", nchain = 3)
