#--
# Functions I cannibalized from the
#  c++ code in gdm.
#--


# This function creates the basis function
#  for the beta diversity analysis.
dospline <- function(dVal = NULL, my_knots = NULL){
  
  if(is.null(my_knots)){
    knots <- c(
      min(dVal),
      median(dVal),
      max(dVal)
    )
  } else {
    knots <- my_knots
  }
  response_mat <- matrix(
    NA,
    ncol = 3,
    nrow = length(dVal)
  )
  for(i in 1:3){
    if(i == 1){
      q1 <- q2 <- knots[i]
      q3 <- knots[i+1]
    }
    if(i == 2){
      q1 <- knots[i-1]
      q2 <- knots[i]
      q3 <- knots[i+1]
    }
    if( i == 3){
      q1 <- knots[i-1]
      q2 <- q3 <- knots[i]
    }
    response <- rep(NA, length(dVal))
    
    if(any(dVal <= q1)){
      response[
        which(dVal <= q1)
        ] <- 0
    }
    if(any(dVal >= q3)){
      response[
        which(dVal >= q3)
        ] <- 1
    }
    if(any(( q1 < dVal ) & ( dVal < q2 ))){
      tmp_loc <- which(( q1 < dVal ) & ( dVal < q2 ))
      response[tmp_loc] <- 
        (
          (
            (dVal[tmp_loc]-q1) * 
              (dVal[tmp_loc]-q1)
          )/
            (
              (q2-q1)*
                (q3-q1)
            )
        )
    }
    if(any(is.na(response))){
      tmp_loc <- is.na(response)
      response[tmp_loc] <- 
        (
          1.0 - 
            ( 
              ( 
                (q3-dVal[tmp_loc]) *
                  (q3-dVal[tmp_loc])
              ) / 
                (
                  (q3-q2) *
                    (q3-q1)
                ) 
            ) 
        ) 
      
    }
    response_mat[,i] <- response
  }
  return(response_mat)
}

# This one generates predictions from the mcmc steps
#  of the model.
spline_pred <- function(
  knots,
  mcmc,
  my_probs = c(0.025,0.5,0.975),
  predlength = 200
  ){
  dVal <- xVal <- seq(
    knots[1],
    knots[3],
    length.out = predlength
  )
  if(is.matrix(mcmc)){
    if(ncol(mcmc)!=3){
      stop("mcmc needs 3 columns (one for each spline).")
    }
  }
  if(!is.matrix(mcmc)){
    mcmc <- as.matrix(mcmc)
  }
  if(nrow(mcmc) != 3){
    mcmc <- t(mcmc)
  }
  # get splines again
  
  preds <- dospline(dVal = dVal, my_knots = knots) %*% mcmc
  
  to_return <- list(
    x = xVal,
    y = t(
          apply(
            preds,
            1,
            quantile,
            probs = my_probs
          )
        )
  )
  
  return(to_return)
}

dmat <- data.frame(
  site = paste0("s",1:100),
  x1 = rnorm(100)
)

site_column <- "site"

dmat <- dmat[!duplicated(dmat),]

# Make the i spline design matrix. It also returns
#  a bunch of other useful features for the analysis
# Returns named list
#
#  dmat = design matrix of all variables to be
#           used for the beta diversity analysis
#
# metadata = a data.frame that contains info on
#  what covariates are in the design
#  matrix, their order, how many splines,
#  and the knots use to generate the isplines.
#  Will come in handy when making predictions.
#
# site_ids = a data.frame that has the site names
#  that are getting compared in each row of the
#  design matrix. Also has their indices, which
#  can be used in the JAGS model to fit the
#  beta diversity model.
make_spline_matrix <- function(
  dmat,
  site_column,
  x_column = NULL,
  y_column = NULL){
  # Assuming first one is geographic distance
  
  unq_sites <- length(unique(dmat[,site_column]))
  if(unq_sites != nrow(dmat)){
    stop("sites in site_column must be unique")
  }
  
  site_mat <- matrix(
    rep(1:unq_sites, each = unq_sites ),
    ncol = unq_sites,
    nrow= unq_sites
    )
  # index the upper triangle
  idx <-which(upper.tri(site_mat), arr.ind = TRUE)
  # order the matrix appropriately
  idx <- idx[order(idx[,1], idx[,2]),]
  
  # create a site pair matrix for each covariate
  # start with geographic distance

  if(all(!is.null(c(x_column, y_column)))){
    d_list <- vector(
      "list",
      length = ncol(dmat)-2
      )
    names(d_list) <- c(
      "geographic_distance",
      colnames(dmat)[-which(colnames(dmat) %in% c(
        site_column, x_column, y_column
      ))])
  } else {
    d_list <- vector(
      "list",
      length = ncol(dmat) - 1
    )
    names(d_list) <- c(
      colnames(dmat)[-which(colnames(dmat) %in% c(
        site_column, x_column, y_column
      ))])
  }
  for(i in 1:length(d_list)){
    if(names(d_list)[i] == "geographic_distance"){
      tmp <- data.frame(
        siteA_x = dmat[idx[,1], x_column],
        siteA_y = dmat[idx[,1], y_column],
        siteB_x = dmat[idx[,2], x_column],
        siteB_y = dmat[idx[,2], y_column]
      )
      # calculate distance between sites
      my_dist <- sqrt(
        (tmp$siteA_x - tmp$siteB_x)^2 + 
        (tmp$siteA_y - tmp$siteB_y)^2
      )
      d_list[[i]] <- list(
       df = tmp,
       knots = quantile(
         my_dist,
         probs = c(0,0.5,1)
       )
      )
        
    } else {
      tmp <- data.frame(
        siteA_ = dmat[idx[,1], names(d_list)[i]],
        siteB_ = dmat[idx[,2], names(d_list)[i]]
      )
      colnames(tmp) <- paste0(
        colnames(tmp), names(d_list)[i]
      )
      d_list[[i]] <- list(
        df = tmp,
        knots = quantile(
          unlist(tmp),
          probs = c(0,0.5,1)
        )
      )
    }
  }
  # Now we can construct the splines for each site pair.
  spline_list <- vector(
    "list",
    length = length(d_list)
  )
  names(spline_list) <- names(d_list)
  for(i in 1:length(d_list)){
    if("geographic_distance" == names(d_list)[i]){
      tmp <- d_list[[i]]$df
      dx <- abs(tmp$siteA_x - tmp$siteB_x)
      dy <- abs(tmp$siteA_y - tmp$siteB_y)
      dval2 <- sqrt(
        dx^2 + dy^2
      )
      # all zeroes for geographic distance
      d1a <- matrix(
        0,
        ncol = 3,
        nrow = length(dval2)
      )
      d2a <- dospline(
        dval2, d_list[[i]]$knots
      )
      spline_list[[i]] <- abs(
        d2a - d1a
      )
    } else {
      d1a <- dospline(
        d_list[[i]]$df[,1],
        d_list[[i]]$knots
      )
      d2a <- dospline(
        d_list[[i]]$df[,2],
        d_list[[i]]$knots
      )
      spline_list[[i]] <- abs(
        d2a - d1a
      )
    }
  }
  # make the matrix
  splinemat <- do.call("cbind", spline_list)
  # get knots & add number
  #  of splines
  all_knots <- matrix(
    sapply(
      d_list, "[[", 2
    ),
    nrow = 3,
    ncol = length(d_list)
  )
  # transpose to make site by knot matrix
  all_knots <- t(all_knots)
  # make it a data.frame with 
  #  extra info!
  all_knots <- data.frame(
    covariate = names(d_list),
    nspline = 3,
    min = all_knots[,1],
    median = all_knots[,2],
    max = all_knots[,3]
  )
  
  to_return <- list(
    dmat = splinemat,
    metadata = all_knots,
    dmat_site_ids = data.frame(
      siteA = unique(dmat[,site_column])[idx[,1]],
      siteA_id = idx[,1],
      siteB = unique(dmat[,site_column])[idx[,2]],
      siteB_id = idx[,2]
    )
  )
  return(to_return)
  
}

# ecological distance prediction
ed_pred <- function(
  my_splines,
  mcmc_mat,
  my_data,
  predlength = 200
){
  # steps. 
  
  # 1. Come up with the prediction based on all of the
  #  data
  
  all_preds <- cbind(1, my_splines) %*% t(mcmc_mat)
  #all_preds <- 1 - exp(-all_preds)
  
  median(all_min)
  
  # ged median estimate
  tmp <- apply(
    all_preds,
    1,
    quantile,
    probs = c(0.025,0.5,0.975)
  )
  #tmp <- 1 - exp(-tmp)
    #plogis(tmp)
  
  # 2. get min and max of distances and do something like this
  overlayX <- seq(from = min(tmp[2,]), to = max(tmp[2,]), 
                  length = predlength)
  
  
  
  
  plot(my_data ~ tmp[2,],
    #c(my_data[,1]/my_data[,2]) ~ tmp[2,],
       ylim = c(0,0.7),
       #xlim = c(0,0.7),
       ylab = "Observed dissimilarity",
       xlab = "Predicted ecological distance",
       bty = "l",
       las = 1
       )
  abline(a=0, b = 1, col = "red", lwd = 3)
  overlayY <-  1 - exp(-overlayX) # use logit instead
  lines(overlayY ~ overlayX, col = "blue", lwd = 8)
  
  # 3. Modify the above stuff to allow for 95% CI
  
  # 4. package up for plotting and return.
  
  
}

# Not from gdm, just my way to line up species info
zmat_long <- function(x, site_ids){
  if(is.list(site_ids) & "site_ids" %in% names(site_ids)){
    ids <- site_ids$site_ids
  } else if (
    !is.data.frame(site_ids) |
    !any(c("siteA_id", "siteB_id") %in% colnames(site_ids))){
    stop("site_ids must come from output of spline_dm()")
  } else {
    ids <- site_ids
  }
  to_return <- matrix(
    NA,
    ncol = 2,
    nrow = nrow(ids)
  )
  # fill this matrix.
  for(i in 1:nrow(to_return)){
    # get number of dissimilar species in site pair
    to_return[i,1] <- sum(
      (1 - x[ids$siteA_id[i],]) *
           x[ids$siteB_id[i],]
      )
    # and then the total number of unique species in site pair
    to_return[i,2] <- ncol(x) - sum(
      (1 - x[ids$siteA_id[i],]) *
      (1 - x[ids$siteB_id[i],])
    )
  }
  return(to_return)
}
# params = parameters
# X = design matrix
# y = response variable (a matrix).
loglikelihood <- function(params,X,y){
  lp <- X %*% c(params[1],exp(params[-1]))
  p <- make.link(link = "logit")
  lb <- dbinom(y[,1], y[,2], p$linkinv(lp),log = TRUE)
  ll.contr<-sum(lb)
  return(ll.contr)
}

gradient <- function(params,X,y){
  eta <- X %*% c(params[1],exp(params[-1]))
  np<-length(params)
  ns<-nrow(X)
  deri<- matrix(NA,ns,np)
  p <- make.link(link = "logit")
  for(i in 1:ns){
    mu <- p$linkinv(eta[i])
    dldm <- ((y[i,1]/y[i,2])/mu) - ((1-(y[i,1]/y[i,2]))/(1-mu))
    dmde <- p$mu.eta(eta[i])*y[i,2]
    dedb <- X[i,]
    dbdg <- c(1,exp(params[-1]))
    deri[i,] <- dldm * dmde * dedb * dbdg #chain-rule baby!
  }
  sum_deri <- apply(deri,2,sum)
  return(sum_deri)
}


fit_gdm <- function(X, y, est.var=TRUE, trace=FALSE,
                      control=bbgdm_control()){
  
  loglike <- function(x) {
    -loglikelihood(x, X, y)
  }
  grad <- function(x) {
    -gradient(x, X, y)
  }
    if(trace)control$trace <- 1
    method <- control$method
    hessian <- control$hessian
    init.par <- control$start
    fsmaxit <- control$fsmaxit
    fstol <- control$fstol
    control$method <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- control$cores <- NULL
    if(is.null(init.par)){
      fm.b <- glm.fit(X,y,family=binomial(link="logit"))
      init.par <- c(fm.b$coefficients)
    }
    if(est.var)hessian <- TRUE
    fit <- optim(par = init.par, fn = loglike, gr =grad,
                 method = method, hessian = hessian, control = control)
  invisible(fit)
  fit$par <- c(fit$par[1],exp(fit$par[-1]))
  names(fit$par) <- colnames(X)
  var <- NULL
  if (est.var) {
    var <- solve(fit$hessian)
    colnames(var) <- rownames(var) <- names(fit$par)
  }
  out <- list(coef = fit$par, logl = fit$value,
              counts=fit$counts, fitted = pi, var=var,
              X=X, y=y,control=control)
  class(out) <- "gdm"
  return(out)
}

bbgdm_control <- function (method = "BFGS", maxit = 1000, hessian = FALSE,
                           trace = FALSE, start = NULL, fsmaxit = 20, fstol = 1e-05,cores=1,
                           ...)
{
  rval <- list(method = method, maxit = maxit, hessian = hessian, trace = trace,
               start = start, fsmaxit = fsmaxit, fstol = fstol, cores=cores)
  rval <- c(rval, list(...))
  if (!is.null(rval$fnscale))
    warning("fnscale must not be modified")
  rval$fnscale <- 1
  if (is.null(rval$reltol))
    rval$reltol <- 1e-05
  rval
}

