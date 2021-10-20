#' @title bbgdm.fit objects
#' @rdname bbgdm.fit
#' @name bbgdm.fit
#' @description function for fitting a \code{bbgdm} model. This is essential a single gdm
#' and builds the core of \code{bbgdm}. Like \code{\link[stats]{glm.fit}} a number of options
#' can be defined. Like: optimisation method \code{\link[stats]{optim}} or \code{\link[stats]{nlminb}},
#' link function \code{\link[stats]{make.link}}.
#' @param \dots for \code{bbgdm.fit()} or more \code{bbgdm.control()}.
#' @param X Model matirx of predictors, see \code{\link[stats]{model.matrix}}.
#' @param y Model response, see \code{\link[stats]{model.response}}.
#' @param wt weights for model.
#' @param link character link functions. default is 'logit', can call other
#' link functions \code{\link[stats]{make.link}} or a 'negexp' custom link function (e.g., Ferrier etal., 2007).
#' @param optim.meth optimisation method options avaliable are 'optim' and 'nlmnib',
#' calls either method \code{\link[stats]{optim}} or \code{\link[stats]{nlminb}}.
#' @param est.var logical if true estimated parameter variance using optimiser.
#' @param trace options looks at \code{\link[stats]{optim}} for details.
#' @param prior numeric vector of starting values for intercept and splines.
#' @param control control option from optim see \code{\link[bbgdm]{bbgdm.control}} or \code{\link[stats]{optim}}.
#' @return fit single gdm
#' @export
#' @author Skipton Woolley

bbgdm.fit <- function(X, y, wt=NULL, link, optim.meth="optim", est.var=TRUE, trace=FALSE,prior=FALSE,
                      control=bbgdm.control()){
  
  loglike <- function(x) {
    -loglikelihood(x, X, y, wt, link)
  }
  grad <- function(x) {
    -gradient(x, X, y, wt, link)
  }
  if(is.null(wt)){
    if(!is.null(dim(y))) wt <- rep(1,nrow(y))
    else wt <- rep(1,length(y))
  }
  if(length(wt)!=nrow(X))stop('Weights are not the same size as model matrix')
  if(optim.meth=="optim"){
    if(trace)control$trace <- 1
    method <- control$method
    hessian <- control$hessian
    init.par <- control$start
    fsmaxit <- control$fsmaxit
    fstol <- control$fstol
    control$method <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- control$cores <- NULL
    if(is.null(init.par)){
      if(link=='negexp')fm.b <- glm.fit(X,y,family=binomial(link=negexp()))
      else fm.b <- glm.fit(X,y,family=binomial(link=link))
      init.par <- c(fm.b$coefficients)
      if(prior)init.par<-rbeta(length(init.par),2,2)
    }
    if(est.var)hessian <- TRUE
    fit <- optim(par = init.par, fn = loglike, gr =grad,
                 method = method, hessian = hessian, control = control)
  }
  if (optim.meth=="nlmnib"){
    init.par <- control$start
    if(is.null(init.par)){
      if(link=='negexp')fm.b <- glm.fit(X,y,family=binomial(link=negexp()))
      else fm.b <- glm.fit(X,y,family=binomial(link=link))
      init.par <- c(fm.b$coefficients)
      if(prior)init.par<-rbeta(length(init.par),2,2)
    }
    fit <- nlminb(start=init.par,objective=loglike,gradient=grad,control = list(trace=trace))
    fit$value <- fit$objective
    fit$counts <- fit$evaluations
  }
  invisible(fit)
  fit$par <- c(fit$par[1],exp(fit$par[-1]))
  names(fit$par) <- colnames(X)
  var <- NULL
  if (est.var) {
    if(optim.meth=='optim') var <- solve(fit$hessian)
    if(optim.meth=='nlmnib') var <- solve(numDeriv::hessian(loglike,init.par))
    colnames(var) <- rownames(var) <- names(fit$par)
  }
  out <- list(coef = fit$par, logl = fit$value,
              counts=fit$counts, fitted = pi, var=var,
              X=X, y=y,control=control)
  class(out) <- "gdm"
  return(out)
}

#'@rdname bbgdm.fit
#'@name bbgdm.control
#'@param method characters string specifying the method argument passed to optim.
#'@param maxit  integer specifying the maxit argument (maximal number of iterations)
#'passed to optim.
#'@param hessian	logical. Should the numerical Hessian matrix from the optim output
#' be used for estimation of the covariance matrix? By default the analytical solution
#' is employed. For details see below.
#'@param start	an optional vector with starting values for all parameters.
#'@param fsmaxit	integer specifying maximal number of additional (quasi) Fisher scoring
#' iterations. For details see below.
#'@param fstol	numeric tolerance for convergence in (quasi) Fisher scoring.
#'For details see \code{\link[stats]{optim}}.
#'@param cores the number of cores to use in fitting.
#'@export

bbgdm.control <- function (method = "BFGS", maxit = 1000, hessian = FALSE,
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

loglikelihood <- function(params,X,y,wt,link){
  lp <- X %*% c(params[1],exp(params[-1]))
  if(link=='negexp') p <- bbgdm::negexp()
  else  p <- make.link(link = link)
  lb <- dbinom(y[,1], y[,2], p$linkinv(lp),log = TRUE)*wt
  ll.contr<-sum(lb)
  return(ll.contr)
}

gradient <- function(params,X,y,wt,link){
  eta <- X %*% c(params[1],exp(params[-1]))
  np<-length(params)
  ns<-nrow(X)
  deri<- matrix(NA,ns,np)
  if(link=='negexp')  p<- bbgdm::negexp()
  else p <- make.link(link = link)
  for(i in 1:ns){
    mu <- p$linkinv(eta[i])
    dldm <- ((y[i,1]/y[i,2])/mu) - ((1-(y[i,1]/y[i,2]))/(1-mu))
    dmde <- p$mu.eta(eta[i])*y[i,2]
    dedb <- X[i,]
    dbdg <- c(1,exp(params[-1]))
    deri[i,] <- wt[i] * dldm * dmde * dedb * dbdg #chain-rule baby!
  }
  sum_deri <- apply(deri,2,sum)
  return(sum_deri)
}

#' @rdname bbgdm.fit
#' @name negexp
#' @export

negexp<- function()
{
  linkfun <- function(mu) -log(1-mu)
  linkinv <- function(eta) 1-exp(-eta)
  mu.eta <- function(eta) exp(-eta)
  valideta <- function(eta) all(is.finite(eta))
  link <- paste0("negexp")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, name = link),
            class = "link-glm")
}


bbgdm <- function(form, sp.dat, env.dat, family="binomial",link='logit',
                  dism_metric="number_non_shared", nboot=100,
                  spline_type="ispline",spline_df=2,spline_knots=1,
                  geo=FALSE,geo.type='euclidean',coord.names=c("X","Y"),
                  optim.meth="nlmnib", est.var=FALSE, trace=FALSE,prior=FALSE,control=bbgdm.control()){
  
  cat(family,"regression is on the way. \n")
  if (family!='binomial') {
    print(family)
    stop("'family' not recognized")
  }
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  
  if(dism_metric=="number_non_shared") left <- "cbind(nonsharedspp_ij,sumspp_ij)"
  if(dism_metric=="bray_curtis") left <- "cbind(dissimilarity,100)"
  if(geo) { form <- update.formula(form, ~ X + Y + .)
  if(!all(coord.names %in% colnames(env.dat)))
  {
    stop(cat("Coordinates are missing, add coordinate data as columns called 'X' & 'Y'\n"))
  }
  }
  env.dat <- model.frame(form, as.data.frame(env.dat))
  mean.env.dat <- sd.env.dat <- NA
  env.dat <- model.frame(as.data.frame(env.dat))
  dissim_dat <- dissim_table(sp.dat,env.dat,dism_metric=dism_metric,spline_type=spline_type,spline_df=spline_df,spline_knots=spline_knots,
                             geo=geo,geo.type=geo.type)
  dissim_dat_table <- as.data.frame(dissim_dat$diff_table)
  dissim_dat_params <- dissim_dat$diff_table_params
  if(dism_metric=="number_non_shared") preds <- colnames(dissim_dat_table[,3:ncol(dissim_dat_table)])
  else preds <- colnames(dissim_dat_table[,2:ncol(dissim_dat_table)])
  form <- update.formula(form, paste0(left," ~ ",paste(preds, collapse=" + ")))
  temp <- model.frame(form, as.data.frame(dissim_dat_table))
  y <- model.response(temp)
  X <- model.matrix(form, as.data.frame(dissim_dat_table))
  mod  <- bbgdm.fit(X,y,link=link,optim.meth=optim.meth,est.var=TRUE,trace=trace,prior=prior,control=control)
  Nsite <- nrow(sp.dat)
  nreps <- nboot
  cl <- parallel::makeCluster(control$cores)
  # mods <- lapply(1:nreps,bb_apply,Nsite,X,y,link,optim.meth,est.var,trace,prior,control)
  mods <- surveillance::plapply(1:nreps, bb_apply, Nsite,X,y,link,optim.meth,est.var,trace,prior,control,.parallel = cl)
  
  #summary stats
  all.stats.ll <- plyr::ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained))
  median.ll <- apply( plyr::ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained)),2,median,na.rm=T)
  quantiles.ll <- apply( plyr::ldply(mods, function(x) c(ll=x$logl,AIC=x$AIC,BIC=x$BIC,x$null.deviance,x$gdm.deviance,x$deviance.explained)),2,function(x)quantile(x,c(.05,.95),na.rm=T))
  all.coefs.se <-  plyr::ldply(mods, function(x) c(x$coef))
  median.coefs.se <- apply(plyr::ldply(mods, function(x) c(x$coef)),2,median,na.rm=T)
  quantiles.coefs.se <- apply(plyr::ldply(mods, function(x) c(x$coef)),2,function(x)quantile(x,c(.05,.95),na.rm=T))
  
  bbgdm.results <- list()
  bbgdm.results$starting_gdm <- mod
  # bbgdm.results$bb_gdms <- mods //hopefully this will mean less memory for each model.
  bbgdm.results$all.stats.ll <- all.stats.ll
  bbgdm.results$median.ll <- median.ll
  bbgdm.results$quantiles.ll <- quantiles.ll
  bbgdm.results$all.coefs.se <- all.coefs.se
  bbgdm.results$median.coefs.se <- median.coefs.se
  bbgdm.results$quantiles.coefs.se <- quantiles.coefs.se
  bbgdm.results$nboots <-nboot
  bbgdm.results$formula <- form
  bbgdm.results$dism_metric <- dism_metric
  bbgdm.results$sp.dat <- sp.dat
  bbgdm.results$env.dat <- env.dat
  bbgdm.results$X <- X
  bbgdm.results$y <- y
  bbgdm.results$dissim_dat <- dissim_dat_table
  bbgdm.results$dissim_dat_params <- dissim_dat_params
  bbgdm.results$family <- as.character(family)[1]
  bbgdm.results$geo <- geo
  bbgdm.results$link <- link
  if(geo){
    bbgdm.results$geo.type <- geo.type
  }
  class(bbgdm.results) <- "bbgdm"
  return(bbgdm.results)
}

print.bbgdm <- function (x, ...) {
  cat(' A Bayesian Bootstrap GDM fitted against:\n',
      nrow(x$sp.dat),'sites,\n',
      ncol(x$sp.dat),'species and \n',
      nrow(x$X), 'dissimilarities used as observations in the model.\n\n')
  
  cat(' A total of',x$nboots, 'Bayesian Bootstraps were ran.\n\n')
  
  cat(' Spline base parameter estimates are: \n', paste(names(x$median.coefs.se),round(x$median.coefs.se,4), collapse="\n "),
      '\n')
}

plot.bbgdm <- function(x, ...){
  if(nrow(x$X)>200) pred_sample <- 200
  else pred_sample <- nrow(x$X)
  link <-x$link
  if(link=='negexp') link.fun <- bbgdm::negexp()
  else link.fun <- make.link(link=link)
  X <- x$X
  Y <- x$starting_gdm$y
  bb.eta <- X%*%x$median.coefs.se
  bb.eta.05 <- X%*%x$quantiles.coefs.se[1,]
  bb.eta.95 <- X%*%x$quantiles.coefs.se[2,]
  bb.pred <- link.fun$linkinv(bb.eta)
  plot(bb.eta,Y[,1]/Y[,2], type = "n",ylim=c(0,1),...)
  points(bb.eta,Y[,1]/Y[,2], pch = 20, cex = 0.25)
  y.pred <- Y[sample(pred_sample),]
  y.pred <- y.pred[order(y.pred[,2], y.pred[,1]),]
  overlayX.bb <- seq(from = min(bb.eta), to = max(bb.eta),
                     length = pred_sample)
  overlayY.bb <- link.fun$linkinv(overlayX.bb)
  lines(overlayX.bb, overlayY.bb,...)
}

bb_apply <- function(x,Nsite,X,y,link,optim.meth,est.var,trace,prior,control){
  # if(nboot>1) pb <- txtProgressBar(min = 1, max = nreps, style = 3, char = '~')
  # for (ii in 1:nreps){
  w <- gtools::rdirichlet(Nsite, rep(1/Nsite,Nsite))
  wij <- w%*%t(w)
  wij <- wij[upper.tri(wij)]
  x <- bbgdm.fit(X,y,wt=wij,link=link,optim.meth=optim.meth,est.var=est.var,trace=trace,prior=prior,control=control)
  return(x)
  # if(nboot>1) setTxtProgressBar(pb, ii)
}
