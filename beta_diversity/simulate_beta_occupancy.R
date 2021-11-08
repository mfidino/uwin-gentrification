set.seed(122)

library(mvtnorm)
library(runjags)

nsite <- 50
nspecies <- 30
nrep <- 4

source("./beta_diversity/gdm_functions.R")


# generate community averages
params <- list(
  # occupancy average
  bmu = c(-0.7, 0.5),
  # occupancy sd
  bsd = sqrt(c(1,2)),
  # detection average
  rmu = c(0, -0.4),
  rsd = c(0.5,0.5)
)

# Get species specific estimates
psi_species <- rmvnorm(
  nspecies,
  params$bmu,
  diag(params$bsd)
)



rho_species <- rmvnorm(
  nspecies,
  params$rmu,
  diag(params$rsd)
)

# get environmental covariate
x <-cbind(
  1,
  sort(rnorm(nsite))
)

# and get  occupancy
psi <- plogis(x %*% t(psi_species))

# and detection
rho <-plogis(x %*% t(rho_species))

# estimate species presence

# get matrices of same dimensions as psi
z <- y <- psi

# fill in the z matrix
z[] <- rbinom(length(z), 1, psi)

# add observation error
y[] <- rbinom(length(z), nrep, rho * z) 


# Calculate environmental distances
ed <- data.frame(
  site = paste0("s", 1:nsite),
  lat = runif(nsite, -47, -44),
  long = runif(nsite, 111, 117),
  x = x[,2]
)

# make spline matrix
design_matrix_beta <- make_spline_matrix(ed, "site", "long", "lat")

# get observed species

data_list <- list(
  # Number of sites
  nsite = nsite,
  nspecies = nspecies,
  z = z,
  design_matrix_psi = x,
  design_matrix_alpha = x,
  npar_psi = 2,
  npar_alpha = 2,
  ones = rep(1, nsite),
  CONSTANT = 10000
)


data_list <- list(
  # Number of sites
  nsite = nsite,
  # Number of repeat samples
  nsurvey = nrep,
  # number of unique site pairs for beta diversity model
  n = nrow(design_matrix_beta$dmat),
  # number of species
  nspecies = nspecies,
  # design matrices for all models
  design_matrix_psi = x,
  design_matrix_alpha = x,
  design_matrix_rho = x,
  design_matrix_beta = design_matrix_beta$dmat[,4:6],
  # nested indexing for beta diversity model
  siteA_id = design_matrix_beta$dmat_site_ids$siteA_id,
  siteB_id = design_matrix_beta$dmat_site_ids$siteB_id,
  # number of parameters for each model
  npar_psi = 2,
  npar_rho = 2,
  npar_alpha = 2,
  npar_beta = 3,
  # the observed data
  y = y,
  # ones trick for alpha
  alpha_ones = rep(1, nsite),
  # ones trick for beta
  beta_ones = rep(1, nrow(design_matrix_beta$dmat)),
  CONSTANT = 10000
)

# get starts for y1 and y2
z_guess <- y
z_guess[z_guess > 0] <- 1
tmp_y <- zmat_long(z_guess, design_matrix_beta$dmat_site_ids)

my_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = z_guess,
      beta_psi = matrix(
        rnorm(data_list$npar_psi * data_list$nspecies),
        nrow  = data_list$nspecies,
        ncol = data_list$npar_psi),
      beta_rho = matrix(
        rnorm(data_list$npar_psi * data_list$nspecies),
        nrow  = data_list$nspecies,
        ncol = data_list$npar_psi),
      beta_alpha = rnorm(data_list$npar_alpha),
      beta_log = rnorm(data_list$npar_beta),
      b0 = rnorm(1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

my_start <- Sys.time()
m1 <- run.jags(
  model = "./JAGS/gdm_occupancy.R",
  monitor = c(
    "beta_psi", "beta_rho",
    "beta_psi_mu", "beta_rho_mu",
    #"beta_alpha", "beta_beta",
    "sd_psi", "sd_rho", "b0","z"),
  data = data_list,
  n.chains = 5,
  modules = "glm",
  method = "parallel",
  inits = my_inits,
  adapt = 1000,
  burnin = 10000,
  sample = 5000,
  thin = 1
)
my_end <- Sys.time()

saveRDS(m1, "msom_impute.rds")
msum1 <- summary(
  m1,
  vars = c("beta_psi", "beta_rho",
           "beta_psi_mu", "beta_rho_mu",
           "sd_psi", "sd_rho")#, "beta_alpha", "beta_beta")
)
round(msum1,2)
msum2 <- summary(
  m1,
  vars = c( "beta_alpha", "beta_beta", "b0")
)


# Compare models, alpha diversity

to_plot <- list(
  bayes = NULL,
  zx = NULL,
  yx = NULL
)

# bayesian estimate
mcmc <- do.call("rbind", m1$mcmc)

mcmc <- mcmc[,grep("alpha", colnames(mcmc))]

pred_x <- cbind(1, seq(-2.75,2.75,0.05))

pred_y <- mcmc %*% t(pred_x)
pred_y <- t(exp(apply(pred_y,2, quantile, prob  = c(0.025,0.5,0.975))))
to_plot$bayes <- data.frame(
  med = pred_y[,2],
  lo = pred_y[,1],
  hi = pred_y[,3]
)


alpha_truedf <- data.frame(
  r = rowSums(z),
  x = x[,2]
)
alpha_true <- glm(r ~ x, data = alpha_truedf, family = "poisson")

summary(alpha_true)

# get predictions
true_pred <- data.frame(
  x = seq(-2.75, 2.75, 0.05)
)

tmp <- predict(
  alpha_true,
  newdata = true_pred,
  type = "link",
  se.fit = TRUE
)

# get 95% CI
to_plot$zx <- data.frame(
  med = exp(tmp$fit),
  lo = exp(qnorm(0.025, tmp$fit, tmp$se.fit)),
  hi = exp(qnorm(0.975, tmp$fit, tmp$se.fit))
)


# and with observed richness
alpha_obsdf <- data.frame(
  y = rowSums(y>0),
  x = x[,2]
)

alpha_obs <- glm(y~x, data = alpha_obsdf, family = "poisson")

tmp <- predict(
  alpha_obs,
  newdata = true_pred,
  type = "link",
  se.fit = TRUE
)
to_plot$yx <- data.frame(
  med = exp(tmp$fit),
  lo = exp(qnorm(0.025, tmp$fit, tmp$se.fit)),
  hi = exp(qnorm(0.975, tmp$fit, tmp$se.fit))
)

# plot them out

tiff("./beta_diversity/alpha_comparison.tiff",
     width = 6, height = 6, units = "in", res = 300, compression = "lzw")
plot(1~1, type = "n", ylim = c(0,25), xlim = c(-2.75,2.75),
     xlab = "Environmental gradient",
     ylab = "Species richness", las = 1, bty = "l",
     xaxs = "i", yaxs = "i")


poly_plot <- function(my_data, my_col, my_lty){
  polygon(
    x = c(seq(-2.75,2.75, 0.05), rev(seq(-2.75,2.75, 0.05))),
    y = c(my_data$lo, rev(my_data$hi)),
    col = scales::alpha(my_col, alpha = 0.15),
    border = NA
  )
  lines(y = my_data$med, x = seq(-2.75,2.75,0.05), lwd = 3,
        col = my_col, lty = my_lty)
  
}

poly_plot(to_plot$bayes, "#7fc97f",1)
poly_plot(to_plot$zx, "#beaed4",2)
poly_plot(to_plot$yx, "#fdc086",3)

legend("topleft", legend = c("Detection corrected", "Truth", "Observed richness"),
       title = "Models", col = c("#7fc97f","#beaed4","#fdc086"),
       lty = c(1,2,3), lwd = 3, bty = "n")
dev.off()

# Do the same thing but with the beta diversity data

beta_plot <- list(
  bayes = NULL,
  zx = NULL,
  yx = NULL
)


# beta diversity prediction along environmental gradient for bayes

mcmc <- do.call("rbind", m1$mcmc)


my_mcmc <- mcmc[,c("beta_beta[4]","beta_beta[5]","beta_beta[6]")] 

beta_plot$bayes <- spline_pred(
  as.numeric(design_matrix_beta$metadata[2,3:5]),
  my_mcmc
)


# beta diversity prediction along environmental gradient for zx


boot_mat <- matrix(NA, ncol = 7, nrow = 1225)
pb <- txtProgressBar(max = 1225)
for(i in 1:1225){
  setTxtProgressBar(pb, i)
  to_swap <- c(
    data_list$siteA_id[i], data_list$siteB_id[i]
  )
  tmp_z <- z
  tmp_z[to_swap[1],] <- z[to_swap[2],]
  tmp_z[to_swap[2],] <- z[to_swap[1],]
  true_beta <- zmat_long(tmp_z, design_matrix_beta$dmat_site_ids)
  dm_mod <- cbind(1, data_list$design_matrix_beta)
  boot_mat[i,] <- fit_gdm(X = dm_mod, y = true_beta)$coef
}

beta_plot$zx <- spline_pred(
  as.numeric(design_matrix_beta$metadata[2,3:5]),
  boot_mat[,5:7]
)


# and do it for the observed data

boot_mat2 <- matrix(NA, ncol = 7, nrow = 1225)
pb <- txtProgressBar(min = 1,max = 1225)
obs_z <- data_list$y
obs_z[obs_z>0] <- 1
for(i in 1:1225){
  setTxtProgressBar(pb, i)
  to_swap <- c(
    data_list$siteA_id[i], data_list$siteB_id[i]
  )
  tmp_z <- obs_z
  tmp_z[to_swap[1],] <- z[to_swap[2],]
  tmp_z[to_swap[2],] <- z[to_swap[1],]
  true_beta <- zmat_long(tmp_z, design_matrix_beta$dmat_site_ids)
  dm_mod <- cbind(1, data_list$design_matrix_beta)
  response <- try(fit_gdm(X = dm_mod, y = true_beta), silent = TRUE)
  if(class(response) != "try-error"){
    boot_mat2[i,] <- response$coef
  }
}
if(any(is.na(boot_mat2))){
  to_go <- which(is.na(boot_mat2[,1]))
  boot_mat2 <- boot_mat2[-to_go,]
}

beta_plot$yx <- spline_pred(
  as.numeric(design_matrix_beta$metadata[2,3:5]),
  boot_mat2[,5:7]
)




tiff("./beta_diversity/beta_comparison.tiff",
     width = 6, height = 6, units = "in", res = 300, compression = "lzw")
plot(1~1, type = "n", ylim = c(0,1.75), xlim = c(-2.75,2.75),
     xlab = "Environmental gradient",
     ylab = "f(Environmental gradient)", las = 1, bty = "l",
     xaxs = "i", yaxs = "i")


poly_plot <- function(my_data, my_col, my_lty){
  polygon(
    x = c(my_data$x, rev(my_data$x)),
    y = c(my_data$y[,1], rev(my_data$y[,3])),
    col = scales::alpha(my_col, alpha = 0.15),
    border = NA
  )
  lines(y = my_data$y[,2], x = my_data$x, lwd = 3,
        col = my_col, lty = my_lty)
  
}

poly_plot(beta_plot$bayes, "#7fc97f",1)
poly_plot(beta_plot$zx, "#beaed4",2)
poly_plot(beta_plot$yx, "#fdc086",3)

legend("topleft", legend = c("Detection corrected", "Truth", "Observed richness"),
       title = "Models", col = c("#7fc97f","#beaed4","#fdc086"),
       lty = c(1,2,3), lwd = 3, bty = "n")
dev.off()


hm <- fit_gdm(X = dm_mod, y = true_beta)

round(apply(boot_mat, 2, quantile, prob = c(0.025,0.5,0.975)),2)

round(msum2,2)[,1:4]

# beta diversity prediction along environmental gradient for yx




library(nnls)
library(glmnet)
# get z matrix
mcmc <- do.call("rbind", m1$mcmc)
mcmc <- mcmc[,grep("z", colnames(mcmc))]
hm <- matrix(NA, ncol = 7, nrow = nrow(mcmc))
beta_d <- matrix(NA, ncol = data_list$n, nrow = nrow(mcmc))

pb <- txtProgressBar(max = nrow(mcmc))
for(i in 1:nrow(mcmc)){
setTxtProgressBar(pb, i)
    tmp <- matrix(
    mcmc[i,],
    nrow = data_list$nsite,
    ncol = data_list$nspecies
  )
  # calcluate metrics
  tmp_y <- matrix(NA,ncol = 2, nrow =  data_list$n)
  for(j in 1:data_list$n){
    tmp_y[j,1] <- sum(
      (1 - tmp[data_list$siteA_id[j],]) *
        tmp[data_list$siteB_id[j],]
    )
    # Get total richness between site pairs
    tmp_y[j,2] <- data_list$nspecies - sum(
      (1 - z[data_list$siteA_id[j],]) *
        (1 - z[data_list$siteB_id[j],])
    )
  }
  beta_d[i,] <- tmp_y[,1] / tmp_y[,2]
  longshot <- glmnet(y = 
    cbind(tmp_y[,1], tmp_y[,2] - tmp_y[,1]),
    
    x = data_list$design_matrix_beta,
    lower.limits = c(-Inf, rep(0, 6)),
    family = "binomial"
  )
    hm[i,] <- nnls(
      cbind(1, data_list$design_matrix_beta), tmp_y[,1]/tmp_y[,2]
    )$x
}

hey <- glm(cbind(tmp_y[,1], tmp_y[,2] - tmp_y[,1])
           ~ data_list$design_matrix_beta,
           family = "binomial")

med_beta <- apply(beta_d, 2, median)

apply(hm, 2, function(x) sum(x == 0)) / nrow(hm)


set.seed(1)
n <- 100; p <- 10
x <- matrix(rnorm(n * p), nrow = n)
y <- x %*% matrix(rep(c(1, -1), length.out = p), ncol = 1) + rnorm(n)
