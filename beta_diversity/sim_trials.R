set.seed(122)

library(mvtnorm)
library(runjags)

nsite <- 50
nspecies <- 10
nrep <- 4

psi_species <- cbind(
  -2,
  seq(-3,3, length.out = nspecies)
)

rho_species <- cbind(
  -2,
  seq(-1,1,length.out = 10)
)

x <-cbind(
  1,
  seq(-2.5,2.5, length.out = nsite)
)


# and get  occupancy
psi <- plogis(x %*% t(psi_species))

# and detection
rho <-plogis(x %*% t(rho_species))

z <- y <- psi

# fill in the z matrix
z[] <- rbinom(length(z), 1, psi)

# add observation error
y[] <- rbinom(length(z), nrep, rho * z) 


ed <- data.frame(
  site = paste0("s", 1:nsite),
  lat = runif(nsite, -47, -44),
  long = runif(nsite, 111, 117),
  x = x[,2]
)


# make spline matrix
design_matrix_beta <- make_spline_matrix(ed, "site", "long", "lat")


data_list <- list(
  # Number of sites
  nsite = nsite,
  # Number of repeat samples
  nsurvey = 4,
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
    "beta_alpha", "beta_beta",
    "sd_psi", "sd_rho", "b0"),
  data = data_list,
  n.chains = 4,
  modules = "glm",
  method = "parallel",
  inits = my_inits,
  adapt = 1000,
  burnin = 5000,
  sample = 5000,
  thin = 1
)
my_end <- Sys.time()