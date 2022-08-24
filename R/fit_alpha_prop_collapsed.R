library(dplyr)
library(runjags)

sp_rich <- read.csv("./results/alpha_for_stage_two_collapsed.csv")

cat("Loading functions...\n")
# load functions used to clean data
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)

for(fn in functions_to_load){
  source(fn)
}


my_site <- read.csv("./data/cleaned_data/covariates/site_covariates.csv")


my_site$imp <- my_site$mean_19 / 100


sp_rich <-  dplyr::inner_join(
    sp_rich,
    my_site[,c("City", "Site", "imp")],
    by = c("City", "Site")
  )



data_list <- list(
  alpha_z = sp_rich$mu,
  alpha_sd_known = sp_rich$sd,

  city_vec_alpha = as.numeric(factor(sp_rich$City)),
  npar_alpha = 4,
  design_matrix_alpha = cbind(
    1, sp_rich$gentrifying,  sp_rich$imp,
    sp_rich$gentrifying * sp_rich$imp),
  ncity = length(unique(sp_rich$City)),
  ndata_alpha = nrow(sp_rich)
)

inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
     alpha = matrix(
       rnorm(data_list$ncity * data_list$npar_alpha),
       ncol = data_list$npar_alpha,
       nrow = data_list$ncity
     ),
     resid = rnorm(data_list$ndata_alpha),
     theta = rnorm(data_list$ncity),
     theta_mu = rnorm(1),
     theta_tau = rgamma(1,1,1),
     alpha_mu = rnorm(data_list$npar_alpha),
     alpha_tau = rgamma(data_list$npar_alpha,1,1),
     resid_tau = rgamma(1,1,1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
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

m1 <- run.jags(
  "./JAGS/impute_alpha.R",
  monitor= c("alpha", "alpha_mu", "alpha_sd", "resid_sd"),
  n.chains = 4,
  burnin = 10000,
  sample = 20000,
  adapt = 1000,
  thin = 2,
  inits = inits,
  modules = "glm",
  method= "parallel",
  data = data_list
)

msum <- summary(m1)
range(msum[,11])

round(summary(m1, vars = "alpha_mu"),2)

saveRDS(m1, "./mcmc_output/alpha_output/alpha_mcmc.RDS")
