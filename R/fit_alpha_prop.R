sp_rich <- read.csv("./results/alpha_for_stage_two.csv")
cat("Loading functions...\n")
# load functions used to clean data
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)

for(fn in functions_to_load){
  source(fn)
}

sp_rich$Season <- order_seasons(sp_rich$Season)

data_list <- list(
  alpha_z = sp_rich$mu,
  alpha_sd_known = sp_rich$sd,
  site_vec_alpha = as.numeric(factor(paste0(sp_rich$City,"-", sp_rich$Site))),
  city_vec_alpha = as.numeric(factor(sp_rich$City)),
  npar_alpha = 3,
  design_matrix_alpha = cbind(1, sp_rich$gentrifying, sp_rich$mean_19),
  ncity = length(unique(sp_rich$City)),
  ndata_alpha = nrow(sp_rich)
)
data_list$nsite <- max(data_list$site_vec_alpha)

m1 <- run.jags(
  "./JAGS/impute_alpha_no_among.R",
  monitor= c("alpha", "alpha_mu", "alpha_sd", "alpha_site_sd", "resid_sd"),
  n.chains = 4,
  burnin = 10000,
  sample = 20000,
  adapt = 1000,
  modules = "glm",
  method= "parallel",
  data = data_list
)

msum <- summary(m1)

round(msum, 2)
