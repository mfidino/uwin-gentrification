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

sp_rich$Season <- factor(
  sp_rich$Season,
  order_seasons(sp_rich$Season)
)

test <- dplyr::inner_join(
  sp_rich,
  covars,
  by = c("City", "Site")
)

data_list <- list(
  alpha_z = sp_rich$mu,
  alpha_sd_known = sp_rich$sd,
  site_vec_alpha = as.numeric(factor(paste0(sp_rich$City,"-", sp_rich$Site))),
  city_vec_alpha = as.numeric(factor(sp_rich$City)),
  npar_alpha = 4,
  design_matrix_alpha = cbind(
    1, sp_rich$gentrifying, sp_rich$mean_19,
    sp_rich$gentrifying * sp_rich$mean_19),
  ncity = length(unique(sp_rich$City)),
  ndata_alpha = nrow(sp_rich)
)
data_list$nsite <- max(data_list$site_vec_alpha)

mm <- covars %>% 
  group_by(City) %>% 
  summarise(gen = mean(gentrifying),
            imp = mean(mean_19))


m1 <- run.jags(
  "./JAGS/impute_alpha_no_among.R",
  monitor= c("alpha", "alpha_mu", "alpha_sd", "alpha_site_sd", "resid_sd"),
  n.chains = 3,
  burnin = 7500,
  sample = 10000,
  adapt = 1000,
  thin = 10,
  modules = "glm",
  method= "parallel",
  data = data_list
)

msum <- summary(m1)



yo <- do.call("rbind", m1$mcmc)

mr <- round(msum, 2)

mr[grep(",4\\]", row.names(mr)),]

