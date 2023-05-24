library(runjags)
library(dplyr)
library(gtools)

# There are quite a few files that got spun up when we imputed the beta
#  diversity stuff

# beta diversity and richness
beta_est <- read.csv(
  "./results/beta_summary_for_analysis_collapsed_vegan.csv"
)
beta_est <- beta_est[
  gtools::mixedorder(beta_est$loc),
]

# site ids for each site comparison
site_ids <- read.csv(
  "./mcmc_output/beta_output/dmat_site_ids_collapsed.csv"
)

# The site spline matrix. 
site_spline <- read.csv(
  "./mcmc_output/beta_output/site_splines.csv"
)

# the knots of the splines
knots <- read.csv(
  "./mcmc_output/beta_output/knots.csv"
)

# First thing to do is make the spline design matrix.
# -2 because two site columns
spline_dm <- matrix(
  NA,
  ncol = ncol(site_spline)-2,
  nrow = nrow(beta_est)
)

site_columns <- which(colnames(site_spline) %in% c("siteA","siteB"))
pb <- txtProgressBar(max = nrow(spline_dm))
for(i in 1:nrow(site_ids)){
  setTxtProgressBar(pb, i)
  # they will either be in the siteA or siteB, so 
  #  we'll need to check both and fill it in as we
  #  go
  loc1 <- which(
    site_spline$siteA == site_ids$siteA[i]   &
    site_spline$siteB == site_ids$siteB[i]
  )
  loc2 <- which(
    site_spline$siteA == site_ids$siteB[i]   &
    site_spline$siteB == site_ids$siteA[i]
  )
  if(length(loc1)>0 & length(loc2)>0){
    stop("mixed up site ids, investigate")
  }
  if(length(loc1)>0 & length(loc2)==0){
    spline_dm[i,] <- as.numeric(site_spline[loc1,-site_columns])
  }
  if(length(loc1)==0 & length(loc2)>0){
    spline_dm[i,] <- as.numeric(site_spline[loc2,-site_columns])
  }
  
}

# Going to add the gentrification sites as dummy variables
#  in their own design matrix. We are fitting those params
#  on the log scale and then will exponentiate them. So
#  we need to read in the site covariates.
site_covs <- read.csv(
  "./data/cleaned_data/covariates/site_covariates.csv"
)

# We can combine this with the site_ids to get which sites
#  were gentrifying.

gent_dm <- dplyr::inner_join(
  site_ids,
  data.frame(
    gentA = site_covs$gentrifying,
    siteA = paste0(site_covs$City,"-",site_covs$Site)
  ),
  by = "siteA"
)

# and get site B too
gent_dm <- dplyr::inner_join(
  gent_dm,
  data.frame(
    gentB = site_covs$gentrifying,
    siteB = paste0(site_covs$City,"-",site_covs$Site)
  ),
  by = "siteB"
)

# there are three comporasions. Non-gent to non-gent (control),
# non-gent vs gent, and gent vs gent.

gent_dm <- cbind(
  1,
  rowSums(gent_dm[,c("gentA","gentB")]) == 1,
  rowSums(gent_dm[,c("gentA","gentB")]) == 2
)

# add that onto the spline
spline_dm <- cbind(1,spline_dm, gent_dm[,2])

data_list <- list(
  # estimated dissimilarity 
  dissim = beta_est$mu_beta,
  # variance of dissimilarity estimate
  var_dissim = beta_est$var_beta,
  # number of data points
  ndata = nrow(beta_est),
  # log mean for richness, used in log-normal in model.
  total_rich = beta_est$mu_rich,
  # log_mu = convert_to_logmean(
  #   beta_est$mu_rich,
  #   sqrt(beta_est$var_rich)
  # ),
  # # log precision for richness, used in log-normal in model.
  # log_tau = convert_to_logsd(
  #   beta_est$mu_rich,
  #   sqrt(beta_est$var_rich)
  # )^-2,
  # gentrification design matrix
  #gent_dm = gent_dm,
  # vectors to denote the site random effects
  #siteA_vec = site_ids$siteA_id,
  #siteB_vec = site_ids$siteB_id,
  # vector to denote city id
  city_vec = site_ids$City_id,
  # design matrix for spline terms
  spline_matrix = as.matrix(spline_dm),
  # number of unique sites
  nsite = nrow(site_covs),
  # number of columns in gentrification design matrix
  #npar_gent = ncol(gent_dm),
  # number of columns in spline design matrix
  npar_spline = ncol(spline_dm),
  # number of cities
  ncity = max(site_ids$City_id)
)
saveRDS(data_list, "./mcmc_output/beta_output/data_list.RDS")
# function for initial values
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      # total_rich = abs(
      #   rlnorm(
      #     data_list$ndata,
      #     mean(data_list$log_mu),
      #     1 / mean(data_list$log_tau)
      #   )
      # ),
      #gamma = matrix(
      #  rnorm(
      #    data_list$ncity * data_list$npar_gent,
      #    -3
      #  ),
      #  nrow = data_list$ncity,
      #  ncol = data_list$npar_gent
      #),
      #site_mu = rnorm(1, -5),
      #site_tau = 1 / rgamma(1,1,1),
      #site_re = rnorm(data_list$nsite,-5),
      beta_exp = matrix(
        rgamma(
          data_list$ncity * data_list$npar_spline,
          1,1
        ),
        nrow = data_list$ncity,
        ncol = data_list$npar_spline
      ),
      #beta_a = rgamma(data_list$npar_spline,1,1),
      #beta_b = rgamma(data_list$npar_spline,1,1),
      beta_mu = abs(rnorm(data_list$npar_spline,0,0.5)),
      tau = 1 / rgamma(data_list$npar_spline,1,1),
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

mout <- run.jags(
  model = "./JAGS/beta_model_collapsed_norm.R",
  monitor = c(
    "beta_mu", "beta_exp",
    "beta_sd"#,"beta_a","beta_b"
  ),
  data = data_list,
  inits = inits,
  adapt = 1000,
  burnin = 2000, #20000,
  sample = 20000,
  n.chains = 4,
  thin = 3,
  modules = "glm",
  method = "parallel",
  jags.refresh = 60
)

saveRDS(mout, "./mcmc_output/beta_output/beta_results_collapsed_norm_vegan.RDS")

msum<- summary(mout)


round(summary(mout, vars = "beta_mu" ),2)
range(msum[,11])

which(msum[,11]>1.1)

round(exp(msum[1:8,1:3]),2)
round(tail(msum, 50),2)

my_mc <- mout$mcmc

cnames <- row.names(msum)
for( i in 1:200){
  
  my_range <- lapply(
    my_mc,
    function(x) range(x[,i])
  )
  my_range <- range(unlist(my_range))
  jpeg(
    paste0(
      "./mcmc_plots/",
      cnames[i],".jpg"
    )
  )
  plot(as.numeric(my_mc[[1]][,i]), ylim = my_range, type = "l",
       xlab = "mcmc step", ylab = "value", main = 
         paste0(cnames[i], ifelse(msum[i,11]<1.1, " converged", " NO CONVERGE")),
       col = scales::alpha("black", 0.25))
  lines(as.numeric(my_mc[[2]][,i]), col =scales::alpha("red", 0.25))
  lines(as.numeric(my_mc[[3]][,i]), col =scales::alpha("blue", 0.25))
  lines(as.numeric(my_mc[[4]][,i]), col =scales::alpha("green", 0.25))
  dev.off()
}

