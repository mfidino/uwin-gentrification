sp_rich <- read.csv("./results/alpha_for_stage_two.csv")
cat("Loading functions...\n")
# load functions used to clean data
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)



within_covs <- read.csv(
  "./data/cleaned_data/covariates/site_covariates.csv",
  stringsAsFactors = FALSE
)

# see which cities to drop
to_go <- within_covs %>%
  dplyr::group_by(City) %>%
  dplyr::summarise(
    count = sum(gentrifying)
  ) %>%
  dplyr::filter( count < 4) %>%
  data.frame()

sp_rich <- sp_rich[-which(sp_rich$City %in% to_go$City),]

for(fn in functions_to_load){
  source(fn)
}

sp_rich$Season <- factor(
  sp_rich$Season,
  order_seasons(sp_rich$Season)
)

# create autologistic term
tmp <- sp_rich

scs <- dplyr::distinct(tmp[,c("Site","City","Season")])
scs$si <- paste0(scs$City,"-",scs$Site)

scs <- split(
  scs,
  factor(scs$City)
)
tmp$Season <- factor(tmp$Season, order_seasons(tmp$Season))
all_seas <- unique(tmp$Season)
all_seas <- order_seasons(all_seas)

first_sites <- vector("list", length = length(scs))
for(i in 1:length(scs)){
  my_seas <- unique(scs[[i]]$Season)
  my_seas <- all_seas[
    which(all_seas == my_seas[1]):
      which(all_seas == my_seas[length(my_seas)])
  ]
  
  first_sites[[i]] <- data.frame(
    sites = scs[[i]]$Site[scs[[i]]$Season == my_seas[1]],
    season = my_seas[1],
    city = unique(scs[[i]]$City)
  )
  for(j in 2:length(my_seas)){
    a1 <- scs[[i]]$Site[scs[[i]]$Season == my_seas[j]]
    if(length(a1) == 0) next
    a2 <- scs[[i]]$Site[scs[[i]]$Season == my_seas[j-1]]
    tmp_si <- which(
      !a1 %in% a2
    )
    if(length(tmp_si)>0){
      first_sites[[i]] <- rbind(
        first_sites[[i]],
        data.frame(
          sites = scs[[i]]$Site[scs[[i]]$Season == my_seas[j]][tmp_si],
          season = my_seas[j],
          city = unique(scs[[i]]$City)
        )
      )
    }
  }
}
# these are the "first" instances of sampling at a location
# that may need to be linked to sit
first_sites <- do.call("rbind", first_sites)
first_sites$season <- factor(
  first_sites$season,
  order_seasons(first_sites$season)
)
first_sites <- first_sites[order(first_sites$city, first_sites$season, first_sites$sites),]
row.names(first_sites) <- NULL
# reorder data so these are all at the very beginning
my_rows <- vector("list", length = nrow(first_sites))
for(i in 1:nrow(first_sites)){
  my_rows[[i]] <-  which(
    tmp$Site == first_sites$sites[i] &
      tmp$City == first_sites$city[i]  &
      tmp$Season == first_sites$season[i]
  )
}

new_tmp <- vector("list", length = length(my_rows))
for(i in 1:length(my_rows)){
  new_tmp[[i]] <- tmp[my_rows[[i]],]
}
new_tmp <- do.call("rbind", new_tmp)
new_tmp$starter <- TRUE
tmp$starter <- FALSE
tmp <- rbind(
  new_tmp,
  tmp[-which(
    paste0(
      tmp$Site,"-",tmp$City,"-",tmp$Season
    ) %in% 
      paste0(
      new_tmp$Site,"-",new_tmp$City,"-",new_tmp$Season
      )
),]
)
row.names(tmp) <- NULL

# make a unique identifier for each sample
tmp$rowID <- 1:nrow(tmp)

# vector that incidates where the previous sample was
tmp$last_sample_vec <- NA
sp_recs <- dplyr::distinct(tmp[,c("Site","City")])
sp_recs_list <- vector("list", length = nrow(sp_recs))
tmp$Season <- factor(tmp$Season, order_seasons(tmp$Season))
# vector that indicates where the previous sample was
tmp$last_sample_vec <- NA


pb <- txtProgressBar(max = nrow(sp_recs))
for(i in 1:nrow(sp_recs)){
  setTxtProgressBar(pb, i)
  small_dat <- tmp[
      tmp$Site == sp_recs$Site[i] &
      tmp$City == sp_recs$City[i],
  ]
  small_dat <- small_dat[order(small_dat$Season),]
  small_dat$Season_id <- as.numeric(small_dat$Season)
  if(all(diff(small_dat$Season_id) == 1)){
    tmp$last_sample_vec[small_dat$rowID[-1]] <- 
      small_dat$rowID[1:(nrow(small_dat)-1)]
  } else {
    s_groups <- which(small_dat$starter)
    # check if any of the s_groups only have 0 trailing seasons
    
    for(j in 1:length(s_groups)){
      if(j < length(s_groups)){
        tdat <- small_dat[s_groups[j]:(s_groups[j+1]-1),]
      } else {
        tdat <- small_dat[s_groups[j]:nrow(small_dat),]
      }
      if(nrow(tdat) == 1) next
      if(all(diff(tdat$Season_id) == 1)){
        tmp$last_sample_vec[tdat$rowID[-1]] <- 
          tdat$rowID[1:(nrow(tdat)-1)]
      } else {
        stop("investigate")
      }
    }
  }
}

# check to see if it is correct
my_seas <- levels(tmp$Season)
for(i in rev(1:nrow(tmp))){
  if(tmp$starter[i]) next
  my_eval <- 
    tmp$Site[tmp$last_sample_vec[i]] == tmp$Site[i] &
    tmp$Season[tmp$last_sample_vec[i]] ==
    my_seas[which(my_seas == tmp$Season[i])-1]
  if(!my_eval){
    stop("last_sample_vec wrong, investigate.")
  }
}

sp_rich <- tmp

# convert mean and sd to log mean and sd
convert_to_logmean <- function(mean, sd) {
  log(mean^2 / sqrt(sd^2 + mean^2))
}

convert_to_logsd <- function(mean, sd) {
  sqrt(log(1 + (sd^2 / mean^2)))
}

sp_rich$log_mu <- convert_to_logmean(
  sp_rich$mu,
  sp_rich$sd
)

sp_rich$log_sd <- convert_to_logsd(
  sp_rich$mu,
  sp_rich$sd
)


# get mean impervious cover
mean_imp <- within_covs %>% 
  dplyr::group_by(
    City
  ) %>% 
  dplyr::summarise(
    mean_imp = mean(mean_19)
  ) %>% 
  data.frame()
# scale it
mean_imp <- mean_imp[mean_imp$City %in% 
                       sp_rich$City,]
mean_imp$mean_imp <- scale(mean_imp$mean_imp)


data_list <- list(
  alpha_z = sp_rich$mu,
  alpha_sd_known = sp_rich$sd,
  last_sample_vec = sp_rich$last_sample_vec,
  alpha_first = sum(sp_rich$starter),
  alpha_resample = sum(!sp_rich$starter),
  city_vec_alpha = as.numeric(factor(sp_rich$City)),
  npar_alpha = 4,
  npar_among = 2,
  design_matrix_alpha = cbind(
    1, sp_rich$gentrifying, sp_rich$mean_19,
    sp_rich$gentrifying * sp_rich$mean_19),
  design_matrix_among = cbind(
    1, mean_imp$mean_imp
  ),
  ncity = length(unique(sp_rich$City)),
  ndata_alpha = nrow(sp_rich),
  log_mu = sp_rich$log_mu,
  log_sd = sp_rich$log_sd,
  nsite = length(
    unique(
      paste0(sp_rich$City, "-", sp_rich$Site)
    )
  )
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
     # alpha_mu = matrix(
     #   rnorm(data_list$npar_alpha * data_list$npar_among),
     #   nrow = data_list$npar_alpha,
     #   ncol = data_list$npar_among
     # ),
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
  "./JAGS/impute_alpha_ar1.R",
  monitor= c("alpha", "alpha_mu", "alpha_sd", "resid_sd",
             "theta", "theta_mu", "theta_sd"),
  n.chains = 3,
  burnin = 2500,
  sample = 5000,
  adapt = 1000,
  thin = 1,
  inits = inits,
  modules = "glm",
  method= "parallel",
  data = data_list
)

msum <- summary(m1)

summary(m1, vars = "alpha_mu")


yo <- do.call("rbind", m1$mcmc)

mr <- round(msum, 2)

yo <- mr[grep(",2\\]", row.names(mr)),]


g_val <- within_covs %>% 
  dplyr::group_by(City) %>% 
  dplyr::summarise(gent = mean(gentrifying),
                   imp = mean(mean_19),
                   gi = mean(gentrifying * mean_19),
                   count = sum(gentrifying))

plot(yo[,2] ~ g_val$gent)
plot(yo[1:16,2] ~ data_list$design_matrix_among[,2])
plot(yo[,2] ~ g_val$gi)



