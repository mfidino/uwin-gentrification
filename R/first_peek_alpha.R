# just a first little look into the alpha diversity model.
library(runjags)

source("./R/prep_data_occupancy.R")

alpha <- tmp %>% 
  dplyr::group_by(
    Site, City, Season
  ) %>% 
  dplyr::summarise(
    n = sum(Y>0) * any(J>0)
  )

alpha$Season <- factor(
  alpha$Season,
  order_seasons(alpha$Season)
)

alpha <- alpha[order(alpha$City, alpha$Season, alpha$Site),]

covars <- read.csv(
  "./data/cleaned_data/covariates/site_covariates.csv"
)

t2 <- readRDS("./data/census_data/gent_sites_proportion.rds")
t2 <- do.call("rbind", t2)

covars <- dplyr::inner_join(
  covars,
  t2,
  by = c("Site", "City")
)


covars <- covars %>% 
  dplyr::group_by(City) %>% 
  dplyr::mutate(
    gentrifying = gentrifying,
    mean_19 = (mean_19 - mean(mean_19))/100,
    prop_gent = prop_gent - mean(prop_gent)
  )

#covars$mean_19 <- covars$mean_19 / sd(covars$mean_19)
#$prop_gent <- covars$prop_gent / sd(covars$prop_gent)

alpha <- dplyr::inner_join(
  alpha,
  covars,
  by = c("City","Site")
)


test <- dplyr::inner_join(
  alpha,
  hm[,c("Site", "City","Season","mu")],
  by= c("Site", "City", "Season")
)

alpha$City <- factor(alpha$City)

alpha$site_re <- factor(
  paste0(alpha$City,"-",alpha$Site)
)


library(lme4)

m1 <- glmer(n ~ mean_19 *gentrifying + (1 + mean_19 * gentrifying|City),
            family = poisson, data = alpha)

m1 <- glmer(n ~ mean_19 +gentrifying + (1 + mean_19 + gentrifying|City),
            family = poisson, data = alpha)

data_list <- list(
  design_matrix_alpha = cbind(
    1, alpha$prop_gent, alpha$mean_19, 
    alpha$prop_gent * alpha$mean_19),
  site_vec_alpha = as.numeric(alpha$site_re),
  city_vec_alpha = as.numeric(alpha$City),
  ndata_alpha = nrow(alpha),
  npar_alpha = 4,
  ncity = nrow(city_map),
  nsite = length(levels(alpha$site_re)),
  alphaz = alpha$n
)

longshot2 <- run.jags(
  model = "./JAGS/impute_alpha.R",
  monitor = c("alpha_mu", "alpha", "alpha_sd", "alpha_site_sd"),
  data = data_list,
  adapt = 1000,
  burnin = 20000,
  sample = 20000,
  method = "parallel",
  n.chains = 3
)

msum2 <- summary(longshot2)

round(msum2,2)


# try with among city covars

my_segs <- list.files(
  "./data/cleaned_data/covariates/",
  pattern = "segregation",
  full.names = TRUE
)

seg_data <- lapply(my_segs, read.csv)
names(seg_data) <- c("education", "income", "race")

seg_data <- lapply(
  seg_data,
  function(x) x[order(x$city),]
)
seg_data <- cbind(
  seg_data$education$median,
  seg_data$income$median,
  seg_data$race$median
)

seg_data <- apply(
  seg_data,
  2, 
  scale
)

seg_data <- cbind(1, seg_data)

data_list <- list(
  design_matrix_alpha = cbind(
    1, alpha$gentrifying, alpha$mean_19,
    alpha$gentrifying * alpha$mean_19),
  site_vec_alpha = as.numeric(alpha$site_re),
  city_vec_alpha = as.numeric(alpha$City),
  ndata_alpha = nrow(alpha),
  npar_alpha = 4,
  ncity = nrow(city_map),
  nsite = length(levels(alpha$site_re)),
  alphaz = alpha$n,
  design_matrix_among = seg_data,
  npar_seg = 4
)

longshot3 <- run.jags(
  model = "./JAGS/impute_alpha2.R",
  monitor = c("alpha_mu", "alpha", "alpha_sd", "alpha_site_sd"),
  data = data_list,
  adapt = 1000,
  burnin = 20000,
  sample = 20000,
  method = "parallel",
  n.chains = 3
)


msum3 <- summary(longshot3)

round(msum3,2)


