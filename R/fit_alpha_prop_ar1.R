library(dplyr)
library(runjags)

sp_rich <- read.csv("./results/alpha_for_stage_two_collapsed.csv")

cat("Loading functions...\n")
# load functions used to clean data
functions_to_load <- list.files(
  "./R/functions/",
  full.names = TRUE
)

my_site <- read.csv("./data/cleaned_data/covariates/site_covariates.csv")


my_site$imp <- my_site$mean_19 / 100

my_site$nvuln <- !my_site$vulnerable
sp_rich <-  dplyr::inner_join(
    sp_rich,
    my_site[,c("City", "Site", "imp", "nvuln")],
    by = c("City", "Site")
  )

gent_prop <- read.csv("gent_prop1k.csv")

sp_rich <- dplyr::inner_join(
  sp_rich,
  gent_prop,
  by = c("Site", "City")
)

sp_rich <- sp_rich %>% group_by(City) %>% 
  mutate(sp = prop_gent - mean(prop_gent))

rc <-read.csv("./data/cleaned_data/covariates/racial_segregations.csv")
rc <- rc[order(rc$city),]

rc$med_scale <- scale(rc$median)
# see which cities to drop
# to_go <- within_covs %>%
#   dplyr::group_by(City) %>%
#   dplyr::summarise(
#     count = sum(gentrifying)
#   ) %>%
#   dplyr::filter( count < 4) %>%
#   data.frame()
# 
# sp_rich <- sp_rich[-which(sp_rich$City %in% to_go$City),]

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
first_sites <- dplyr::bind_rows(first_sites)
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
new_tmp <- dplyr::bind_rows(new_tmp)
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

bb <- scale(my_results$prop_gent)

data_list <- list(
  alpha_z = sp_rich$mu,
  alpha_sd_known = sp_rich$sd,
  #last_sample_vec = sp_rich$last_sample_vec,
  #alpha_first = sum(sp_rich$starter),
  #alpha_resample = sum(!sp_rich$starter),
  city_vec_alpha = as.numeric(factor(sp_rich$City)),
  npar_alpha = 4,
  #npar_among = 2,
  design_matrix_alpha = cbind(
    1, sp_rich$gentrifying,  sp_rich$imp,
    sp_rich$gentrifying * sp_rich$imp),
  #design_matrix_among = cbind(
  #  1, as.numeric(bb)
  #),
  ncity = length(unique(sp_rich$City)),
  ndata_alpha = nrow(sp_rich)#,
  #log_mu = sp_rich$log_mu,
  #log_sd = sp_rich$log_sd
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
     #alpha_mu = matrix(
     #  rnorm(data_list$npar_alpha * data_list$npar_among),
     #  nrow = data_list$npar_alpha,
     #  ncol = data_list$npar_among
     #),
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
  monitor= c("alpha", "alpha_mu", "alpha_sd", "resid_sd",
             "theta", "theta_mu", "theta_sd"),
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
which(msum[,11]>1.1)
round(summary(m1, vars = "alpha_mu"),2)


yo <- do.call("rbind", m1$mcmc)

mr <- round(msum, 2)

yo <- mr[grep(",2\\]", row.names(mr)),]


plot(yo[,2] ~ g_val$gent)
cor(yo[,2], rc$median)
plot(yo[1:data_list$ncity,2] ~ my_results$prop_gent)
abline(a = 1.4, b = -0.16)
plot(yo[,2] ~ g_val$gi)

# make a prediction




xx <- seq(0, 0.9, length.out = 400)

my_mc <- do.call("rbind", m1$mcmc)

my_mc <- split_mcmc(my_mc)




# non-gentrifying
p1 <- my_mc$alpha_mu %*% t(cbind(1,0,xx,0))
p1 <- exp(p1)
p1 <- apply(p1, 2, quantile, probs  = c(0.025,0.5,0.975))
p1 <- t(p1)

# gentrifying
p2 <- my_mc$alpha_mu %*% t(cbind(1,1,xx,1*xx))
p2 <- exp(p2)
p2 <- apply(p2 - p1, 2, quantile, probs  = c(0.025,0.5,0.975))
p2 <- t(p2)

plot()


city_gent <- city_ngent <- vector("list",
                                  length = data_list$ncity)


city_map <- data.frame(
  City = my_results$City
)
for(i in 1:data_list$ncity){
  pp <-  my_mc$alpha[,i,] %*% t(cbind(1,0,xx,0))
  pp <- exp(pp)
  city_ngent[[i]] <- apply(pp, 2, quantile, probs = c(0.025,0.5,0.975))
  
  pp <-  my_mc$alpha[,i,] %*% t(cbind(1,1,xx,xx))
  pp <- exp(pp)
  city_gent[[i]] <- apply(pp, 2, quantile, probs = c(0.025,0.5,0.975))
}



{
bbplot::blank(ylim = c(-4,4), xlim = c(0, 1), bty = "l",
              xaxs = "i", yaxs = "i")
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 1, line = 0.8)
bbplot::axis_text("Impervious cover (proportion)",1, line = 2)
bbplot::axis_text(side = 2, line = 0.8, las = 1)
bbplot::axis_text("Species richness",2, line = 2)


bbplot::ribbon(x = xx, y = p2[,-2], col = "purple", alpha = 0.3)
bbplot::ribbon(x = xx, y = p1[,-2], col = "brown", alpha = 0.3)
lines(x = xx, y = p1[,2], col = "brown", lwd = 3, lty = 2)
lines(x = xx, y = p2[,2], col = "purple", lwd = 3)
legend(
  "topright",
  c("Non-gentrifying", "Gentrifying"),
  col = c("brown","purple"),
  bty = "n",
  lwd = 4)
}

# non-gentrifying for each city

pdf("all_cities3.pdf")
for(i in 1:length(city_gent)){
  
  pp <- sp_rich[sp_rich$City == city_map$City[i],]
  pp$mulog <- convert_to_logmean(pp$mu, pp$sd)
  pp$sdlog <- convert_to_logsd(pp$mu, pp$sd)
  
  pp <- split(
    pp,
    factor(pp$Site)
    
  )
  for(j in 1:length(pp)){
    tt <- pp[[j]]
    tmptt <- vector("list", length = nrow(tt))
    for(k in 1:length(tmptt)){
      tmptt[[k]] <- rlnorm(
        10000,
        tt$mulog[k],
        tt$sdlog[k]
      )
    }
    pp[[j]] <- data.frame(
      imp = unique(pp[[j]]$imp),
      gentrifying = unique(pp[[j]]$gentrifying),
      median = median(unlist(tmptt)),
      lo = quantile(unlist(tmptt), probs = 0.025),
      hi = quantile(unlist(tmptt), probs = 0.975)
    )
  }
  pp <- do.call("rbind", pp)
  
  bbplot::blank(ylim = c(0,10), xlim = c(0, 0.9), bty = "l",
                xaxs = "i", yaxs = "i")
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.8)
  bbplot::axis_text("Impervious cover (proportion)",1, line = 2)
  bbplot::axis_text(side = 2, line = 0.8, las = 1)
  bbplot::axis_text("Species richness",2, line = 2)
  bbplot::axis_text(city_map$City[i], 3)
  
  bbplot::ribbon(x = xx, y = t(city_gent[[i]][-2,]),
                 col = "purple", alpha = 0.3)
  bbplot::ribbon(x = xx, y = t(city_ngent[[i]][-2,]),
                 col = "brown", alpha = 0.3)
  lines(x = xx, y = city_ngent[[i]][2,], col = "brown", lwd = 3, lty = 2)
  lines(x = xx, y = city_gent[[i]][2,], col = "purple", lwd = 3)
  points(
    x = pp$imp[pp$gentrifying],
    y = pp$median[pp$gentrifying],
    pch = 24, bg = "purple", cex = 1.5
  )
  points(
    x = pp$imp[!pp$gentrifying],
    y = pp$median[!pp$gentrifying],
    pch = 22, bg = "burlywood4"
  )
  legend(
    "topright",
    c("Non-gentrifying", "Gentrifying"),
    col = c("brown","purple"),
    bty = "n",
    lwd = 4)
}
dev.off()
{
  bbplot::blank(ylim = c(0,10), xlim = c(0, 1), bty = "l",
                xaxs = "i", yaxs = "i")
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.8)
  bbplot::axis_text("Impervious cover (proportion)",1, line = 2)
  bbplot::axis_text(side = 2, line = 0.8, las = 1)
  bbplot::axis_text("Species richness",2, line = 2)
  
  
  
  bbplot::ribbon(x = xx, y = p2[,-2], col = "purple", alpha = 0.3)
  bbplot::ribbon(x = xx, y = p1[,-2], col = "brown", alpha = 0.3)
  lines(x = xx, y = p1[,2], col = "brown", lwd = 3, lty = 2)
  lines(x = xx, y = p2[,2], col = "purple", lwd = 3)
  legend(
    "topright",
    c("Non-gentrifying", "Gentrifying"),
    col = c("brown","purple"),
    bty = "n",
    lwd = 4)
}

# plot out results as a function of sample size

# gent effect
my_ss <- table(sp_rich$City)

gr <- apply(
  my_mc$alpha[,,2],
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
) %>% 
  t

range(gr)
range(my_ss)

sigs <- gr[,1]> 0

{
bbplot::blank(
  ylim = c(-0.5, 1), xlim = c(50,750),
  bty = "l"
)
bbplot::axis_blank(1)
bbplot::axis_blank(2)
bbplot::axis_text(side = 1, line = 0.7)
bbplot::axis_text(side = 2, las = 1, line = 0.7)
bbplot::axis_text("Sample size", side = 1, line = 2.75)
bbplot::axis_text("Regression coefficient", side = 2, line = 3)
bbplot::axis_text("Change in species richness from gentrification")
for(i in 1:nrow(gr)){
  lines(
    x = rep(my_ss[i], 2),
    y = gr[i, -2],
    lwd = 2
  )
  points(
    x = my_ss[i],
    y = gr[i,2],
    pch = 21,
    cex = 1.25,
    bg = ifelse(
      sigs[i],
      "black", "gray"
    )
  )
}
abline(h = 0, lty = 2)
}


gr <- apply(
  my_mc$alpha[,,3],
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
) %>% 
  t

range(gr)
range(my_ss)

sigs <- gr[,3]< 0

{
  bbplot::blank(
    ylim = c(-2, 0.5), xlim = c(50,750),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.7)
  bbplot::axis_text(side = 2, las = 1, line = 0.7)
  bbplot::axis_text("Sample size", side = 1, line = 2.75)
  bbplot::axis_text("Regression coefficient", side = 2, line = 3)
  bbplot::axis_text("Change in species richness from impervious cover")
  for(i in 1:nrow(gr)){
    lines(
      x = rep(my_ss[i], 2),
      y = gr[i, -2],
      lwd = 2
    )
    points(
      x = my_ss[i],
      y = gr[i,2],
      pch = 21,
      cex = 1.25,
      bg = ifelse(
        sigs[i],
        "black", "gray"
      )
    )
  }
  abline(h = 0, lty = 2)
}



gr <- apply(
  my_mc$alpha[,,4],
  2,
  quantile,
  probs = c(0.025,0.5,0.975)
) %>% 
  t

range(gr)
range(my_ss)

sigs <- gr[,3]< 0

{
  bbplot::blank(
    ylim = c(-2.5, 1), xlim = c(50,750),
    bty = "l"
  )
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  bbplot::axis_text(side = 1, line = 0.7)
  bbplot::axis_text(side = 2, las = 1, line = 0.7)
  bbplot::axis_text("Sample size", side = 1, line = 2.75)
  bbplot::axis_text("Regression coefficient", side = 2, line = 3)
  bbplot::axis_text("Change in species richness from impervious cover x gentrification")
  for(i in 1:nrow(gr)){
    lines(
      x = rep(my_ss[i], 2),
      y = gr[i, -2],
      lwd = 2
    )
    points(
      x = my_ss[i],
      y = gr[i,2],
      pch = 21,
      cex = 1.25,
      bg = ifelse(
        sigs[i],
        "black", "gray"
      )
    )
  }
  abline(h = 0, lty = 2)
}

# non-gentrifying
betas <- array(
  NA,
  dim = c(24000, 4, 3)
)

gent_tp <- c(0.05, 0.15, 0.25)
gent_tp <- (gent_tp - mean(my_results$prop_gent)) / sd(my_results$prop_gent)

for(i in 1:4){
  betas[,i,] <- my_mc$alpha_mu[,i,] %*% rbind(1,gent_tp)
}

gent_plot <- vector("list", length = 3)

for(i in 1:3){
  p1 <- betas[,,i] %*% t(cbind(1,0,xx,0))
  #p1 <- exp(p1)
  p1 <- apply(p1, 2, quantile, probs  = c(0.025,0.5,0.975))
  p1 <- exp(t(p1))
  gent_plot[[i]] <- list(
    no_gent = p1
  )
  
  p2 <- betas[,,i] %*% t(cbind(1,1,xx,1*xx))
  #p2 <- exp(p2)
  p2 <- apply(p2, 2, quantile, probs  = c(0.025,0.5,0.975))
  p2 <- exp(t(p2))
  gent_plot[[i]]$gent <- p2
  
}

main = c(
  "Proportion of city gentrified = 0.05",
  "Proportion of city gentrified = 0.15",
  "Proportion of city gentrified = 0.25"
)
windows(12,4)
par(mfrow = c(1,3))
for(i in 1:3){
  {
    p1 <- gent_plot[[i]]$no_gent
    p2 <- gent_plot[[i]]$gent
    bbplot::blank(ylim = c(0,10), xlim = c(0, 1), bty = "l",
                  xaxs = "i", yaxs = "i")
    bbplot::axis_blank(1)
    bbplot::axis_blank(2)
    bbplot::axis_text(side = 1, line = 0.8)
    bbplot::axis_text("Impervious cover (proportion)",1, line = 2)
    bbplot::axis_text(side = 2, line = 0.8, las = 1)
    bbplot::axis_text("Species richness",2, line = 2)
    bbplot::axis_text(
      main[i],3, line = 2
    )
    
    bbplot::ribbon(x = xx, y = p2[,-2], col = "purple", alpha = 0.3)
    bbplot::ribbon(x = xx, y = p1[,-2], col = "brown", alpha = 0.3)
    lines(x = xx, y = p1[,2], col = "brown", lwd = 3, lty = 2)
    lines(x = xx, y = p2[,2], col = "purple", lwd = 3)
    legend(
      "topright",
      c("Non-gentrifying", "Gentrifying"),
      col = c("brown","purple"),
      bty = "n",
      lwd = 4)
  }
}


# get log scale range

