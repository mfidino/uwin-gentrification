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


sp_rich <-  dplyr::inner_join(
    sp_rich,
    my_site[,c("City", "Site", "imp")],
    by = c("City", "Site")
  )


for(fn in functions_to_load){
  source(fn)
}

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
  n.chains = 2,
  burnin = 10000,
  sample = 40000,
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

