library(runjags)

# bring in the chains, and the data
sp_rich <- read.csv(
  "./results/alpha_for_stage_two_collapsed.csv"
)

# load functions for splitting mcmc
source(
  "./R/alpha_beta_functions.R"
)

# covariates
my_site <- read.csv(
  "./data/cleaned_data/covariates/site_covariates.csv"
)

# divide by 100
my_site$imp <- my_site$mean_19 / 100

# join to data
sp_rich <-  dplyr::inner_join(
  sp_rich,
  my_site[,c("City", "Site", "imp")],
  by = c("City", "Site")
)
# for plotting
city_map <- data.frame(
  city = unique(sp_rich$City)
)
# pretty name
city_map$pname <- c(
  "Athens, GA",
  "Bay Area, CA",
  "Boston, MA",
  "Chicago, IL",
  "Denver, CO",
  "Houston, TX",
  "Indianapolis, IN",
  "Iowa City, IA",
  "Jackson, MS",
  "Little Rock, AR",
  "Madison, WI",
  "Metro LA, CA",
  "National Capital",
  "Phoenix, AZ",
  "Portland, OR",
  "Rochester, NY",
  "Sanford, FL",
  "Salt Lake City, UT",
  "Seattle, WA",
  "St. Louis, MO",
  "Tacoma, WA",
  "Urbana, IL",
  "Wilmington, DE"
)



my_mc <- readRDS(
  "./mcmc_output/alpha_output/alpha_mcmc.RDS"
)

my_mc <- do.call(
  "rbind",
  my_mc$mcmc
)

my_mc <- split_mcmc(
  my_mc
)

xx <- seq(0, 0.8, length.out = 300)

# get the average across the study, but set it up as the 
#  difference
# non-gentrifying
p1 <- my_mc$alpha_mu %*% t(cbind(1,0,xx,0)) 
p1 <- exp(p1)

# gentrifying
p2 <- my_mc$alpha_mu %*% t(cbind(1,1,xx,xx))
p2 <- exp(p2)
# difference
p2 <- apply(p2 - p1, 2, quantile, probs  = c(0.025,0.5,0.975))
mu_diff <- t(p2)



# predictions for each city
city_gent <- city_ngent <- vector(
  "list",
  length = nrow(city_map)
)
names(city_gent) <- city_map$city

pb <- txtProgressBar(max = nrow(city_map))
for(i in 1:nrow(city_map)){
  setTxtProgressBar(pb, i)
  prop_gent <- prop.table(
    table(
      sp_rich$gentrifying[sp_rich$City == city_map$city[i]]
    ))
  pp <-  my_mc$alpha[,i,] %*% t(cbind(1,0,xx,0))
  pp <- exp(pp)
  
  pp2 <-  my_mc$alpha[,i,] %*% t(cbind(1,1,xx,xx))
  pp2 <- exp(pp2)
  city_gent[[i]] <- t(
    apply(
      (pp2 * prop_gent[1]) + (pp * prop_gent[2]),
      2,
      quantile,
      probs = c(0.025,0.5,0.975)
    )
  )
}

# get which cities we found a difference
gent_1 <- t(
  apply(
    my_mc$alpha[,,2],
    2,
    quantile, 
    probs = c(0.025,0.975)
  )
)

# sig effect
gent_1 <- apply(
  sign(gent_1),
  1,
  function(x) length(unique(x)) == 1
)
sum(gent_1)

housing <- t(
  apply(
    my_mc$alpha[,,3],
    2,
    quantile, 
    probs = c(0.025,0.975)
  )
)

housing <- apply(
  sign(housing),
  1,
  function(x) length(unique(x)) == 1
)
sum(
  housing
)

gent_2 <- t(
  apply(
    my_mc$alpha[,,4],
    2,
    quantile, 
    probs = c(0.025,0.975)
  )
)

# sig effect
gent_2 <- apply(
  sign(gent_2),
  1,
  function(x) length(unique(x)) == 1
)
sum(gent_2)

# getting some colors together
c1 <- col2rgb("chocolate1")
c2 <- col2rgb("cornflowerblue")
nullc <- col2rgb("white")

combo <- rgb(t((c1 + c2) / 2), maxColorValue = 255)

tmp1 <- (((c1 +nullc) / 2 ) + ((c2 + nullc) /2)) / 2

my_cols <- c(
  rgb(t(c2), maxColorValue = 255),
  combo,
  rgb(t(tmp1), maxColorValue = 255),
  rgb(t(c1), maxColorValue = 255)
  
)
names(my_cols) <- c("just_imp", "both", "none", "just_gen")


windows(12,8)
# try out a bivariate legend
svg("./plots/alpha_impervious.svg",
     width = 12, height = 8)
{

m <- matrix(
  c(1:5,24,6:23),
  ncol = 6,
  nrow = 4,
  byrow = TRUE
)
layout(m)
par(mar = c(1,1,1,1), oma = c(6,6,0,0))


# bbplot::blank(ylim = c(0,1.2), xlim = c(0,1), xaxs = "i", yaxs = "i")
# rect(0.4,0.4,0.7,0.7, col = my_cols[3])
# rect(0.7,0.7,1,1, col = my_cols[2])
# rect(0.4,0.7,0.7,1, col = my_cols[4])
# rect(0.7,0.4,1,0.7, col = my_cols[1])
# text(x = 0.55, y = 0.3,  "N", adj = 0.5, cex = 2.25)
# text(x = 0.85, y = 0.3,  "Y", adj = 0.5, cex = 2.25)
# text(x = 0.3, y = 0.55, "N", adj = 0.5, cex = 2.25)
# text(x = 0.3, y = 0.85, "Y", adj = 0.5, cex = 2.25)
# text("Impervious", x = 0.67, y = 0.12, ad = 0.5, cex = 1.67)
# text("Gentrification", x = 0.1, y = 0.62, ad = 0.5, cex = 1.67, srt = 90)

mean_response <- apply(
  my_mc$alpha[,,3],
  2, 
  median
)

biggest_diff <- rep(NA, 23)
for(i in 1:23){
  biggest_diff[i] <- city_gent[[i]][1,2] - city_gent[[i]][300,2]
}
my_order <- order(biggest_diff, decreasing = TRUE)
for(i in 1:24){
  
  if(i == 24){
   bbplot::blank()
    par(xpd = NA)
    legend(
      x = 0.57,
      y = 1.6,
      legend = c(
        "Median estimate",
        "Imp. cover 95% CI\nexcludes zero",
        "Imp. cover 95% CI\nincludes zero",
        "Site estimate",
        "Site 95% CI"
      ),
      lty = c(1,NA,NA,NA,NA,NA),
      pch = c(NA, 22,22,21,NA),
      pt.cex = c(NA, 4,4,2,NA),
      pt.bg = c(
        NA,
        my_cols["just_imp"],
        my_cols["none"],
        scales::alpha("black", 0.25),NA),
      lwd = c(3,NA,NA,NA,NA),
      bty = "n",
      y.intersp = 1.3,
      cex = 1.2
    )
    tmp <- 0.64
    lines(x = rep(0.65,2), y = c(tmp, tmp+0.15),
          col = scales::alpha("black", 0.7),
          lwd = 2)
    bbplot::axis_text(
      "Impervious cover (proportion)",
      1,
      outer = TRUE,
      line = 3.25,
      at = NA,
      cex = 1.5
    )
    bbplot::axis_text(
      "Species richness",
      2,
      outer = TRUE,
      line = 3.25,
      at = NA,
      cex = 1.5
    )
    next
  }
  my_city <- city_map$city[my_order[i]]


  
  if(housing[my_order[i]]){
    col_choice <- my_cols["just_imp"]
  }

  if(!housing[my_order[i]]){
    col_choice <- my_cols["none"]
  }
  
  bbplot::blank(ylim = c(0,22), xlim = c(0, 0.8), bty = "l",
                main = city_map$pname[my_order[i]])
  bbplot::axis_blank(1)
  bbplot::axis_blank(2)
  if(i %in% m[,1]){
    bbplot::axis_text(side = 2, las = 1, line = 0.5)
  }
  if(i %in% m[4,]){
    bbplot::axis_text(side = 1, line = 0.75)
  }
  
  my_dots <- sp_rich[sp_rich$City == city_map$city[my_order[i]],c("mu","sd","imp")]
  
  my_dots$lower <- qnorm(
    c(0.025),
    my_dots$mu,
    my_dots$sd
  )
  my_dots$lower[my_dots$lower<0] <- 0
  my_dots$upper<- qnorm(
    c(0.975),
    my_dots$mu,
    my_dots$sd
  )
  my_dots$upper[my_dots$upper>22] <- 22

  
  bbplot::ribbon(x = xx, y = city_gent[[my_city]][,-2], alpha = 1, col = col_choice)
  

  lines(x = xx, y = city_gent[[my_city]][,2], lwd = 4)
  #lines(x = xx, y = city_gent[[i]][,2], lwd = 1, col = "white")
  for(j in 1:nrow(my_dots)){
    lines(x = rep(my_dots$imp[j],2),
         y = as.vector(my_dots[j,c("lower","upper")]),
         col = scales::alpha("black", 0.25))
    points(x = my_dots$imp[j], y = my_dots$mu[j], pch = 19,
           col = scales::alpha("black", 0.25))
  }

}




}
dev.off()

# non-gentrifying for each city

  
  for(i in 1:length(city_gent)){
    lines(
      x = xx, y = city_gent[[i]][,2],
      lwd = 2, col = "gray50"
    )
  }
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

