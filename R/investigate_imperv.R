
# check out the impervious data to see what is going on

library(dplyr)
library(tidyr)


cc <- read.csv(
  "./data/imperv.csv"
)

nsite <- nrow(cc)
yrs <- c(1,4,6,8,11,16,19) - 1
yrs <- as.numeric(scale(yrs))
yrs2 <- yrs^2

my_pred <- data.frame(
  yrs = 0:18
)
dd <- c(1,4,6,8,11,16,19) - 1
my_pred$yrs <-( my_pred$yrs - mean(dd)) / sd(dd)
my_pred$yrs2 <- my_pred$yrs^2
results <- vector("list", length = nsite)
smat <- data.frame(
  site = cc$site,
  yr = FALSE,
  yr2 = FALSE
)

for(i in 1:nsite){
  dat <- as.numeric(cc[i,grep("mean", colnames(cc))])
  results[[i]] <- lm(dat ~ yrs + yrs2)
  g1 <- summary(results[[i]])$coefficients[2,4]<0.05
  g1 <- ifelse(is.na(g1), FALSE, g1)
  if(g1){
    smat$yr[i] <- TRUE
  }
  g2 <- summary(results[[i]])$coefficients[3,4]<0.05
  g2 <- ifelse(is.na(g2), FALSE, g2)
  if(g2){
    smat$yr2[i] <- TRUE
  }
}


# plot them out
lin_year <- which(smat$yr & !smat$yr2)
pdf("linyear_plot.pdf")
for(i in 1:length(lin_year)){
  tp <- predict(
    results[[lin_year[i]]],
    newdata = my_pred
  )
  plot(tp ~ c(0:18), type = "l", las = 1, bty = "l", ylim = c(0,100))
  points(results[[lin_year[i]]]$model$dat ~ dd, pch = 19)
}
dev.off()


# plot them out
c_year <- which(smat$yr2)
pdf("cyear_plot.pdf")
for(i in 1:length(c_year)){
  tp <- predict(
    results[[c_year[i]]],
    newdata = my_pred
  )
  plot(tp ~ c(0:18), type = "l", las = 1, bty = "l", ylim = c(0,100))
  points(results[[c_year[i]]]$model$dat ~ dd, pch = 19)
}
dev.off()



