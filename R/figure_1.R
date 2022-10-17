library(pals)
library(bbplot)
library(vegan)



m1 <- matrix(
  c(0,0,1,0,1,
    1,1,0,1,1),
  ncol = 5,
  nrow = 2, byrow = TRUE
)
cmats <- list(
  plot1  = matrix(
    c(0,0,1,0,1,
      1,1,0,1,1),
    ncol = 5,
    nrow = 2, 
    byrow = TRUE
  ),
  plot2 = matrix(
    c(0,1,1,0,0,
      0,0,0,1,1),
    ncol = 5,
    nrow = 2,byrow = TRUE
  ),
  plot3 = matrix(
    c(0,0,1,1,1,
      1,1,1,1,1),
    ncol = 5,
    nrow = 2,
    byrow = TRUE
  ),
  plot4 = matrix(
    c(1,0,1,0,1,
      1,0,1,0,1),
    ncol = 5,
    nrow = 2,
    byrow = TRUE
  )
)

sp <- rev(21:25)
my_cols <- pals::brewer.seqseq2(9)[c(7,3,8,6,9)]


ab <- function(x){
  alpha <- diff(rowSums(x))
  
  beta <- vegan::vegdist(
    x,
    binary =TRUE
  )
  return(c('alpha' = alpha, 'beta' = beta))
}

ab_vals <- lapply(
  cmats,
  ab
)

sp_loc <- function(x){
  ind_rich <- as.character(rowSums(x))
  my_res <- vector("list", 2)
  for(i in 1:2){
  my_res[[i]] <- switch(ind_rich[i],
         "1" = data.frame(x = 5, y = 5),
         "2" = data.frame(x = c(-1,1), y = c(5,5)),
         "3" = data.frame(x = c(-1,0,1), y = c(5,7.5,5)),
         "4" = data.frame(x = c(-1,0,1, 0), y = c(6.5,9,6.5,4)),
         "5" = data.frame(
           x = c(-1,0,1, 0, 0), y = c(6.5,9,6.5,4, 6.5))
  )
  }
  return(my_res)
}
gplot <- function(){
  bbplot::blank(xlim = c(0,10), xaxs = "i", ylim = c(0,10), yaxs = "i")
  lines(x = c(0.25, 4.75), y = c(2.25,2.25), lwd = 2)
  lines(x = c(5.25, 9.75), y = c(2.25,2.25), lwd = 2)
}

windows(6,4)
m <- matrix(
  c(1,1,5,6,
    2,2,7,8,
    3,3,9,10,
    4,4,11,12),
  ncol = 4,
  nrow = 4,
  byrow = TRUE
)
tiff("tester.tiff", height = 4, width = 6, units = "in",
     res = 1200, compression = "lzw")
layout(m)
par(
  mar = c(0.5,0.5,0.5,0.5),
  xpd = NA
)

for(i in 1:4){
  gplot()
  u <- par("usr")
  text(x = 0.25, y = 9.75, paste0(LETTERS[i],")"))
  g1 <- sp_loc(cmats[[i]])
  g1[[1]][,1] <- g1[[1]][,1] + 2.5 
  g1[[2]][,1] <- g1[[2]][,1] + 7.5
  points(
    g1[[1]], 
    pch = sp[cmats[[i]][1,] == 1],
    bg = "black",
    cex = 3.5
  )
  points(
    g1[[1]], 
    pch = sp[cmats[[i]][1,] == 1],
    bg = my_cols[cmats[[i]][1,] == 1],
    cex = 3
  )
  points(
    g1[[2]], 
    pch = sp[cmats[[i]][2,] == 1],
    bg = "black",
    cex = 3.5
  )
  points(
    g1[[2]], 
    pch = sp[cmats[[i]][2,] == 1],
    bg = my_cols[cmats[[i]][2,] == 1],
    cex = 3
  )
  if(i ==4){
    text(2.5, y = 2, "Not gentrified", pos = 1, cex = 1.5)
    text(7.5, y = 2, "Gentrified", pos = 1, cex = 1.5)
  }
}
dev.off()
tp1 <- ab()



#centre points
centre_x = 5
centre_y = 5
#radius
r = 2.25

deg2rad <- function(d){
  return((d*pi)/180)
} #Converts Degrees to radians
X_coord <- function(r=2.25,centre_x=2.5,angle) #Finds Xcoordinate on the circumference 
{
  return(r*cos(deg2rad(angle)) + centre_x)
}
Y_coord <- function(r=2.25,centre_y=2.5,angle) #Finds Ycoordinate on the circumference 
{
  return(r*sin(deg2rad(angle)) + centre_x)
}

# series of angles after dividing the circle in to 5 
angles <- list()
for(i in 1:5)
{
  angles[i] <- 72*i
}
angles <- unlist(angles) #flattening the list 

for(i in seq_along(angles)){
  print(i)
  print(angles[i])
  if(i == 1)
  {
    coordinates <- 
      cbind(c(
        x = X_coord(angle = angles[i]),
        y = Y_coord(angle = angles[i]))
      )
  }
  else{
    coordinates <- cbind(coordinates,cbind(c(
      x = X_coord(angle = angles[i]),
      y = Y_coord(angle = angles[i]))))
  }
}
plot(xlim = c(0,30), ylim = c(0,30),x = coordinates[1,], y=coordinates[2,], asp  =1)

polygon(x = coordinates[1,c(1,3,5,2,4,1)],                           
        y=coordinates[2,c(1,3,5,2,4,1)],                             
        col = "#1b98e0",                                             
        border = "red",                                              
        lwd = 5)
