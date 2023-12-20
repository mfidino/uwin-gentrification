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
  plot1 = matrix(
    c(0,1,1,0,0,
      0,0,0,1,1),
    ncol = 5,
    nrow = 2,byrow = TRUE
  ),
  plot2  = matrix(
    c(0,0,1,0,1,
      1,1,0,1,1),
    ncol = 5,
    nrow = 2, 
    byrow = TRUE
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

sp <- 21:25 
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
    # anchoring in place now by plotting all
    #  points, then only coloring some in.
    my_res[[i]] <- data.frame(
      x = c(-1,0,1,-0.5, 0.5),
      y = c(4,4,4,6.65,7.35)
    )
  }
  return(my_res)
}
gplot <- function(){
  bbplot::blank(xlim = c(0,10), xaxs = "i", ylim = c(0,10), yaxs = "i")
  lines(x = c(0.25, 4.75), y = c(2.25,2.25), lwd = 2)
  lines(x = c(5.25, 9.75), y = c(2.25,2.25), lwd = 2)
}

#windows(6,4)
m <- matrix(
  c(
    1,5,5,
    2,5,5,
    3,5,5,
    4,5,5
  ),
  ncol = 3,
  nrow = 4,
  byrow = TRUE
)

tiff("./plots/fig_1_gentrification.tiff", height = 4, width = 6, units = "in",
     res = 1200, compression = "lzw")
{
  layout(m)
  
  par(
    mar = c(0.5,0.5,0.5,0.5),
    oma = c(0,0,4,0),
    xpd = NA
  )
  
  plot_titles <- paste0(
    LETTERS[1:4],")", c(" Distinct", " Increased", " Nested", " No difference")
  )

  for(i in 1:4){
    gplot()
    u <- par("usr")
    text(x = 0.15, y = 9.75, plot_titles[i], pos = 4, cex = 1.4)
    g1 <- sp_loc(cmats[[i]])
    g1[[1]][,1] <- g1[[1]][,1] + 2.5 
    g1[[2]][,1] <- g1[[2]][,1] + 7.5
    # left side
    points(
      g1[[1]], 
      pch = sp,#sp[cmats[[i]][1,] == 1],
      bg = ifelse(
        cmats[[i]][1,],
        "black",
        "gray70"
      ),
      col = ifelse(
        cmats[[i]][1,],
        "black",
        "gray70"
      ),
      cex = 3,
      lwd = 1.1
    )
    tmp_cols <- my_cols
    tmp_cols[cmats[[i]][1,] == 0] <- "white"
    points(
      g1[[1]], 
      pch = sp,#sp[cmats[[i]][1,] == 1],
      bg = tmp_cols,
      col = ifelse(
        cmats[[i]][1,],
        "black",
        "gray70"
      ),
      cex = 2.5
    )
    # right side
    points(
      g1[[2]], 
      pch = sp,#sp[cmats[[i]][1,] == 1],
      bg = ifelse(
        cmats[[i]][2,],
        "black",
        "gray70"
      ),
      col = ifelse(
        cmats[[i]][2,],
        "black",
        "gray70"
      ),
      cex = 3,
      lwd = 1.1
    )
    tmp_cols <- my_cols
    tmp_cols[cmats[[i]][2,] == 0] <- "white"
    points(
      g1[[2]], 
      pch = sp,#sp[cmats[[i]][1,] == 1],
      bg = tmp_cols,
      col = ifelse(
        cmats[[i]][2,],
        "black",
        "gray70"
      ),
      cex = 2.5
    )
   
    if(i ==1){
      text(2.5, y = 18.3-1.3, "Site", pos = 1, cex = 1.4)
      text(7.5, y = 18.3-1.3, "Site", pos = 1, cex = 1.4)
      text(2.5, y = 16-1.3, "not gentrified", pos = 1, cex = 1.4)
      text(7.5, y = 16-1.3, "gentrified", pos = 1, cex = 1.4)
    }

  }
  
  par(mar = c(3.5,5.5,0,4))
  
  bbplot::blank(xlim = c(-0.1, 2), ylim = c(-0.1, 1.2), bty = 'n',
                xaxs = "i", yaxs= "i")
  u <- par("usr")
  bbplot::axis_text(
    "E) Difference in alpha and beta diversity",
    line = 0.58,
    cex = 0.88,
    at = 0.95
    )
  box(which = "plot", bty = "l", lwd = 2)
  
  bbplot::axis_blank(
    side = 1,
    minor = FALSE,
    lwd = 2
  )
  
  bbplot::axis_text(
    text = c(0,2),
    side = 1,
    at = c(0,2),
    line = 0.75
  )
  bbplot::axis_text(
    text = c(0,1),
    side = 2,
    at = c(0,1),
    line = 0.75,
    las = 1
  )
  bbplot::axis_blank(
    at = c(0,0.5,1),
    side = 2,
    lwd = 2,
    minor = FALSE,
    tck = -0.03
  )
  
  bbplot::axis_text(
    text = bquote(Delta ~"alpha"),
    side = 1,
    line = 2.25,
    cex = 1.3
  )
  
  bbplot::axis_text(
    text ="beta",
    side = 2,
    line = 2.25,
    cex = 1.3
  )
  
  
  plot_titles2 <- c("Distinct", "Increased", "Nested", "No difference")
  for(i in 1:4){
    points(
      x = ab_vals[[i]][1], #- 0.01,
      y = ab_vals[[i]][2] - 0.0025,
      pch = 21,
      col = "gray",
      cex = 4
      
    )
    text(
      LETTERS[i],
      x = ab_vals[[i]][1],
      y = ab_vals[[i]][2],
      adj = 0.5,
      cex = 1.7
    )
    text(
      plot_titles2[i],
      x = ab_vals[[i]][1] + c(0.07,0 , 0, 0.07)[i],
      y = ab_vals[[i]][2] + c(-0.007,0.03,0.03, -0.007)[i],
      pos = c(4,3,3,4)[i],
      cex = 1.7
    )
  }
  
}



dev.off()

