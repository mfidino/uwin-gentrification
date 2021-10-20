# Max entropy will depend on number of categories.
# highest when equal representation. lowest when
# only 1 category.

entropy <- function(p){
  if(any(p == 0)){
    p <- p[!p==0]
  }
  sum(p * log(1/p))
}

# highest when equal split between highest and lowest
#  categories. Minimum when there is equal division
#  among all categories. Maxes out at 1.
ordinal_variation <- function(c){
  if(sum(c) != 1){
    stop("c must sum to 1")
  }
  tmp<- rep(NA, length(c)-1)
  for(i in 1:length(tmp)){
    tmp[i] <- 4 * sum(c[1:i]) * (1 - sum(c[1:i]))
  }
  (1 / (length(c)-1)) * sum(tmp)
}


w <- floor(rnorm(100, 100, 30))

tmp <- matrix(rnorm(500), ncol = 5, nrow = 100)
tmp[,1] <- 0
tmp <- t(apply(tmp, 1, function(x) exp(x) / sum(exp(x))))

my_ht <- rep(NA, 100)
for(i in 1:100){
  my_ht[i] <- entropy(tmp[i,])
}
my_ht <- apply(tmp, 1, entropy)

my_hm <- sweep(tmp, 1, w, "*")

my_hm <- colSums(my_hm) / sum(my_hm)
my_hm <- entropy(my_hm)
information_index <- function(w, hm, ht){
  sum(
    (w * (hm - ht))/(sum(w)*hm)
  )
}

information_index(w=w, my_hm, my_ht)
