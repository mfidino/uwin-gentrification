grab <- function(x, y){
  x[,grep(y, colnames(x))]
}
