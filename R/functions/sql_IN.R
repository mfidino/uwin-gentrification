sql_IN <- function(x){
  to_return <- paste0("'",x,"'")
  to_return <- paste0(to_return, collapse = ", ")
  to_return <- paste0("(", to_return, ")")
  return(to_return)
}
