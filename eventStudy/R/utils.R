#' @export
max2 <- function(x, na.rm){
  y <- suppressWarnings(max(x, na.rm = na.rm))
  if(y==-Inf){
    y <- 0
  }
  return(y)
}
