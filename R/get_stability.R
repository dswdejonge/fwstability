getStability <- function(JM) {
  stability <- max(Re(eigen(JM)$values))
  return(stability)
}
