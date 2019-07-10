#' Get the mathematical stability of a matrix.
#'
#' This function allows you to find the stability of a state matrix by
#' finding the maximum value of the real part of its eigenvalues.
#' @param JM A Jacobian matrix.
#' @return a numeric value, a negative value indicates a stable matrix.
#' @export
#' @examples
#' getStability(JM)
getStability <- function(JM) {
  stability <- max(Re(eigen(JM)$values))
  return(stability)
}
