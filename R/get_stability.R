getMaxReEV <- function(JM) {
  result <- max(Re(eigen(JM)$values))
  return(result)
}

getCriticalDiagonal <- function(JM) {
  result <- diag(JM) - getMaxReEV(JM)
  return(result)
}

getStepSize <- function(criticalDiagonal, mortalities) {
  scalars <- - criticalDiagonal / mortalities
  stepsize <- max(scalars, na.rm = TRUE) / 100
  return(stepsize)
}

getScalarStability <- function(JM, mortalities, stepsize, to_scale) {
  s <- stepsize
  # iterate to find scalar that produces stable matrix
  repeat{
    diag(JM)[to_scale] <- -mortalities[to_scale] * s
    maxEV <- getMaxReEV(JM)
    if(maxEV < 0){break}
    s <- s + stepsize
  }
  stability <- s
  return(stability)
}

#' Get the mathematical stability of a matrix.
#'
#' This function finds the stability of a state matrix either by
#' finding the maximum value of the real part of its eigenvalues (requires a quantified
#' diagonal) or as the scalar of natural mortality rates that results in a stable matrix
#' (requires mortality rate estimates).
#' @param JM A named Jacobian matrix, with the effect of one compartment (rows)
#' on another compartment (columns). (required)
#' @param method Either "eigenvalue" (default) or "scalar". The method "eigenvalue" finds
#' stability as the maximum real part of the eigenvalues calculated from the Jacobian
#' matrix. The "eigenvalue" method relies on the quantification of the diagonal in the
#' Jacobian matrix. The "scalar" method finds stability as the scalar of natural mortality
#' rates needed to acquire a stable matrix. A critical matrix i.e. a matrix that is on the
#' stability threshold is calculated by setting the diagonal to the given mortality values
#' and subtracting the maximum real part of the eigenvalues from each diagonal value.
#' A stepsize is determined by estimating the scalars needed to acquire the critical matrix from
#' the given mortality values, and dividing the largest scalar by 100.
#' Subsequently, the provided natural mortality is iteratively scaled with the determined
#' stepsize until the matrix becomes stable, i.e. the maximum real part of the eigenvalues
#' is negative. If dead compartments exist, their diagonal values are
#' not scaled, because dead compartments do not have mortality rates. Therefore, the diagonal
#' values of dead compartments in the Jacobian matrix should be quantified correctly,
#' or set to zero (assuming there is no intracompartmental feedback). (required)
#' @param mortalities A named numeric vector containing mortality of the faunal compartments
#' (per unit time). Can be for example be calculated as mortality rate as biomass
#' per surface area per unit time divided by the biomass per surface area, or as the
#' inverse of the natural lifespan of the species. The values and names must be in the same
#' order as the Jacobian matrix, and the values for dead compartments should be
#' set to NA. (required if method is "scalar")
#' @param dead Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna). (optional)
#' @return This function returns a numeric value. For the "eigenvalue" method
#' a negative value indicates a stable matrix. For the "scalar" method the value represents
#' the fraction self-dampening effect needed for system stability.
#' @export
#' @examples
#' getStability(JM)
getStability <- function(JM, method = "eigenvalue",
                         mortalities = NULL, dead = NULL) {
  # Get indices compartments to scale (excluding dead compartments)
  to_scale <- 1:dim(JM)[1]
  if(!is.null(dead)) {
    to_scale <- to_scale[-which(rownames(JM) %in% dead)]
  }

  if(method == "eigenvalue") {
    stability <- getMaxReEV(JM)
  } else if (method == "scalar") {
    diag(JM)[to_scale] <- -mortalities[to_scale]
    stepsize <- getStepSize(getCriticalDiagonal(JM), mortalities)
    stability <- getScalarStability(JM, mortalities, stepsize, to_scale)
  }
  return(stability)
}
