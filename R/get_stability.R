getMaxReEV <- function(JM) {
  result <- max(Re(eigen(JM)$values))
  return(result)
}

getCriticalDiagonal <- function(JM) {
  result <- diag(JM) - getMaxReEV(JM)
  return(result)
}

getStepSize <- function(criticalDiagonal, MR) {
  scalars <- - criticalDiagonal / MR
  stepsize <- max(scalars, na.rm = TRUE) / 100
  return(stepsize)
}

getScalarStability <- function(JM, MR, stepsize, to_scale) {
  s <- stepsize
  # iterate to find scalar that produces stable matrix
  repeat{
    diag(JM)[to_scale] <- -MR[to_scale] * s
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
#' @details The interpretation of the "eigenvalue" method relies on the quantification of the diagonal in the
#' Jacobian matrix. If there has been no thought on the quantification of the diagonal, the
#' "eigenvalue" method might not be informative. The diagonal can be quantified with upperbound
#' values for self-dampening by setting the argument \code{diagonal} to "model" in the function
#' \code{getJacobian}
#' \cr \cr
#' The "scalar" method relies on the estimation of natural mortality rates (t-1).
#' A critical matrix, i.e. a matrix that is on the stability threshold, is calculated by setting the
#' diagonal to the given mortality values and subtracting the maximum real part of the eigenvalues
#' from each diagonal value.
#' A stepsize is then determined by estimating the scalars needed to acquire the critical matrix from
#' the given mortality values, and dividing the largest scalar by 100.
#' Subsequently, the provided natural mortality is iteratively scaled with the determined
#' stepsize until the matrix becomes stable, i.e. the maximum real part of the eigenvalues
#' is negative.
#' If dead compartments exist, their diagonal values are not scaled, because dead compartments do not
#' have mortality rates. Therefore, the diagonal values of dead compartments in the Jacobian matrix should
#' be quantified correctly, or set to zero (assuming there is no intracompartmental feedback).
#' \cr \cr
#' Mortality rates per unit time can for example be calculated as the mortality flux (biomass
#' per surface area per unit time) divided by the biomass (biomass) per surface area, or as the
#' inverse of the natural lifespan of the species. The mortality flux can be found with the
#' function \code{getMortalityRates}, which assumes natural mortality equal production
#' (AE x GE x Consumption) minus predation (flux to all other faunal compartments).
#' @param JM (required) A square named Jacobian matrix with numeric values representing the effect of one compartment (rows)
#' on another compartment (columns).
#' @param method (required) Either "eigenvalue" (default) or "scalar".
#' \itemize{
#' \item{
#' The method "eigenvalue" finds stability as the maximum real part of the eigenvalues
#' calculated from the Jacobian matrix.
#' }
#' \item{
#' The "scalar" method finds stability as the scalar of natural mortality
#' rates needed to acquire a stable matrix.
#' }
#' }
#' @param MR (required if method is "scalar")
#' A named numeric vector containing mortality of the faunal compartments
#' (per unit time, t-1). The values and names must be in the same
#' order as the Jacobian matrix, and the values for dead compartments should be
#' set to NA.
#' @param dead (optional if method is "scalar") Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna).
#' @return This function returns a numeric value. For the "eigenvalue" method
#' a negative value indicates a stable matrix. For the "scalar" method the value represents
#' the fraction self-dampening effect needed for system stability.
#' @export
#' @examples
#' getStability(JM)
getStability <- function(JM, method = "eigenvalue",
                         MR = NULL, dead = NULL) {
  # Check data format
  checkMformat(JM)
  checkNamingFormat(matrices = list(JM), vectors = list(MR))
  checkStabilityMethod(method, JM, MR)
  checkMortalityFormat(MR, dead)
  if(!is.null(dead) & FALSE %in% (dead %in% rownames(JM))) {
    stop("the names of the dead compartments are unknown")
  }
  # Warnings
  if(method == "eigenvalue" && (!is.null(MR) | !is.null(dead))) {
    warning("given mortality values or dead compartments are irrelevant for the eigenvalue method")
  }


  # Get indices compartments to scale (excluding dead compartments)
  to_scale <- 1:dim(JM)[1]
  if(!is.null(dead)) {
    to_scale <- to_scale[-which(rownames(JM) %in% dead)]
  }

  if(method == "eigenvalue") {
    stability <- getMaxReEV(JM)
  } else if (method == "scalar") {
    diag(JM)[to_scale] <- -MR[to_scale]
    stepsize <- getStepSize(getCriticalDiagonal(JM), MR)
    stability <- getScalarStability(JM, MR, stepsize, to_scale)
  }
  return(stability)
}
