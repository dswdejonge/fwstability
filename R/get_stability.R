getMaxReEV <- function(JM) {
  result <- max(Re(eigen(JM)$values))
  return(result)
}

getCriticalDiagonal <- function(JM) {
  result <- diag(JM) - getMaxReEV(JM)
  return(result)
}

getStepSize <- function(criticalDiagonal, MR, iters) {
  scalars <- - criticalDiagonal / MR
  stepsize <- max(scalars, na.rm = TRUE) / iters
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

getInitialStability <- function(JM) {
  # Arnoldi, J.F., Loreau, M., Haegeman, B., 2016.
  # Resilience, reactivity and variability: A mathematical comparison
  # of ecological stability measures.
  # J. Theor. Biol. 389, 47–59.
  # https://doi.org/10.1016/j.jtbi.2015.10.012
  stability <- 0.5*getMaxReEV(JM + t(JM))
}

#' Get the mathematical stability of a matrix.
#'
#' This function can find the initial and asymptotic stability of a Jacobian matrix.
#' @details
#' Initial stability is the initial reaction of a system to a small perturbation.
#' Initial stability is found by equation 4 from Arnoldi et al. (2016),
#' with the difference that th function in this package does not change the sign of the value.
#' Therefore, a negative initial stability means the system is reactive i.e. the perturbations are amplified.
#' \cr
#' \cr
#' Asymptotic stability is the long term reaction of a system to a small perturbation.
#' Asymptotic stability can be found either as the maximum value of the real part of its eigenvalues (requires a quantified
#' diagonal) (May 1972) or as the scalar of natural mortality rates that results in a stable matrix
#' (requires mortality rate estimates) (De Ruiter et al. 1995).
#' \cr
#' \cr
#' The interpretation of the "eigenvalue" method relies on the quantification of the diagonal in the
#' Jacobian matrix. If there has been no thought on the quantification of the diagonal, the
#' "eigenvalue" method might not be informative. The diagonal can be quantified with upperbound
#' values for self-dampening by setting the argument \code{diagonal} to "model" in the function
#' \code{getJacobian}
#' \cr
#' \cr
#' The "scalar" method for asymptotic stability relies on the estimation of natural mortality rates
#' (unit per time, t-1).
#' A critical matrix, i.e. a matrix that is on the stability threshold, is calculated by setting the
#' diagonal to the given mortality values and subtracting the maximum real part of the eigenvalues
#' from each diagonal value.
#' A stepsize is then determined by estimating the scalars needed to acquire the critical matrix from
#' the given mortality values, and dividing the largest scalar by 100.
#' Subsequently, the provided natural mortality is iteratively scaled with the determined
#' stepsize until the matrix becomes stable, i.e. the maximum real part of the eigenvalues
#' is negative.
#' If dead compartments exist, their diagonal values are not scaled, because dead compartments do not
#' have natural mortality rates. Therefore, the diagonal values of dead compartments in the Jacobian matrix should
#' be quantified correctly, or set to zero (assuming there is no intracompartmental detritus feedback).
#' \cr \cr
#' Mortality rates per unit time can for example be calculated as the mortality flux (biomass
#' per surface area per unit time) divided by the biomass per surface area, or as the
#' inverse of the natural lifespan of the species. Mortality fluxes can be found with the
#' function \code{getMortalityRates} if you have the relevant data;
#' the function assumes natural mortality equals production
#' (AE x GE x Consumption) minus predation (flux to all other faunal compartments).
#' @param JM (required) A square named Jacobian matrix with numeric values representing the effect of one compartment (rows)
#' on another compartment (columns).
#' @param method (required) Either "eigenvalue" (default), "scalar", or "initial".
#' \itemize{
#' \item{
#' The method "eigenvalue" finds asymptotic stability as the maximum real part of the eigenvalues
#' calculated from the Jacobian matrix.
#' }
#' \item{
#' The "scalar" method finds asymptotic stability as the scalar of natural mortality
#' rates needed to acquire a stable matrix.
#' }
#' \item{
#' The "initial" method finds initial stability as equation 4 from Arnoldi et al. (2016);
#' initial stability is half the maximum real part of the eigenvalues derived from the
#' matrix obtained by addition of the Jacobian matrix to its transpose.
#' }
#' }
#' @param MR (required if method is "scalar")
#' A named numeric vector containing mortality of the faunal compartments
#' (per unit time, t-1). The values and names must be in the same
#' order as the Jacobian matrix, and the values for dead compartments should be
#' set to NA.
#' @param iters (required if method is "scalar") Default is 100.
#' The max number of iterations to find stability with the iterative scalar method.
#' A greater number of iters will result in a higher resolution of stability value.
#' @param dead_names (optional if method is "scalar") Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna).
#' @return This function returns a numeric value. For the "eigenvalue" method
#' a negative value indicates a stable matrix. For the "scalar" method the value represents
#' the fraction self-dampening effect needed for system stability.
#' @references \itemize{
#' \item{
#' Arnoldi, J.F., Loreau, M., Haegeman, B., 2016. Resilience, reactivity and variability: A mathematical comparison of ecological stability measures. J. Theor. Biol. 389, 47–59. https://doi.org/10.1016/j.jtbi.2015.10.012
#' }
#' \item{
#' de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of Interaction Strengths, and Stability in Real Ecosystems. Science (80-. ). 269, 1257–1260. https://doi.org/10.1126/science.269.5228.1257
#' }
#' \item{
#' May, R.M., 1972. Will a large network be stable? Nature 238, 37–38. https://doi.org/10.1038/238413a0
#' }
#' }
#' @export
#' @examples
#' \dontrun{getStability(JM)}
getStability <- function(JM, method = "eigenvalue",
                         MR = NULL, dead_names = NULL, iters = 100) {
  # Errors: check data format
  checkMformat(JM)
  checkNamingFormat(matrices = list(JM), vectors = list(MR))
  checkStabilityMethod(method, JM, MR)
  checkMortalityFormat(MR, dead_names)
  if(!is.null(dead_names) & FALSE %in% (dead_names %in% rownames(JM))) {
    stop("the names of the dead compartments are unknown")
  }
  # Warnings
  if(method == "eigenvalue" && (!is.null(MR) | !is.null(dead_names))) {
    warning("given mortality values or dead compartments are irrelevant for the eigenvalue method")
  } else if(method == "initial" && (!is.null(MR) | !is.null(dead_names))) {
    warning("given mortality values or dead compartments are irrelevant for the initial stability")
  }

  if(method == "eigenvalue") {
    stability <- getMaxReEV(JM)
  } else if (method == "scalar") {
    # Get indices compartments to scale (excluding dead compartments)
    to_scale <- 1:dim(JM)[1]
    if(!is.null(dead_names)) {
      to_scale <- to_scale[-which(rownames(JM) %in% dead_names)]
    }
    diag(JM)[to_scale] <- -MR[to_scale]
    stepsize <- getStepSize(getCriticalDiagonal(JM), MR, iters)
    stability <- getScalarStability(JM, MR, stepsize, to_scale)
  } else if(method == "initial") {
    stability <- getInitialStability(JM)
  }
  return(stability)
}
