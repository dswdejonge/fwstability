#' Calculate diagonal value for species
#'
#' This function finds the diagonal values for species to be used in the Jacobian matrix.
#' The calculation is based on their non-predatory mortality rate and biomass, and is
#' an implementation of the equation from Neutel & Thorne (2014). \cr
#' References: \cr
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' @param MR A named numeric vector with non-predatory mortality rates for all
#' compartments (biomass per unit time or biomass per unit time per surface area).
#' (required)
#' @param BM A named numeric vector with biomasses of all compartments
#' (just biomass or biomass per unit area), must be in the same order as MR. (required)
#' @return This function returns a named numeric vector with diagonal values for the
#' species in the food web (per unit time). It is important to review the units of the
#' input data. If MR is biomass per unit time then BM must be just biomass. If MR is
#' biomass per unit time per surface area then BM must be biomass per surface area.
#' @export
#' @examples
getDiagonalSpecies <- function(MR, BM) {
  if(length(MR) != length(BM)) {
    stop("input vectors have unequal lengths")
  } else if(!is.numeric(MR) | !is.numeric(BM)) {
    stop("input vectors must be numeric")
  }
  result <- -MR/BM
  #result <- -getMortalityRate(flow_solutions, BM, dead)
  #result <- result[-dead]
  return(result)
}

#' Calculate diagonal for dead (detritus) compartments
#'
#' This function finds diagonal values for detritus compartments to be used in the Jacobian
#' matrix. The calculation is based on the total amount of assimilated detritus in all
#' consumers of detritus and the biomass of the detritus compartment, and is an
#' implementation of the equation from Neutel & Thorne (2014). \cr
#' References: \cr
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' @param CR A named numeric vector with detritus consumption rates of all compartments
#' feeding on the detritus compartment. (required)
#' @param AE A named numeric vector with assimilation efficiencies of all compartments
#' feeding on the detritus compartment. Must be in the same order as CR. Must be a
#' fraction i.e. between 0 and 1 (required)
#' @param BM The biomass of the detritus compartment.
#' @return This function returns a value which is the diagonal value for the
#' detritus compartment in the food web (per unit time). It is important to review the units of the
#' input data. If CR is biomass per unit time then BM must be just biomass. If CR is
#' biomass per unit time per surface area then BM must be biomass per surface area.
getDiagonalDetritus <- function(CR, BM, AE){
  result <- -sum(AE * CR) / BM
  #detritus_ingestion <- colSums(t(FM[dead,-dead])*AE[-dead])
  #result <- -detritus_ingestion/BM[dead]
  return(result)
}

#diagonal <- diag(FM)
#aii <- diagonalSpecies(flow_solutions = pars$X, BM, dead)
#add <- diagonalDetritus(FM[-externals, -externals], BM, AE, dead)
#diagonal[names(aii)] <- aii
#diagonal[names(add)] <- add
#diag(JM2) <- diagonal[-externals]
#JM2[is.na(JM2)] <- 0
