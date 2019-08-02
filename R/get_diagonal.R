#' Calculate diagonal value for species
#'
#' This function find the diagonal value for species based on their non-predatory
#' mortality rate and biomass. Based on the equation from Neutel & Thorne (2014). \cr
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

# equation 14; diagonal for detritus
# the total assimilated detritus is all consumers,
# divided by the biomass of detritus
#diagonalDetritus <- function(FM, BM, AE, dead){
#  detritus_ingestion <- colSums(t(FM[dead,-dead])*AE[-dead])
#  result <- -detritus_ingestion/BM[dead]
#}

#diagonal <- diag(FM)
#aii <- diagonalSpecies(flow_solutions = pars$X, BM, dead)
#add <- diagonalDetritus(FM[-externals, -externals], BM, AE, dead)
#diagonal[names(aii)] <- aii
#diagonal[names(add)] <- add
#diag(JM2) <- diagonal[-externals]
#JM2[is.na(JM2)] <- 0
