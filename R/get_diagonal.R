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
getDiagonalSpecies <- function(MR, BM) {
  if(length(MR) != length(BM)) {
    stop("input vectors have unequal lengths")
  } else if(!is.numeric(MR) | !is.numeric(BM)) {
    stop("input vectors must be numeric")
  } else if(is.null(names(MR)) | is.null(names(BM))) {
    stop("input vectors must be named")
  } else if(!all(names(MR) == names(BM))) {
    stop("names of vectors do not match")
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
#' consumers of detritus and the biomass of the detritus compartments, and is an
#' implementation of the equation from Neutel & Thorne (2014). \cr
#' References: \cr
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' @param FM A named square flowmatrix, source compartments as rows,
#' sink compartments as columns. Should NOT contain external compartments. (required)
#' @param BM A named numeric vector with biomasses of all compartments, must be in the same
#' order as the flow matrix after externals are excluded. (required)
#' @param AE A named numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix after externals
#' are excluded. Must be a fraction i.e. between 0 and 1 (required)
#' @param dead Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna). (required)
#' @return This function returns a named numeric vector with the diagonal values for the
#' detritus compartments in the food web (per unit time). It is important to review the units of the
#' input data. If the FM is biomass per unit time then BM must be just biomass. If FM is
#' biomass per unit time per surface area then BM must be biomass per surface area.
getDiagonalDetritus <- function(FM, BM, AE, dead){
  dead_i <- which(rownames(FM) %in% dead)
  if(length(dead) == 1) {
    det_assimilation <- sum(t(FM[dead_i, -dead_i] * AE[-dead_i]), na.rm = T)
  } else {
    det_assimilation <- colSums(t(FM[dead_i, -dead_i] * AE[-dead_i]))
  }
  result <- -det_assimilation/BM[dead_i]
  return(result)
}


#' Calculate Jacobian diagonal from flux values
#'
#' This function finds diagonal values to be used in the Jacobian matrix based on flux
#' values. Diagonal values for species depend on their non-predatory mortality rate and
#' biomass. Diagonal values for dead (detritus) compartments are based on the total amount
#' of assimilated detritus in all consumers of detritus and the biomass of the detritus
#' compartments. This function implements the equations from from Neutel & Thorne (2014).
#' The input should be named according to the food web compartment names. \cr
#' References: \cr
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' @param MR A named numeric vector with non-predatory mortality rates for all
#' compartments (biomass per unit time or biomass per unit time per surface area).
#' (required)
#' @param BM A named numeric vector with biomasses of all compartments
#' (just biomass or biomass per unit area), must be in the same order as MR. (required)
#' @param dead Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna). (optional)
#' @param FM A named square flowmatrix, source compartments as rows,
#' sink compartments as columns. Should NOT contain external compartments. Order of names
#' must be equal to MR. (required when dead compartments are included)
#' @param AE A named numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as MR. Must be a fraction i.e. between 0 and 1.
#' (required when dead compartments are included)
#' @return This function returns a named numeric vector with diagonal values for the
#' species in the food web (per unit time). It is important to review the units of the
#' input data. If MR is biomass per unit time then BM must be just biomass. If MR is
#' biomass per unit time per surface area then BM must be biomass per surface area.
#' @export
getDiagonal <- function(MR, BM, dead = NULL, FM = NULL, AE = NULL) {

  if(!is.null(dead) & is.null(FM) | is.null(AE)) {
    stop("please provide all required data to calculate dead diagonal values")
  }

  # Get diagonal, and find detritus values if necessary
  if(!is.null(dead)) {
    diagonal <- diag(FM)
    add <- getDiagonalDetritus(FM = FM, BM = BM, AE = AE, dead = dead)
  } else {
    diagonal <- MR
  }

  # Find species values and add to diagonal
  aii <- getDiagonalSpecies(MR = MR, BM = BM)
  diagonal[names(aii)] <- aii
  # Add detritus values to diagonal if necessary
  if(!is.null(dead)) {
    diagonal[names(add)] <- add
  }

  return(diagonal)
}

