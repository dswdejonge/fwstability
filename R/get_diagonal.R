# Output is t-1
getMortalityRates <- function(FM, AE, GE, BM) {
  m <- AE*GE*Q - Predation
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
    det_assimilation <- colSums(t(FM[dead_i, -dead_i] * AE[-dead_i]), na.rm = T)
  }
  result <- -det_assimilation/BM[dead_i]
  return(result)
}


#' Calculate Jacobian diagonal from flux values
#'
#' This function finds diagonal values to be used in the Jacobian matrix based on flux
#' values.
#' @details Diagonal values for species are equal to their natural death rate (unit t-1).
#' Natural death rates can be provided directly as \code{MR},
#' or assumed to be non-predatory mortality calculated from consumption rates and
#' physiological parameters.
#' Diagonal values for dead (detritus) compartments are based on the total amount
#' of assimilated detritus in all consumers of detritus and the biomass of the detritus
#' compartments. This function implements the equations from from Neutel & Thorne (2014).
#' The input should be named according to the food web compartment names.
#' It is important to review the units of the input data.
#' @references Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' @param MR (required) A named numeric vector with non-predatory mortality rates for all
#' compartments (per unit time, t-1).
#' If not provided, the input of \code{BM}, \code{FM}, \code{AE}, and \code{GE} are required.
#' @param dead (optional) Character vector with all names of detritus and nutrient
#' compartments.
#' @param BM (required when \code{MR} is NULL, and when dead compartments are included)
#' A named numeric vector with biomasses of all compartments.
#' Must be in the same order as MR.
#' @param FM (required if \code{MR} is NULL and when dead compartments are included)
#' A named square flowmatrix, source compartments as rows, sink compartments as columns.
#' Should NOT contain external compartments.
#' Order of names must be equal to MR.
#' @param AE (required if \code{MR} is NULL and when dead compartments are included)
#' A named numeric vector with assimilation efficiencies of all compartments.
#' Must be in the same order as MR.
#' Must be a fraction i.e. between 0 and 1.
#' @param GE (required if \code{MR} is NULL)
#' A named numeric vector with growth (production) efficiencies of all compartments.
#' Must be a fraction i.e. between 0 and 1.
#' @return This function returns a named numeric vector with diagonal values for the
#' species in the food web (per unit time, t-1).
#' @export
getDiagonal <- function(MR, dead = NULL, BM = NULL, FM = NULL, AE = NULL) {

  # Species
  checkNamingFormat(vectors = list(MR))
  aii <- -MR

  # Dead compartments
  if(!is.null(dead)) {
    if(is.null(BM) | is.null(FM) | is.null(AE)) {
      stop("please provide all required data to calculate dead diagonal values")
    }
    add <- getDiagonalDetritus(FM, BM, AE, dead)
    aii[names(add)] <- add
  }
  return(aii)
}

