# Output is same as flowmatrix
getMortalityFlow <- function(FM, AE, GE, dead_names) {
  dead_i <- which(rownames(FM) %in% dead_names)
  production <- AE*GE*colSums(FM, na.rm = T)
  predation <- rowSums(FM[,-dead_i], na.rm = T)
  m <- production - predation
  return(m)
}

#' Calculate mortality rates from energy-fluxes
#'
#' This function can be used to calculate mortality rates (unit t-1) based on a Flowmatrix,
#' conversion efficiencies, and biomasses.
#' @param FM (required) A named square flowmatrix, source compartments as rows,
#' sink compartments as columns. Should NOT contain external compartments.
#' @param AE (required) A named numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix after externals
#' are excluded. Must be a fraction i.e. between 0 and 1.
#' @param GE (required) A named numeric vector with growth (production) efficiencies of all compartments.
#' Must be a fraction i.e. between 0 and 1.
#' @param BM (required) A named numeric vector with biomasses of all compartments, must be in the same
#' order as the flow matrix after externals are excluded.
#' @param dead_names (optional) Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna).
#' @return Returns a names vector with mortality rates in the unit per time, t-1.
#' @export
getMortalityRates <- function(FM, AE, GE, BM, dead_names) {
  checkMformat(FM)
  checkBMformat(BM)
  m <- getMortalityFlow(FM, AE, GE, dead_names)
  MR <- m / BM
  return(MR)
}
