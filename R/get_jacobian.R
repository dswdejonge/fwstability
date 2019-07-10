#' Per capita effect on resource population.
#'
#' This function finds the per capita effect of a consumer on the growth
#' of the resource population. It uses functions from de Ruiter et al.
#' (1995), Neutel et al. (2002) and Neutel & Thorne (2014). The flow
#' from prey to predator is divided by the biomass of the consumer.
#' The flow from consumer to detritus, plus the defecation of predators of
#' the consumer into the detritus pool, minus the uptake of detritus, is
#' divided by the biomass of the consumer.
#' @param FM A square flowmatrix, source compartments as rows,
#' sink compartments as columns.
#' @param BM Numeric vector with biomasses of all compartments, must be in the same order as the flow matrix.
#' @param AE Numeric vector with assimilation efficiencies of all
#' compartments,must be in the same order as the flow matrix.
#' @param dead Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna).
#' @return A matrix containing interaction strengths - measured as the
#' effect of the consumer (rows) on the resource (columns) - for all
#' interactions in the food web.
#' @export
#' @examples
#' effectOnResource(Flowmatrix(model), BM = c(10,20,30), AE = c(0.2, 0.6, 0.5), dead = "DETRITUS")
effectOnResource <- function(FM, BM, AE, dead){
  # As a matrix is by default divided by row (which now contain the sources,
  # but we want to divide by consumer), the matrix must be transposed first.
  result <- -(t(FM) / BM)
  # eq. 12
  dead_interactions <- expand.grid(rownames(FM)[-dead], rownames(FM)[dead])
  for(i in 1:dim(dead_interactions)[1]){
    consumer <- as.character(dead_interactions[i,1])
    resource <- as.character(dead_interactions[i,2])
    result[consumer, resource] <-
      FM[consumer, resource] +
      sum(FM[consumer,-dead]*AE[-dead]) -
      FM[resource,consumer]
    # Consumer of detritus is not seen as resource if it deposits detritus.
    result[resource, consumer] <- NA
  }
  return(result)
}

#' Per capita effect on consumer population.
#'
#' This function finds the per capita effect of a resource on the growth
#' of the consumer population. It uses functions from de Ruiter et al.
#' (1995) and Neutel et al. (2002). The flow from resource to consumer is
#' adjusted for conversion efficiency and divided by the biomass of the
#' resource.
#' @param FM A square flowmatrix, source compartments as rows,
#' sink compartments as columns.
#' @param BM Numeric vector with biomasses of all compartments, must be in the same order as the flow matrix.
#' @param AE Numeric vector with assimilation efficiencies of all
#' compartments,must be in the same order as the flow matrix.
#' @param GE Numeric vector with growth efficiencies of all compartments,
#' must be in the same order as the flow matrix.
#' @return A matrix containing interaction strengths - measured as the
#' effect of the resources (rows) on the consumers (columns) - for all
#' interactions in the food web.
#' @export
#' @examples
effectOnConsumer <- function(FM, BM, AE, GE){
  # Conversion efficiencies of the predators must be included, which are in the
  # columns. So, transposition is needed before multiplying for AE and GE.
  # Finally, the matrix is transposed back to its original form.
  result <- t(t(FM / BM)*AE*GE)
  return(result)
}

# equation 13; diagonal for species
# is the negative mortality rate of species
#diagonalSpecies <- function(flow_solutions, BM, dead){
#  result <- -getMortalityRate(flow_solutions, BM, dead)
#  result <- result[-dead]
#  return(result)
#}

# equation 14; diagonal for detritus
# the total assimilated detritus is all consumers,
# divided by the biomass of detritus
#diagonalDetritus <- function(FM, BM, AE, dead){
#  detritus_ingestion <- colSums(t(FM[dead,-dead])*AE[-dead])
#  result <- -detritus_ingestion/BM[dead]
#}

#' Jacobian matrix with interaction strengths
#'
#' This functions calculates interaction strengths from a resolved energy-flux
#' food web model and uses these values as entries for a Jacobian matrix.
#' The diagonal can be either be user-defined, set to zero (default), or calculated
#' from the energy-flux model.
#' @param FM A square flowmatrix, source compartments as rows,
#' sink compartments as columns.
#' @param BM Numeric vector with biomasses of all compartments, must be in the same order as the flow matrix.
#' @param AE Numeric vector with assimilation efficiencies of all
#' compartments,must be in the same order as the flow matrix.
#' AE should be set to NA for non-faunal compartments.
#' @param GE Numeric vector with growth efficiencies of all compartments,
#' must be in the same order as the flow matrix. GE should be set to NA
#' for non-faunal compartments.
#' @param dead Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna).
#' @param externals Character vector with all names of external
#' compartments, i.e. which have no biomass.
#' @param diagonal Either a single value, a numeric vector, or the
#' charcter string "model". A single value with set all diagonal
#' values to this number, a vector will set the diagonal to this
#' user-specified diagonal, "model" will calculate diagonal values from
#' the energy-flux model. Default is a zero diagonal.
#' @return A matrix containing interaction strengths, i.e. the
#' effect of the resources (rows) on the consumers (columns) - for all
#' interactions in the food web.
#' @export
#' @examples
#' getJacobian(FM = Flowmatrix(lim), BM = lim$Components$val)
getJacobian <- function(FM, BM, AE, GE, dead, externals,
                        diagonal = 0) {
  # Remove external compartments, keep internals
  internals <- !(rownames(flow_matrix) %in% externals)
  FM_int <- FM[internals, internals]

  # Get indices of dead compartments
  dead_i <- which(rownames(flow_matrix) %in% dead)

   # Get interaction strengths
  eff.on.pred <- effectOnConsumer(FM_int, BM, AE, GE)
  eff.on.prey <- effectOnResource(FM_int, BM, AE, dead_i)

  # Set interaction strength to 0 if NA
  a <- which(is.na(eff.on.pred))
  if(length(a) > 0){
    eff.on.pred[a] <- 0
  }
  b <- which(is.na(eff.on.prey))
  if(length(b) > 0){
    eff.on.prey[b] <- 0
  }

  JM <- eff.on.pred + eff.on.prey
  diag(JM) <- 0

  #diagonal <- diag(FM)
  #aii <- diagonalSpecies(flow_solutions = pars$X, BM, dead)
  #add <- diagonalDetritus(FM[-externals, -externals], BM, AE, dead)
  #diagonal[names(aii)] <- aii
  #diagonal[names(add)] <- add
  #diag(JM2) <- diagonal[-externals]
  #JM2[is.na(JM2)] <- 0
  return(JM)
}
