#' Per capita effect on resource population.
#'
#' This function finds the per capita effect of a consumer on the growth
#' of the resource population. It implements the functions from de Ruiter et al.
#' (1995), Neutel et al. (2002) and Neutel & Thorne (2014). The flow
#' from prey to predator is divided by the biomass of the consumer.
#' The flow from consumer to detritus, plus the defecation of predators of
#' the consumer into the detritus pool, minus the uptake of detritus, is
#' divided by the biomass of the consumer. \cr
#' References: \cr
#' de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of
#' Interaction Strengths, and Stability in Real Ecosystems. Science. 269,
#' 1257–1260. https://doi.org/10.1126/science.269.5228.1257 \cr
#' Neutel, A.M., Heesterbeek, J.A.P., Ruiter, P.C. De, 2002. Stability in Real Food Webs:
#' Weak Links in Long Loops. Science. 296, 1120–1123. https://doi.org/10.1126/science.1068326 \cr
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' @param FM A square flowmatrix, source compartments as rows,
#' sink compartments as columns. (required)
#' @param BM Numeric vector with biomasses of all compartments,
#' must be in the same order as the flow matrix. (required)
#' @param AE Numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix. (required)
#' @param dead Integer vector with indices of all detritus and nutrient
#' compartments (everything that is not fauna). (optional)
#' @return This function returns a matrix containing interaction strengths - measured as the
#' effect of the consumer (rows) on the resource (columns) - for all
#' interactions in the food web.
#' @export
#' @examples
#' effectOnResource(Flowmatrix(model), BM = c(10,20,30), AE = c(0.2, 0.6, 0.5),
#' dead = "DETRITUS")
effectOnResource <- function(FM, BM, AE, dead = NULL){
  # As a matrix is by default divided by row (which now contain the sources,
  # but we want to divide by consumer), the matrix must be transposed first.
  result <- -(t(FM) / BM)
  # eq. 12
  if(!is.null(dead)){
    dead_interactions <- expand.grid(rownames(FM)[-dead], rownames(FM)[dead])
    for(i in 1:dim(dead_interactions)[1]){
      consumer <- as.character(dead_interactions[i,1])
      resource <- as.character(dead_interactions[i,2])

      result[consumer, resource] <-
        (FM[consumer, resource] -
        FM[resource,consumer] +
        sum(FM[consumer,-dead]*(1-AE[-dead]))) /
        BM[consumer]

      # Consumer of detritus is not a resource if it deposits detritus.
      result[resource, consumer] <- NA
    }
  }

  return(result)
}

#' Per capita effect on consumer population.
#'
#' This function finds the per capita effect of a resource on the growth
#' of the consumer population. It uses functions from de Ruiter et al.
#' (1995) and Neutel et al. (2002). The flow from resource to consumer is
#' adjusted for conversion efficiency and divided by the biomass of the
#' resource. \cr
#' References: \cr
#' de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of
#' Interaction Strengths, and Stability in Real Ecosystems. Science. 269,
#' 1257–1260. https://doi.org/10.1126/science.269.5228.1257 \cr
#' Neutel, A.M., Heesterbeek, J.A.P., Ruiter, P.C. De, 2002. Stability in Real Food Webs:
#' Weak Links in Long Loops. Science. 296, 1120–1123. https://doi.org/10.1126/science.1068326 \cr
#' @param FM A square flowmatrix, source compartments as rows,
#' sink compartments as columns. (required)
#' @param BM Numeric vector with biomasses of all compartments, must be in the same
#' order as the flow matrix. (required)
#' @param AE Numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix. (required)
#' @param GE Numeric vector with growth efficiencies of all compartments,
#' must be in the same order as the flow matrix. (required)
#' @return This function returns a matrix containing interaction strengths - measured as the
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
#' from the energy-flux model. This function is an implementation of the equations
#' from the following references: \cr
#' de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of
#' Interaction Strengths, and Stability in Real Ecosystems. Science. 269,
#' 1257–1260. https://doi.org/10.1126/science.269.5228.1257 \cr
#' Neutel, A.M., Heesterbeek, J.A.P., Ruiter, P.C. De, 2002. Stability in Real Food Webs:
#' Weak Links in Long Loops. Science. 296, 1120–1123. https://doi.org/10.1126/science.1068326 \cr
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' @param FM A named square flowmatrix, source compartments as rows,
#' sink compartments as columns. (required)
#' @param BM Named numeric vector with biomasses of all compartments, must be in the same
#' order as the flow matrix. (required)
#' @param AE Named numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix.
#' AE should be set to NA for dead/non-faunal compartments (see argument 'dead' below). (required)
#' @param GE Named numeric vector with growth efficiencies of all compartments,
#' must be in the same order as the flow matrix. GE should be set to NA
#' for dead/non-faunal compartments (see argument 'dead' below). (required)
#' @param diagonal Either a single value, a numeric vector, or the
#' charcter string "model". A single value with set all diagonal
#' values to this number, a vector will set the diagonal to this
#' user-specified diagonal, "model" will calculate diagonal values from
#' the energy-flux model. Default is a zero diagonal. (required)
#' @param dead Character vector with all names of detritus and nutrient
#' compartments (everything that is not fauna). (optional)
#' @param externals Character vector with all names of external
#' compartments, i.e. which have no biomass. (optional)
#' @return This function returns a matrix containing interaction strengths, i.e. the
#' effect of the resources (rows) on the consumers (columns) - for all
#' interactions in the food web.
#' @export
#' @examples
#' getJacobian(FM = Flowmatrix(lim), BM = lim$Components$val)
getJacobian <- function(FM, BM, AE, GE, diagonal = 0,
                        dead = NULL, externals = NULL) {
  # Do checks for required data formats
  if(dim(FM)[1] != dim(FM)[2]) {
    stop("flow matrix is not square")
  } else if((TRUE %in% is.na(BM)) | (TRUE %in% (BM <= 0)) | (!is.numeric(BM))) {
    stop("biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  } else if(is.null(rownames(FM)) | is.null(colnames(FM)) |
            is.null(names(BM)) | is.null(names(AE)) | is.null(names(GE))) {
    stop("all required vectors and matrices must be named")
  } else if(!all(names(BM) == rownames(FM)) | !all(names(BM) == colnames(FM)) |
            !all(names(BM) == names(AE))    | !all(names(BM) == names(GE))) {
    stop("the names and their order must be equal in all named vectors and matrices")
  } else if(FALSE %in% (dead %in% names(BM))) {
    stop("the names of the dead compartments are unknown")
  } else if(FALSE %in% (externals %in% names(BM))) {
    stop("the names of the external compartments are unknown")
  } else if(!is.null(dead)) {
    if(!all(is.na(AE[which(names(AE) == dead)])) |
       !all(is.na(GE[which(names(GE) == dead)]))) {
      AE[which(names(AE) == dead)] <- NA
      GE[which(names(GE) == dead)] <- NA
      warning("physiological values set to NA for dead compartments")
    }
  }

  # Remove external compartments, keep internals
  internals <- !(rownames(FM) %in% externals)
  FM_int <- FM[internals, internals]

  # Get indices of dead compartments
  if(is.null(dead)) {
    dead_i <- NULL
  } else {
    dead_i <- which(rownames(FM) %in% dead)
  }

   # Get interaction strengths
  eff.on.consumer <- effectOnConsumer(FM_int, BM, AE, GE)
  eff.on.resource <- effectOnResource(FM_int, BM, AE, dead_i)

  # Set interaction strength to 0 if NA
  a <- which(is.na(eff.on.consumer))
  if(length(a) > 0){
    eff.on.consumer[a] <- 0
  }
  b <- which(is.na(eff.on.resource))
  if(length(b) > 0){
    eff.on.resource[b] <- 0
  }

  JM <- eff.on.consumer + eff.on.resource
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
