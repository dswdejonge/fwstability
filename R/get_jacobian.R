#' Per capita effect on resource population.
#'
#' This function finds the per capita effect of a consumer on the growth
#' of the resource population. It implements the functions from de Ruiter et al.
#' (1995), Neutel et al. (2002) and Neutel & Thorne (2014). The flow
#' from prey to predator is divided by the biomass of the consumer.
#' The flow from consumer to detritus, plus the defecation of predators of
#' the consumer into the detritus pool, minus the uptake of detritus, is
#' divided by the biomass of the consumer.
#' @references
#' \itemize{
#' \item de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of
#' Interaction Strengths, and Stability in Real Ecosystems. Science. 269,
#' 1257–1260. https://doi.org/10.1126/science.269.5228.1257
#' \item Neutel, A.M., Heesterbeek, J.A.P., Ruiter, P.C. De, 2002. Stability in Real Food Webs:
#' Weak Links in Long Loops. Science. 296, 1120–1123. https://doi.org/10.1126/science.1068326
#' \item Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' }
#' @param FMs (required) A list with two elements \code{original} and \code{netMatrix}, both
#' containing a square flowmatrix, with source compartments as rows, sink compartments as columns.
#' \itemize{
#' \item{
#' the list element \code{original} should contain the original flow matrix as can be obtained
#' with the function \code{getFlowMatrix},
#' i.e. can contain a flow from A to B, and a flow from B to A.
#' }
#' \item{
#' the list element \code{netMatrix} should contain an adjusted version of the original matrix with only
#' net flows, which can be obtained with the function \code{getNetMatrixFM},
#' i.e. a netMatrix flow is the absolute result of flow A to B minus flow B to A in the correct direction.
#' The net flow is only calculated for flows between live compartments, not for dead ones.
#' }
#' }
#' @param BM (required) Numeric vector with biomasses of all compartments,
#' must be in the same order as the flow matrix.
#' @param AE (required) Numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix.
#' @param dead (optional) List with at most two elements named \emph{names} and
#' \emph{frac} containing information on all dead compartments (like detritus and nutrients).
#' \itemize{
#' \item The element \emph{names} is required and contains a character vector with all names of dead
#' compartments.
#' \item The element \emph{frac} is required if there is more than one detritus compartment.
#' It is matrix the same size and order
#' as the \code{FM} matrix, and contains the fraction of each flow that is defecation.
#' If there are multiple defecation compartments, the Flowmatrix \code{FM} combined with \code{frac} is used to
#' calculate the relative distribution of matter into defecation compartments.
#' }
#' @return This function returns a matrix containing the effects of the consumers (rows) on the
#' resources (columns).
#' @export
effectOnResource <- function(FMs, BM, AE, dead = NULL){
  # As a matrix is by default divided by row (which now contain the sources,
  # but we want to divide by consumer), the matrix must be transposed first.
  result <- -(t(FMs$netMatrix) / BM)

  if(!is.null(dead)){
    dead_i <- which(rownames(FMs$original) %in% dead$names)
    dead_interactions <- expand.grid(rownames(FMs$original), rownames(FMs$original)[dead_i])
    # One or multiple defecation compartments?
    if(!is.null(dead$frac)){
      DFM <- FMs$original * dead$frac
      defecation_compartments <- colnames(dead$frac)[which(colSums(dead$frac, na.rm = T) > 0)]
    } else {
      DFM <- FMs$original
      defecation_compartments <- dead$names
    }

    for(i in 1:dim(dead_interactions)[1]){
      consumer <- as.character(dead_interactions[i,1])
      resource <- as.character(dead_interactions[i,2])
      if(consumer == resource) {next}
      # if consumer is not an organism, but detritus, there is detritus-detritus interaction
      if(consumer %in% dead$names){
        det_det_interaction <- TRUE
      }else{
        det_det_interaction <- FALSE
      }
      # direct input to detritus (resource) by consumer
      a <- FMs$original[consumer, resource]
      # uptake of detritus by consumer
      b <- FMs$original[resource, consumer]
      if(det_det_interaction) {b <- 0}
      # consumer-based defecation by other compartments
      c <- FMs$original[consumer,-dead_i]*(1-AE[-dead_i])
      c <- c[which(c > 0)]
      if(!(length(c) > 0)) {
        c <- 0
      }
      # distribution of defecation flow to different dead compartments
      predators <- names(which(FMs$original[consumer,-dead_i] > 0))
      d <- DFM[predators, resource] /
              rowSums(DFM[predators, defecation_compartments,drop=F],  na.rm = T)
      if(!(length(d) > 0)) {
        d <- 0
      }

      value <- (a - b + sum(c * d, na.rm = T)) / BM[consumer]
      result[consumer, resource] <- value

      # Consumer of detritus is not a resource if it deposits detritus.
      if(!det_det_interaction) {result[resource, consumer] <- NA}
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
#' resource.
#' @references
#' \itemize{
#' \item de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of
#' Interaction Strengths, and Stability in Real Ecosystems. Science. 269,
#' 1257–1260. https://doi.org/10.1126/science.269.5228.1257
#' \item Neutel, A.M., Heesterbeek, J.A.P., Ruiter, P.C. De, 2002. Stability in Real Food Webs:
#' Weak Links in Long Loops. Science. 296, 1120–1123. https://doi.org/10.1126/science.1068326
#' }
#' @param FM (required) A square flowmatrix, source compartments as rows,
#' sink compartments as columns.
#' @param BM (required) Numeric vector with biomasses of all compartments, must be in the same
#' order as the flow matrix \code{FM}.
#' @param AE (required) Numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix \code{FM}.
#' @param GE (required) Numeric vector with growth efficiencies of all compartments,
#' must be in the same order as the flow matrix \code{FM}.
#' @return This function returns a matrix containing the effects of resources (rows) on
#' consumers (columns).
#' @seealso \code{\link{effectOnResource}}, \code{\link{getJacobianEnergyFlux}}
#' @export
effectOnConsumer <- function(FM, BM, AE, GE) {
  # Conversion efficiencies of the predators must be included, which are in the
  # columns. So, transposition is needed before multiplying for AE and GE.
  # Finally, the matrix is transposed back to its original form.
  result <- t(t(FM / BM)*AE*GE)
  return(result)
}

removeExternals <- function(externals, FM) {
  internals <- !(rownames(FM) %in% externals)
  FM <- FM[internals, internals]
  return(FM)
}

#' Jacobian matrix with interaction strengths
#'
#' This functions calculates interaction strengths from a resolved energy-flux
#' food web model and uses these values as entries for a Jacobian matrix.
#' The diagonal can be calculated from the model (default), all-zero, or user-defined.
#' This function is an implementation of the equations from the references listed below.
#' @references
#' \itemize{
#' \item de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of
#' Interaction Strengths, and Stability in Real Ecosystems. Science. 269,
#' 1257–1260. https://doi.org/10.1126/science.269.5228.1257
#' \item Neutel, A.M., Heesterbeek, J.A.P., Ruiter, P.C. De, 2002. Stability in Real Food Webs:
#' Weak Links in Long Loops. Science. 296, 1120–1123. https://doi.org/10.1126/science.1068326
#' \item Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles
#' and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' }
#' @param FM (required) A named square flowmatrix, source compartments as rows,
#' sink compartments as columns.
#' @param BM (required) A named numeric vector with biomasses of all compartments, must be in the same
#' order as the flow matrix \code{FM} after externals are excluded (see \code{externals}).
#' @param AE (required) A named numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix \code{FM} after externals
#' are excluded. \code{AE} should be set to NA for dead/non-faunal compartments
#' (see \code{dead}).
#' @param GE (required) A named numeric vector with growth efficiencies of all compartments,
#' must be in the same order as the flow matrix \code{FM} after externals are excluded.
#' \code{GE} should be set to NA for dead/non-faunal compartments (see \code{dead} and
#' \code{externals}).
#' @param diagonal (required) Either the string "model" (default), a single value, or a numeric vector.
#' \itemize{
#' \item The string "model" calculates the diagonal values from flux values.
#' \item A single value will set all diagonal values to this number.
#' \item A vector will set the diagonal to this user-specified diagonal.
#' }
#' @param dead (optional) List with at most two elements named \emph{names} and \emph{frac}
#' containing information on all dead compartments (like detritus and nutrients).
#' \itemize{
#' \item The element \emph{names} is required and contains a character vector with all names of dead
#' compartments.
#' \item The element \emph{frac} is required if there are multiple dead compartments, and
#' is a matrix the same size and order as the \code{FM} matrix with the fraction of each flow that is defecation.
#' If there are multiple defecation compartments, the \code{FM} combined with \code{frac} is used to
#' calculate the relative distribution of matter into the specified defecation
#' compartments.
#' }
#' @param externals (optional) Character vector with all names of external
#' compartments, i.e. which have no biomass, that have to be removed from
#' the flow matrix before calculations.
#' @param MR (optional) Mortality rates for all compartments per unit time (t-1).
#' Default behaviour (MR = NULL) is calculation of mortality rates as (growth - predation)/biomass.
#' Otherwise, you can provide a named numeric vector.
#' @param verbose (optional) Default is TRUE. Wether or not to print messages.
#' @param netMatrix (optional) Boolean. Default is TRUE: the netMatrix is used
#' to calculate interaction strengths. This is only relevant if there are two food web compartments
#' which act both as prey and predators to one another.
#' @details If \code{MR} is set to NULL, it is calculated from the model as the non-predatory mortality:
#' production minus predation divided by biomass.
#' Otherwise \code{MR} can sometimes be extracted from the food web model, for example when
#' natural death results in a flux from the faunal compartment to a carcass compartment.
#' Diving the mortality flux by the biomass of the faunal compartment gives the specific
#' mortality rate per unit time.
#' Mortality rate can also be calculated as the inverse of the natural lifespan of the species
#' (per unit time).
#' @return This function returns a matrix containing interaction strengths, i.e. the
#' effect of the resources (rows) on the consumers (columns) - for all
#' interactions in the food web.
#' @export
getJacobianEnergyFlux <- function(FM, BM, AE, GE, diagonal = "model",
                        dead = NULL, externals = NULL, MR = NULL,
                        verbose = T, netMatrix = TRUE) {

  # Remove externals
  if(!is.null(externals)) {
    checkExternalsFormat(externals, FM)
    FM <- removeExternals(externals, FM)
  }
  # Use netMatrix FM if necessary
  checkMformat(FM)
  checkDeadFormat(dead, FM)
  FMs <- list()
  if(netMatrix){
    FMs$original <- FM
    FMs$netMatrix <- getNetMatrixFM(FM, dead$names)
  } else {
    FMs$original <- FM
    FMs$netMatrix <- FM
  }
  # Do checks for required data formats
  checkMformat(FMs$netMatrix)
  checkNamingFormat(
    matrices = list(FMs$original, FMs$netMatrix, dead$frac),
    vectors = list(BM, AE, GE))
  # Force correct order
  fwnames <- rownames(FM)
  FMs$original <- FMs$original[fwnames, fwnames] # order cols, rows
  FMs$netMatrix <- FMs$netMatrix[fwnames, fwnames] # order cols, rows
  BM <- BM[fwnames]
  AE <- AE[fwnames]
  GE <- GE[fwnames]
  if(!is.null(dead)){
    dead$frac <- dead$frac[fwnames,fwnames]
  }
  if(!is.null(MR)){
    MR <- MR[fwnames]
  }
  checkBMformat(BM)
  checkDiagonalFormat(diagonal, correct_length = length(BM))
  checkCEformat(CE = list(AE, GE))
  # Set AE and GE to NA for dead compartments if necessary
  if(!is.null(dead)) {
    if(!all(is.na(AE[dead$names])) |
       !all(is.na(GE[dead$names]))) {
      AE[dead$names] <- NA
      GE[dead$names] <- NA
      warning("physiological values set to NA for dead compartments")
    }
  }

  # Get interaction strengths
  eff.on.consumer <- effectOnConsumer(FMs$netMatrix, BM, AE, GE)
  eff.on.resource <- effectOnResource(FMs, BM, AE, dead)

  # Set interaction strength to 0 if NA
  a <- which(is.na(eff.on.consumer))
  if(length(a) > 0){
    eff.on.consumer[a] <- 0
  }
  b <- which(is.na(eff.on.resource))
  if(length(b) > 0){
    eff.on.resource[b] <- 0
  }

  # Combine matrices and add diagonal
  JM <- eff.on.consumer + eff.on.resource
  if(all(diagonal == "model")) {
    if(is.null(MR)){
      MR <- getMortalityRates(FM = FM, AE = AE, GE = GE, BM = BM, dead = dead$names)
    }
    diagonal <- getDiagonal(MR = MR, BM = BM, dead = dead$names, FM = FMs$original, AE = AE)
  }
  diag(JM) <- diagonal

  return(JM)
}

getNetMatrixFM <- function(FM, deadnames) {
  netMatrix <- FM - t(FM)
  netMatrix[which(netMatrix < 0)] <- 0
  if(!is.null(deadnames)){
    netMatrix[deadnames,] <- FM[deadnames,]
    netMatrix[,deadnames] <- FM[,deadnames]
  }
  return(netMatrix)
}

#' Jacobian matrix from an energy-flux model
#'
#' This is a wrapper function that reviews the given energy-flux food-web model and redirects
#' the input to the correct function for obtaining interaction strengths in a Jacobian matrix.
#' @param model (required) A named list containing elements with food web model data.
#' One list element named \code{type} denoting the type of input model must exist
#' and should either be the string "EF", or "LIM".
#' \itemize{
#' \item \bold{"EF"} should be used if the model is a quantified energy flux model. Other list elements required
#' for EF models are \code{FM}, \code{BM}, \code{AE}, and \code{GE} to be used by the
#' \code{\link{getJacobianEnergyFlux}} function.
#' \item \bold{"LIM"} should be used if the model is a linear inverse model created
#' with the package LIM and is formatted to contain tags.
#' The other list element required for LIM type models is \code{LIM} containing a read-in LIM.
#' If the LIM is resolved the flow solutions can be provided in the list element \code{web}.
#' If the LIM is not resolved, the function will use the parsimonious (least distance) solution.
#' }
#' @param diagonal (required) is by default calculated from the model (= "model").
#' Can also be set to an all-zero diagonal (= 0) or a numeric vector.
#'
#' @return This function returns a matrix containing interaction strengths, i.e. the
#' effect of the resources (rows) on the consumers (columns) - for all
#' interactions in the food web.
#' If the model input is a LIM the output is a list \code{JM} with element \code{JM} (so \code{JM$JM})
#' that is the Jacobian matrix with interaction strengths, and with element \code{extracted_data}
#' (so \code{JM$extracted_data}) that is another list with all data extracted from the LIM the
#' Jacobian matrix is calculated from.
#'
#' @section Energy Flux models:
#' Interaction strenghts can be calculated from energy flux models as presented in
#' De Ruiter (1995) and Neutel & Thorne (2014). The energy flux model should at least include
#' quantified flows, compartment biomasses, and conversion efficiencies.
#' If you want to calculate the intraspecific interactions, i.e. the diagonal values, from the
#' model, then non-predatory mortality rates should also be included. Beware that these model-derived
#' diagonal values represent an \bold{upperbound}. It is the maximum amount of self-dampening possible
#' in the respective populations.
#' If your model includes dead compartments, like detritus or nutrients, it is highly recommended to
#' also include these in the special argument \code{dead}. Conversion efficiences are not needed for
#' dead compartments, and you can include information on defecation into these compartments.
#' Interaction strengths of dead compartments are calculated in another way than for living compartments.
#' Inclusion of recycling processes may impact your stability results (Wilson & Wolkovich, 2011).\cr
#' The list with model information should contain:
#' \itemize{
#' \item{\code{type} (required) should be "EF" for energy flux models.}
#' \item{\code{FM} (required) is a named square flowmatrix with sources in rows and sinks in columns.}
#' \item{\code{BM} (required) is a named numeric with biomasses for all compartments.}
#' \item{\code{AE} (required) is a named numeric with assimilation efficiencies for living compartments.}
#' \item{\code{GE} (required) is a named numeric with growth efficiencies for living compartments.}
#' \item{\code{dead} (optional) is a list with the elements \code{names} and \code{frac} that
#' contains the names of dead compartments (character vector) and what fraction
#' of each flow comprises defecation (matrix similar to \code{FM} after removal of externals, only required with parallel flows).}
#' \item{\code{externals} (optional) is a character vector with any compartments that should not be
#' considered in the calculations}.
#' \item{\code{MR}} (optional) Mortality rates for all compartments per unit time (t-1).
#' Default behaviour (MR = NULL) is calculation of mortality rates as (growth - predation)/biomass.
#' Otherwise, you can provide a named numeric vector.
#' }
#' More detailed information about argument requirements see \code{?getJacobianEnergyFlux}.
#'
#' @section Linear Inverse models:
#' Linear programming can be used to quantify fluxes. Therefore, a solved linear inverse model (LIM)
#' often already contains the elements needed to calculate interaction strengths, but they need to be
#' extracted and calculated first. A LIM can be created with the R-package LIM. By sticking to certain
#' rules when writing the model input file, the \code{fwstability} package can automatically derive the
#' necessary information for stability analysis. These requirements include the use of tags in the
#' names of compartments, variables, and flows to later identify dead compartments, defecation and
#' mortality flows, and variables representing assimilation and growth. Tags should be used in combination
#' with only the compartmentname to assign the extracted value to the right compartment (exception is
#' for defecation flows). Tags are not case-sensitive and can be placed either in front or after the
#' compartmentname. Default tags are "dead", "ass", "growth", "mort", and "def", but you can also choose
#' your own tag names. More information on setting up the input model with these tags can be found in
#' the vignette. \cr
#' The LIM can be solved and quantified in various ways (e.g. by optimization or maximum likelihood)
#' and added to the model list in the element \code{web} as named numeric.
#' If you do not solve your LIM  before parsing it into this function, the model will be solved using
#' Ldei, which is the parsimonious solution minimizing the sum of squares.\cr
#' The list with model information should contain:
#' \itemize{
#' \item{\code{type} (required) should be "LIM"}
#' \item{\code{LIM} (required) should be a read-in LIM. This can be acquired by Read(<path-to-input-file>) from the
#' LIM R-package.}
#' \item{\code{setup} (optional) can be Setup(Read(<path-to-input-file>)), otherwise set-up within function.}
#' \item{\code{web} (optional) named numeric with flow solutions.}
#' \item{\code{aTag}, \code{gTag}, \code{mTag}, \code{defTag}, \code{deadTag} are optional. Default tags are
#' "ass", "growth", "mort", "def", and "dead" respectively.}
#' }
#' More information on the requirments of the arguments see \code{?extractLIMdata}
#' @seealso \code{\link{getJacobianEnergyFlux}}, \code{\link{extractLIMdata}}.
#' @references \itemize{
#' \item{de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of Interaction Strengths, and Stability in Real Ecosystems. Science (80-. ). 269, 1257–1260. https://doi.org/10.1126/science.269.5228.1257
#' }
#' \item{Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles and the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17, 651–661. https://doi.org/10.1111/ele.12266
#' }
#' \item{Wilson, E.E., Wolkovich, E.M., 2011. Scavenging: How carnivores and carrion structure communities. Trends Ecol. Evol. 26, 129–135. https://doi.org/10.1016/j.tree.2010.12.011
#' }
#' \item{\code{\link{LIM}} package, Soetaert & van Oevelen 2015.}
#' }
#' @export
getJacobian <- function(model = stop("Model input required"),
                        diagonal = "model", verbose = T, netMatrix = T) {
  if(model$type == "EF") {
    JM <- getJacobianEnergyFlux(
      FM = model$FM,
      BM = model$BM,
      AE = model$AE,
      GE = model$GE,
      dead = model$dead,
      externals = model$externals,
      MR = model$MR,
      netMatrix = netMatrix,
      verbose = verbose,
      diagonal = diagonal)
  } else if(model$type == "LIM") {
    if(is.null(model$setup)) {
      model$setup <- Setup(model$LIM)
    }
    if(is.null(model$web)) {
      if(verbose) {
        message("fwstab: No model solutions given, LIM resolved by minimizing sum of squares.")
      }
      model$web <- LIM::Lsei(model$setup, parsimonious = TRUE)$X
    } else if(!is.numeric(model$web) | is.null(names(model$web))) {
      stop("Model solutions in \"web\" must be named numeric vector.")
    }
    extracted_data <- extractLIMdata(model)
    JM <- getJacobianEnergyFlux(
      FM = extracted_data$FM,
      BM = extracted_data$BM,
      AE = extracted_data$AE,
      GE = extracted_data$GE,
      dead = extracted_data$dead,
      externals = extracted_data$externals,
      MR = extracted_data$MR,
      netMatrix = netMatrix,
      verbose = verbose,
      diagonal = diagonal
    )
    JM <- list(JM = JM, extracted_data = extracted_data)
  } else {
    stop("Unknown model input")
  }
  return(JM)
}

#' Normalize the Jacobian matrix
#'
#' This function normalizes the Jacobian matrix by dividing each row by
#' the corresponding absolute diagonal value.
#'
#' @param JM (required) A square named Jacobian matrix with numeric values representing
#' the effect of one compartment (rows) on another compartment (columns).
#' The diagonal must be quantified, and cannot contain zeroes.
#' @param allzero Boolean. Default TRUE. Sets all diagonal values - except for dead compartments -
#' to zero after normalisation. If set to FALSE, the normalized Jacobian will have -1 for each diagonal
#' value.
#' @param dead_names (required if allzero is TRUE and there are dead compartments in the system)
#' Character vector with the names of all dead compartments.
#' @return Returns a normalized Jacobian matrix with either an all-zero (except detritus) diagonal
#' or a diagonal with -1.
#' @export
normalizeJacobian <- function(JM, dead_names = NULL, allzero = TRUE){
  checkMformat(JM)
  if(0 %in% diag(JM)){
    stop("No zeroes may be present on the diagonal when normalizing the Jacobian matrix.")
  }
  JMnorm <- JM / abs(diag(JM))
  if(allzero){
    i <- names(diag(JMnorm)) %in% dead_names
    diag(JMnorm)[!i] <- 0
  }
  return(JMnorm)
}
