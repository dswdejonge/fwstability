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
#' @param FM (required) A square flowmatrix, source compartments as rows,
#' sink compartments as columns.
#' @param BM (required) Numeric vector with biomasses of all compartments,
#' must be in the same order as the flow matrix.
#' @param AE (required) Numeric vector with assimilation efficiencies of all
#' compartments, must be in the same order as the flow matrix.
#' @param dead (optional) List with at most three elements named \emph{names}, \emph{def},
#' and \emph{frac} containing information on all dead compartments (like detritus and nutrients).
#' \itemize{
#' \item The element \emph{names} is required and contains a character vector with all names of dead
#' compartments.
#' \item The element \emph{def} should be either NULL or contain a character vector
#' specifying for each dead compartment if defecation occurs into this compartment
#' (\bold{Def}) or not (\bold{noDef}).
#' If this element is NULL it is assumed no defecation occurs into all dead compartments.
#' If there are multiple defecation compartments, the Flowmatrix \code{FM} is used to
#' calculate the relative distribution of matter into the specified defecation
#' compartments.
#' \item The element \emph{frac} should be either NULL or contain a matrix the same size and order
#' as the \code{FM} matrix with the fraction of each flow that is defecation.
#' This information is only needed if there are multiple parallel flows between two compartments of
#' which only one reflects actual defecation, and the rest is for example mortality.
#' }
#' @return This function returns a matrix containing the effects of the consumers (rows) on the
#' resources (columns).
#' @export
effectOnResource <- function(FM, BM, AE, dead = NULL){
  # As a matrix is by default divided by row (which now contain the sources,
  # but we want to divide by consumer), the matrix must be transposed first.
  result <- -(t(FM) / BM)

  if(!is.null(dead)){
    dead_i <- which(rownames(FM) %in% dead$names)
    dead_interactions <- expand.grid(rownames(FM)[-dead_i], rownames(FM)[dead_i])
    for(i in 1:dim(dead_interactions)[1]){
      consumer <- as.character(dead_interactions[i,1])
      resource <- as.character(dead_interactions[i,2])
      if(is.null(dead$def)) {
        is_defecation_compartment <- FALSE
      } else {
        is_defecation_compartment <- dead$def[which(dead$names == resource)] == "Def"
      }

      a <- FM[consumer, resource]
      b <- FM[resource,consumer]
      if(is_defecation_compartment) {
        c <- sum(FM[consumer,-dead_i]*(1-AE[-dead_i]), na.rm = T)
      } else {
        c <- 0
      }
      if(!is.null(dead$frac)) {
        DFM <- FM * dead$frac
      }else {
        DFM <- FM
      }
      if(length(which(dead$def == "Def")) > 1) {
        defecation_compartments <- dead$names[which(dead$def == "Def")]
        predators <- which(FM[consumer,] > 0)[-dead_i]
        d <- sum(DFM[predators,resource], na.rm = T) /
             sum(DFM[predators,defecation_compartments],  na.rm = T)
        if(is.na(d)) {
          d <- 1
        }
      } else {
        d <- 1
      }

      result[consumer, resource] <- (a - b + c * d) / BM[consumer]

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
  if(!is.null(externals)) {
    if((FALSE %in% (externals %in% rownames(FM))) |
       (FALSE %in% (externals %in% colnames(FM)))) {
      stop("the names of the external compartments are unknown")
    } else {
      # Remove external compartments, keep internals
      internals <- !(rownames(FM) %in% externals)
      FM <- FM[internals, internals]
    }
  }
  return(FM)
}

adjustDeadInput <- function(dead) {
  if(!is.null(dead)) {
    if(!is.list(dead) | is.null(names(dead))) {
      stop("argument \"dead\" must be a named list")
    }
    if(is.null(dead$names)) {
      stop("\"names\" element is required in the \"dead\" list")
    }
    names <- c("names", "def", "frac")
    missing <- c(1:length(names))[-which(names %in% names(dead))]
    if(length(dead) < length(names)) {
      newnames <- c(names(dead), names[missing])
      dead <- c(dead, vector(
        mode = "list", length = length(names) - length(dead)))
      names(dead) <- newnames
    } else if(length(dead) > length(names)) {
      stop(paste("the list \"dead\" should have",length(names),"elements at most"))
    }
    if(!is.null(dead$def) &&
       length(dead$names) !=
              length(which(dead$def == "Def" | dead$def == "noDef"))) {
      stop("the second element of the list \"dead\" may only contain the strings \"Def\" and \"noDef\"")
    }
  }
  return(dead)
}

#' Jacobian matrix with interaction strengths
#'
#' This functions calculates interaction strengths from a resolved energy-flux
#' food web model and uses these values as entries for a Jacobian matrix.
#' The diagonal can be all-zero (default), user-defined, or calculated from the model.
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
#' @param diagonal (required) Either a single value, a numeric vector or the string "model".
#' Default is an all-zero diagonal.
#' \itemize{
#' \item A single value will set all diagonal values to this number.
#' \item A vector will set the diagonal to this user-specified diagonal.
#' \item The string "model" calculates the diagonal values from flux values.
#' This requires the argument \code{MR}.
#' }
#' @param dead (optional) List with at most three elements named \emph{names}, \emph{def},
#' and \emph{frac} containing information on all dead compartments (like detritus and nutrients).
#' \itemize{
#' \item The element \emph{names} is required and contains a character vector with all names of dead
#' compartments.
#' \item The element \emph{def} of the list is optional and can contain a character vector
#' specifying for each dead compartment if defecation occurs into this compartment
#' (\bold{Def}) or not (\bold{noDef}).
#' If this element is omitted it is assumed no defecation occurs into all dead compartments.
#' If there are multiple defecation compartments, the Flowmatrix is used to
#' calculate the relative distribution of matter into the specified defecation
#' compartments.
#' \item The element \emph{frac} is optional and should contain a matrix the same size and order
#' as the \code{FM} matrix with the fraction of each flow that is defecation.
#' This information is only needed if there are multiple parallel flows between two compartments of
#' which only one reflects actual defecation, and the rest is for example mortality.
#' }
#' @param externals (optional) Character vector with all names of external
#' compartments, i.e. which have no biomass, that have to be removed from
#' the flow matrix before calculations.
#' @param MR (required when \code{diagonal} is set to \emph{model})
#' A named numeric vector with non-predatory mortality rates for all
#' compartments (same units as the flow matrix \code{FM}).
#' @details \code{MR} can sometimes be extracted from the food web model, for example when
#' natural death results in a flux from the faunal compartment to a carcass compartment.
#' It can also be calculated as the inverse of the natural lifespan of the species
#' (per unit time) multiplied by the biomass of the compartment.
#' @return This function returns a matrix containing interaction strengths, i.e. the
#' effect of the resources (rows) on the consumers (columns) - for all
#' interactions in the food web.
#' @export
getJacobianEnergyFlux <- function(FM, BM, AE, GE, diagonal = NULL,
                        dead = NULL, externals = NULL, MR = NULL) {

  FM <- removeExternals(externals, FM)
  dead <- adjustDeadInput(dead)
  if(is.null(diagonal)) {
    diagonal <- 0
    message("fwstab: Diagonal by default set to all-zero.")}

  # Do checks for required data formats: throws errors
  if(dim(FM)[1] != dim(FM)[2]) {
    stop("flow matrix is not square")
  } else if(is.null(rownames(FM)) | is.null(colnames(FM)) |
            is.null(names(BM)) | is.null(names(AE)) | is.null(names(GE))) {
    stop("all required vectors and matrices must be named")
  } else if(!all(rownames(FM) == colnames(FM))) {
    stop("row names and column names of flow matrix do not match")
  } else if((TRUE %in% is.na(BM)) | (TRUE %in% (BM <= 0)) | (!is.numeric(BM))) {
    stop("biomass vector contains NA, values equal or smaller than zero, or is non-numeric")
  } else if(!all(names(BM) == rownames(FM)) | !all(names(BM) == colnames(FM)) |
            !all(names(BM) == names(AE))    | !all(names(BM) == names(GE))) {
    stop("the names and their order must be equal in all named vectors and matrices")
  } else if(FALSE %in% (dead$names %in% names(BM))) {
    stop("the names of the dead compartments are unknown")
  } else if(!is.numeric(diagonal) & all(diagonal != "model")) {
    stop("given diagonal not numeric or set to \"model\"")
  } else if(length(diagonal) != 1 & length(diagonal) != length(BM)) {
    stop("given diagonal has incorrect length")
  } else if(any(AE > 1 | AE < 0 | GE > 1 | GE < 0, na.rm = TRUE)) {
    stop("assimilation and growth efficiencies must lie between 0 and 1")
  } else if(!is.null(dead)) {
    if(!all(is.na(AE[names(AE) %in% dead$names])) |
       !all(is.na(GE[names(GE) %in% dead$names]))) {
      AE[names(AE) %in% dead$names] <- NA
      GE[names(GE) %in% dead$names] <- NA
      warning("physiological values set to NA for dead compartments")
    }
  }

  # Get interaction strengths
  eff.on.consumer <- effectOnConsumer(FM, BM, AE, GE)
  eff.on.resource <- effectOnResource(FM, BM, AE, dead)

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
    if(is.null(dead)) {
      diagonal <- getDiagonal(MR = MR, BM = BM)
    } else {
      diagonal <- getDiagonal(MR = MR, BM = BM,
                              dead = dead$names, FM = FM, AE = AE)
    }
  }
  diag(JM) <- diagonal

  return(JM)
}

#' Jacobian matrix from a ODE model
#'
#' This function produces a Jacobian matrix from a model defined as a set of
#' ordinary differential equations.
#' @inheritParams rootSolve::jacobian.full
#' @seealso The package rootSolve for full documentation on \code{\link[rootSolve]{jacobian.full}}.
getJacobianODE <- function(y, func, parms) {
  JM <- rootSolve::jacobian.full(y = y, func = func, parms = parms)
  rownames(JM) <- colnames(JM)
  return(JM)
}

#' Jacobian matrix from a food web model
#'
#' This is a wrapper function that reviews the given food web model and redirects
#' the input to the correct function for obtaining interaction strengths in a Jacobian matrix.
#' @param model (required) A named list containing elements with food web model data.
#' One list element named \code{type} denoting the type of input model must exist
#' and should either be the string "ODE", "EF", or "LIM".
#' \itemize{
#' \item \bold{"ODE"} should be used if the model is a set of
#' ordinary differential equations. Other list element required for ODE type models are
#' \code{func} and \code{y}. The element \code{parms} is optional. All data will be used
#' in the \code{\link{getJacobianODE}} function which relies on the rootSolve package.
#' \item \bold{"EF"} should be used if the model is a quantified energy flux model. Other list elements required
#' for EF models are \code{FM}, \code{BM}, \code{AE}, and \code{GE} to be used by the
#' \code{\link{getJacobianEnergyFlux}} function.
#' \item \bold{"LIM"} should be used if the model is a linear inverse model created
#' with the package LIM and is formatted to contain tags.
#' The other list element required for LIM type models is \code{LIM} containing a read-in LIM.
#' If the LIM is resolved the flow solutions can be provided in the list element \code{web}.
#' If the LIM is not resolved, the function will use the parsimonious (least distance) solution.
#' }
#' @return This function returns a matrix containing interaction strengths, i.e. the
#' effect of the resources (rows) on the consumers (columns) - for all
#' interactions in the food web.
#' @section ODE models:
#' A ODE type model should be set-up in such a way that it can be used by the \code{rootSolve} package.
#' For details, please refer to the documentation of this package.
#' \itemize{
#' \item{The element \code{type} should be "ODE".}
#' \item{The element \code{func} should contain a function that calculates the rate of change for all
#' compartments i.e. a set of ODEs.}
#' \item{The element \code{y} is a named vector with the initial state (i.e. biomass) of compartments,
#' which is needed by the function in \code{func}.
#' In this respect is the data in \emph{y} similar to the BM vector that needs to be provided
#' for the EF type models.}
#' \item{The element \code{parms} is optional and contain parameters that can be used by the function
#' in \code{func}.}
#' }
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
#' \item{\code{diagonal} (required) is by default an all-zero diagonal. Can also be set to other values
#' or calculated from the model by setting it to "model".}
#' \item{\code{dead} (optional) is a list with the elements \code{names}, \code{def}, \code{frac} that
#' contains the names of dead compartments (character vector), wether or not defecation occurs into these
#' dead compartments ("Def" or "noDef" in character vector in same order as names), and what fraction
#' of each flow comprises defecation (matrix similar to \code{FM}, only required with parallel flows).}
#' \item{\code{externals} (optional) is a character vector with any compartments that should not be
#' considered in the calculations}.
#' \item{\code{MR} (required if \code{diagonal} is set to "model") is a named numeric with non-predatory
#' mortality rates in the same units as the flows.}
#' }
#' More detailed information about argument requirements see \code{?getJacobianEnergyFlux}.
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
#' \item{\code{diagonal} (optional) all-zero by default, can also be other numbers or calculated from model
#' by setting this element to "model".}
#' \item{\code{aTag}, \code{gTag}, \code{mTag}, \code{defTag}, \code{deadTag} are optional. Default tags are
#' "ass", "growth", "mort", "def", and "dead" respectively.}
#' }
#' More information on the requirments of the arguments see \code{?extractLIMdata}
#' @seealso \code{\link{getJacobianODE}}, \code{\link{getJacobianEnergyFlux}},
#' \code{\link{extractLIMdata}}.
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
getJacobian <- function(model = stop("Model input required")) {
  if(model$type == "ODE") {
    JM <- getJacobianODE(
      y = model$y,
      func = model$func,
      parms = model$parms)
  } else if(model$type == "EF") {
    JM <- getJacobianEnergyFlux(
      FM = model$FM,
      BM = model$BM,
      AE = model$AE,
      GE = model$GE,
      diagonal = model$diagonal,
      dead = model$dead,
      externals = model$externals,
      MR = model$MR)
  } else if(model$type == "LIM") {
    if(is.null(model$setup)) {
      model$setup <- Setup(model$LIM)
    }
    if(is.null(model$web)) {
      message("fwstab: No model solutions given, LIM resolved by minimizing sum of squares.")
      model$web <- Ldei(model$setup)$X
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
      diagonal = model$diagonal
    )
  } else {
    stop("Unknown model input")
  }
  return(JM)
}
