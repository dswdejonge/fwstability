#' Get a Flowmatrix from a LIM
#'
#' This function extracts flow data from the LIM and produces a Flowmatrix i.e. flows from sources
#' in rows to the sinks in columns.
#' @references LIM package reference
#' @param readLIM (required) A linear inverse model read into the R environment with the R-package LIM.
#' Can be achieved through \code{readLIM <- Read("path_to_input_file")}
#' @param web (optional) The solved (food) web problem, i.e. the values of the unknowns;
#' if not specified, the model is solved first, using Ldei.
#' @details The LIM must be set up in a specific way. sum columns must equal sum rows.
#' @return Returns an energy-flux model in the format needed for the getJacobian function.
#' @export
getFlowMatrix <- function(readLIM, web = NULL) {
  lim <- Setup(readLIM)
  flows <- readLIM$flows[,1:2]
  flowmatrix <- lim$Flowmatrix
  if(is.null(web)) {
    if(!is.null(lim$Cost) || !is.null(lim$Profit)) {
      web <- Linp(lim)$X
    } else {
      web <- Lsei(lim, parsimonious = TRUE)$X
    }
  }
  X <- as.vector(web)
  dup <- which(duplicated(flows))
  for(i in dup) {
    a <- which(flows[,1] == flows[i, 1])
    b <- which(flows[,2] == flows[i, 2])
    same <- c(a, b)[duplicated(c(a, b))]
    X[i] <- sum(X[same])
  }
  Xpos <- pmax(0, X)
  ii <- which(flowmatrix > 0, arr.ind = TRUE)
  flowmatrix[ii] <- Xpos[lim$Flowmatrix[ii]]
  Xneg <- -1 * pmin(0, X)
  if(sum(Xneg) > 0) {
    flowmatrix[ii[, c(2, 1)]] <- flowmatrix[ii[, c(2, 1)]] +
      Xneg[lim$Flowmatrix[ii]]
  }
  return(flowmatrix)
}

#' Calculate value of LIM variables from flow solutions.
#'
#' This function calculates the value of variables as defined in the LIM from the flow solutions.
#' @references LIM package reference.
#' @param readLIM (required) A linear inverse model read into the R environment with the R-package LIM.
#' Can be achieved through \code{readLIM <- Read("path_to_input_file")}
#' @param web (optional) The solved (food) web problem, i.e. the values of the unknowns;
#' if not specified, the model is solved first, using Ldei.
#' @details The optimal solution of the LIM only provides the flow solutions.
#' However, in the LIM variables can be defined as the sum of some values.
#' This function calculates the actual value of these variables based on the flow solutions.
#' @return Returns a named vector with all variable values.
#' @export
getVariables <- function(readLIM, web = NULL) {
  lim <- Setup(readLIM)
  vars <- numeric(lim$NVariables)
  pars <- readLIM$pars$val
  vareq <- readLIM$vars
  if(is.null(web)) {
    if(!is.null(lim$Cost) || !is.null(lim$Profit)) {
      web <- Linp(lim)$X
    } else {
      web <- Lsei(lim, parsimonious = TRUE)$X
    }
  }
  # Possible to remove loop?
  for (i in 1:lim$NVariables) {
    tempvars <- vars
    subset <- vareq[vareq$nr == i,]
    sum <-
      sum(pars[subset$par1]*subset$val, na.rm = TRUE) +
      sum(pars[subset$par2]*subset$val, na.rm = TRUE) +
      sum(pars[subset$par3]*subset$val, na.rm = TRUE) +
      sum(pars[subset$par4]*subset$val, na.rm = TRUE) +
      sum(tempvars[subset$var]*subset$val, na.rm = TRUE) +
      sum(web[subset$flow]*subset$val, na.rm = TRUE)
    vars[i] = sum
  }
  names(vars) <- lim$Variables
  return(vars)
}
s
getTag <- function(vars, tag) {
  x <- vars[which(grepl(tag, names(vars)))]
  names(x) <- gsub(tag, "", names(x))
  return(x)
}

#' Get conversion efficiencies from a LIM.
#'
#' This function calculates assimilation and growth efficiences from a given LIM.
#' @references LIM package reference
#' @param FM (required) A flow matrix with flows from source in rows to sink in columns.
#' @param vars (required) A named vector with the values of variables defined in the LIM.
#' @param lim (required) Setup(Read(lim.input))
#' @param aTag (optional) Tag assigned to the variables containing the assimilated part. Default is set to
#' "ass". Not case sensitive.
#' @param gTag (optional) Tag assigned to the variables containing the growth part. Default is set to
#' "growth". Not case sensitive.
#' @details The LIM must be set up in a specific way. The assimilation efficiency is calculated as
#' the assimilated part (which is defined as variable in the LIM and calculated in the
#' function getVariables) divided by the ingestion (which is the sum of the organism's
#' column in the Flow Matrix FM).
#' by dividing the amount of assimilated material/energy
#' (must be defined in the original LIM model as variable) by the total ingestion of the organism.
#' This function calculates growth (or secondary production) efficiences by dividing the amount
#' of growth by the amount of assimilated material/energy
#' (both must be defined in the original LIM model as variable).
#' The LIM must be set up in a specific way. The assimilation efficiency is calculated as
#' the assimilated part (which is defined as variable in the LIM and calculated in the
#' function getVariables) divided by the ingestion (which is the sum of the organism's
#' column in the Flow Matrix FM).
#' !!! Growth efficiency BAC = growth / inflow
#' @return Returns a named vector with assimilation efficiencies.
#' @export
getCE <- function(FM, vars, lim, aTag = "ass", gTag = "growth") {
  # TODO: how to deal with bacteria
  # TODO: how to deal with words containing more than comp name and tag.
  names(vars) <- toupper(names(vars))
  AE <- rep(NA, length = lim$NComponents)
  names(AE) <- lim$Components$name
  GE <- rep(NA, length = lim$NComponents)
  names(GE) <- lim$Components$name
  AP <- getTag(vars, toupper(aTag))
  GP <- getTag(vars, toupper(gTag))
  AE[names(AP)] <- AP / colSums(FM[,names(AP)], na.rm = TRUE)
  GE[names(GP)] <- GP / AP
  CE <- list(AE = AE, GE = GE)
  return(CE)
}

# model$type
# model$LIM
# model$web
# model$aTag
# model$gTag
extractLIMdata() <- function(model) {
  FM <- getFlowMatrix(model$LIM, web = model$web)

  BM <- model$LIM$comp[,"val"]
  names(BM) <- model$LIM$comp[,"name"]

  CE <- getCE(
    FM = FM,
    vars = getVariables(model$LIM, model$web),
    lim = Setup(model$LIM),
    aTag = model$aTag,
    gTag = model$gTag
  )

  externals <- model$LIM$externnames

  #TODO dead
  #TODO MR

  return(list(
    FM = FM, BM = BM, AE = CE$AE, GE = CE$GE, externals = externals
  ))
}
