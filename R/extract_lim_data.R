getTag <- function(vars, tag) {
  if(!is.atomic(vars) | is.null(names(vars))) {
    stop("getTag only accepts named vectors for argument \"vars\"")
  } else if(!is.character(tag) | length(tag) > 1) {
    stop("getTag only accepts a string for argument \"tag\"")
  }
  names(vars) <- toupper(names(vars))
  tag <- toupper(tag)
  x <- vars[which(grepl(tag, names(vars)))]
  names(x) <- gsub(tag, "", names(x))
  if(length(x) == 0) {
    warning(paste("Provided or default tag",tag,"not found in model."))
  }
  return(x)
}

#' Calculate value of LIM variables from flow solutions.
#'
#' This function calculates the value of variables as defined in the LIM from the flow solutions.
#' @references \code{\link{LIM}} package, Soetaert & van Oevelen 2015.
#' @param readLIM (required) A linear inverse model read into the R environment with the R-package LIM.
#' Can be achieved with \code{readLIM <- Read("path_to_input_file")}
#' @param web (optional) The solved (food) web problem, i.e. the values of the unknowns;
#' if not specified, the model is solved first, using Lsei.
#' @param lim (option) A LIM that has been set up with the function (Setup(readLIM)).
#' @param verbose (optional) Should the function provide messages?
#' @details Variables in the LIM can be defined as the sum of flows, other variables and parameters.
#' This function calculates the value of these variables based on the flow solutions and LIM parameters.
#' @return Returns a named vector with all variables.
#' @export
getVariables <- function(readLIM, web = NULL, lim = NULL, verbose = TRUE) {
  if(is.null(lim)){
    lim <- LIM::Setup(readLIM)
  }
  vars <- numeric(lim$NVariables)
  pars <- readLIM$pars$val
  vareq <- readLIM$vars
  if(is.null(web)) {
    if(verbose) {
      message("fwstab: No model solutions given, LIM resolved by minimizing sum of squares.")
    }
    if(!is.null(lim$Cost) || !is.null(lim$Profit)) {
      web <- LIM::Linp(lim)$X
    } else {
      web <- LIM::Lsei(lim, parsimonious = TRUE)$X
    }
  }
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

#' Get conversion efficiencies from a LIM.
#'
#' This function calculates assimilation and growth efficiencies from a given LIM.
#' @references \code{\link{LIM}} package, Soetaert & van Oevelen 2015.
#' @param FM (required) A flow matrix with flows from source in rows to sink in columns.
#' @param vars (required) A named vector with the values of variables defined in the LIM.
#' @param lim (required) Setup(Read(lim.input))
#' @param aTag (optional) Tag assigned to the variables representing the assimilated amount of material.
#' Default is set to "ass". Not case sensitive.
#' @param gTag (optional) Tag assigned to the variables representing growth.
#' Default is set to "growth". Not case sensitive.
#' @param verbose (optional) Should the function provide messages?
#' @details In order for this function to work, the LIM must be set-up in a specific way.
#' The variables representing total assimilation and growth of an organism must be named as:
#' "compartmentTag" or "tagCompartment" (not case sensitive).
#' All non-dead organisms must have a variable with assimilation and growth.
#' For more information on setting up the LIM, please review the vignette.
#'
#' The assimilation efficiency is calculated as amount of assimilated material divided by the
#' total ingestion. Growth (or secondary production) efficiency is calculated by dividing
#' growth by the amount of assimilated material.
#' @return Returns a list with element "AE" (named vector with assimilation efficiencies) and element
#' "GE" (named vector with growth efficiencies).
#' @export
getCE <- function(FM, vars, lim, aTag = NULL, gTag = NULL, verbose = TRUE) {
  if(is.null(aTag)) {
    aTag <- "ass"
    if(verbose) {
      message("fwstab: Default tag \"ass\" is used to search model for assimilation.")
    }
  }
  if(is.null(gTag)) {
    gTag <- "growth"
    if(verbose) {
      message("fwstab: Default tag \"growth\" is used to search model for secondary production.")
    }
  }
  AE <- rep(NA, length = lim$NComponents)
  names(AE) <- toupper(lim$Components$name)
  GE <- rep(NA, length = lim$NComponents)
  names(GE) <- toupper(lim$Components$name)
  AP <- getTag(vars = vars, tag = aTag)
  GP <- getTag(vars = vars, tag = gTag)
  inAE <- names(AE) %in% names(AP)
  ii <- sort(names(AE)[inAE])
  AE[ii] <- AP[ii] / colSums(FM[,ii], na.rm = TRUE)
  inGE <- names(GE) %in% names(GP)
  ii <- sort(names(GE)[inGE])
  GE[ii] <- GP[ii] / AP[ii]
  CE <- list(AE = AE, GE = GE)
  return(CE)
}

#' Get mortality rates from a LIM.
#'
#' This function calculates the mortality rates (per unit time) of all organisms from a LIM.
#' @references \code{\link{LIM}} package, Soetaert & van Oevelen 2015.
#' @param BM (required) A named numeric vector containing biomasses of all LIM components.
#' @param web (required) A named numeric vector with the flow solutions.
#' @param vars (required) A named numeric vector with the LIM variables.
#' @param mTag (optional) Tag assigned to the flows or variables containing the natural, i.e.
#' non-predatory, mortality. Default is set to "mort". Not case sensitive.
#' @param verbose (optional) Should the function provide messages?
#' @details In order for this function to work, the LIM must be set-up in a specific way.
#' The flows or variables representing total natural mortality of an organism must be named as:
#' "compartmentTag" or "tagCompartment" (not case sensitive).
#' For more information on setting up the LIM, please review the vignette.
#' @return Returns a named vector with mortality rates (per unit time).
#' @seealso \code{getVariables}
#' @export
getMR <- function(BM, web, vars, mTag = NULL, verbose = TRUE) {
  if(is.null(mTag)) {
    mTag <- "mort"
    if(verbose) {
      message("fwstab: Default tag \"mort\" is used to search model for mortality.")
    }
  }
  names(BM) <- toupper(names(BM))
  MR <- rep(NA, length = length(BM))
  names(MR) <- names(BM)
  values <- c(web, vars)
  MP <- getTag(values, mTag)
  inMR <- names(MR) %in% names(MP)
  ii <- names(MR)[inMR]
  MR[ii] <- MP[ii] / BM[ii]
  return(MR)
}

#' Gets list with information on dead compartments.
#'
#' This function determines for each dead compartment whether or not defecation occurs into it, and,
#' if parallel flows occur, what fraction of each flow is defecation.
#' @references \code{\link{LIM}} package, Soetaert & van Oevelen 2015.
#' @param dead (required) A list with the element "names" that contains a character vector with
#' names of all dead compartments.
#' @param readLIM (required) A linear inverse model read into the R environment with the R-package LIM.
#' Can be achieved with \code{readLIM <- Read("path_to_input_file")}
#' @param web (required) A named numeric vector with the flow solutions.
#' @param FM (optional) A flow matrix with flows from source in rows to sink in columns.
#' @param defTag (optional) Tag assigned to the flows representing defecation.
#' Default is set to "def". Not case sensitive.
#' @param verbose (optional) Should the function provide messages?
#' @details In order for this function to work, the LIM must be set-up in a specific way.
#' The flows representing defecation of an organism must be named as:
#' "compartmentTag" or "tagCompartment" (not case sensitive).
#' For more information on setting up the LIM, please review the vignette.
#' @return Returns a list with two elements:
#' \itemize{
#' \item{names: a character vector with names of dead compartments}
#' \item{frac: a matrix the same size as FM containing the fraction of each flow that is defecation}
#' }
#' @seealso \code{getFlowMatrix}
#' @export
getDeadInfo <- function(dead, readLIM, web, FM = NULL, defTag = NULL, verbose = TRUE) {
  if(length(dead$names) == 0){return(NULL)}
  if(is.null(FM)) {
    FM <- getFlowMatrix(readLIM, web)
  }
  if(is.null(defTag)) {
    defTag <- "def"
    if(verbose) {
      message("fwstab: Default tag \"def\" used to search model for defecation.")
    }
  }

  DM <- matrix(0, nrow = length(readLIM$compnames), ncol = length(readLIM$compnames))
  rownames(DM) <- readLIM$compnames
  colnames(DM) <- readLIM$compnames
  # Set defecation flow to one for all flows with the def tag in the name
  all_def_flows <- readLIM$flows[grep(toupper(defTag), toupper(readLIM$flows$name)),]
  DM[cbind(all_def_flows$from, all_def_flows$to)] <- 1
  # Calculate fractions for flows with both mortality and defecation
  flows <- readLIM$flows[,1:2]
  dup <- which(duplicated(flows))
  for(i in dup) {
    if(TRUE %in% (flows[i,] < 0)) {
      next
    } else {
      a <- which(flows[,1] == flows[i, 1])
      b <- which(flows[,2] == flows[i, 2])
      same <- c(a, b)[duplicated(c(a, b))]

      DM[flows[i,"from"],flows[i,"to"]] <-
        sum(getTag(web[same], tag = defTag), na.rm = TRUE) /
        FM[flows[i,"from"],flows[i,"to"]]
    }
  }
  dead$frac <- DM

  return(dead)
}

#' Extract data from formatted LIM.
#'
#' This function automatically extracts (physiological) data from a LIM that is formatted
#' with tags.
#' @references \code{\link{LIM}} package, Soetaert & van Oevelen 2015.
#' @param model (required) A named list containing elements with food web model data.
#' \itemize{
#' \item \code{LIM} (required) should contain a read-in LIM. Achieved through \code{Read(<path-to-input-file>)}
#' \item \code{web} (optional) can contain a named numeric vector with flow solutions.
#' If the LIM is not resolved, the function will use the parsimonious (least distance) solution.
#' \item \code{aTag} (optional) Tag assigned to the variables representing the assimilated amount of material.
#' Default is set to "ass". Not case sensitive.
#' \item \code{gTag} (optional) Tag assigned to the variables representing growth.
#' Default is set to "growth". Not case sensitive.
#' \item \code{mTag} (optional) Tag assigned to the flows or variables containing the natural, i.e.
#' non-predatory, mortality. Default is set to "mort". Not case sensitive.
#' \item \code{deadTag} (optional) Tag assigned to the dead compartments (e.g. detritus, carrion, nutrients).
#' Default is set to "dead". Not case sensitive.
#' \item \code{defTag} (optional) Tag assigned to the flows representing defecation.
#' Default is set to "def". Not case sensitive.
#' \item \code{setup} (optional) Setup(Read(lim.input)). By inclusion the function won't setup the model,
#' thereby potentially saving time.
#' }
#' @param verbose (optional) Should the function provide messages?
#' @details In order for this function to work, the LIM must be set-up in a specific way.
#' The flows and variables representing assimilation, growth, defecation, and mortality of an organism
#' must be named as: "compartmentTag" or "tagCompartment" (not case sensitive).
#' The dead compartments must include the deadTag.
#' For more information on setting up the LIM, please review the vignette.
#' @return Returns a list with flowmatrix \code{FM}, biomasses \code{BM}, conversion efficiencies
#' \code{CE$AE and CE$GE}, names of external compartments \code{externals}, information on dead
#' compartments \code{dead}, and mortality rates \code{MR}.
#' @seealso \code{getJacobian}
#' @export
extractLIMdata <- function(model, verbose = TRUE) {
  FM <- getFlowMatrix(readLIM = model$LIM, web = model$web, lim = model$setup)

  # Stocks can be given directly (stored in comp) or in parameters (stored in params)
  BM <- model$setup$Components[,"val"]
  names(BM) <- toupper(model$setup$Components[,"name"])
  vars <- getVariables(model$LIM, model$web)

  CE <- getCE(
    FM = FM,
    vars = vars,
    lim = model$setup,
    aTag = model$aTag,
    gTag = model$gTag
  )

  externals <- model$LIM$externnames

  if(is.null(model$dead)) {
    fwnames <- toupper(model$LIM$compnames)
    if(is.null(model$deadTag)) {
      deadTag <- "dead"
      if(verbose) {
        message("fwstab: Default tag \"dead\" is used to search model for dead compartments.")
      }
    }
    model$dead <- list(names = c(fwnames[grepl(toupper(deadTag), fwnames)]))
  }
  dead <- getDeadInfo(
    dead = model$dead, readLIM = model$LIM,
    web = model$web, FM = FM, defTag = model$defTag)

  MR <- getMR(BM = BM, web = model$web, vars = vars, mTag = model$mTag)

  remove <- which(BM == 0)
  if(length(remove) > 0) {
    if(verbose) {
      message("fwstab: Internal components with biomass of zero are removed.")
    }
    FM <- FM[-remove, -remove]
    BM <- BM[-remove]
    CE$AE <- CE$AE[-remove]
    CE$GE <- CE$GE[-remove]
    MR <- MR[-remove]
    dead$frac <- dead$frac[-remove, -remove]
    l <- !(dead$names %in% names(BM)[remove])
    dead$names <- dead$names[l]
  }

  return(list(
    FM = FM, BM = BM, AE = CE$AE, GE = CE$GE,
    externals = externals, dead = dead, MR = MR)
  )
}
