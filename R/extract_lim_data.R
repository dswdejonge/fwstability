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
getFlowMatrix <- function(readLIM, web = NULL, lim = NULL) {
  if(is.null(lim)) {
    lim <- Setup(readLIM)
  }
  flows <- readLIM$flows[,1:2]
  flowmatrix <- lim$Flowmatrix
  if(is.null(web)) {
    message("fwstab: No model solutions given, LIM resolved by minimizing sum of squares.")
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
    message("fwstab: No model solutions given, LIM resolved by minimizing sum of squares.")
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
getCE <- function(FM, vars, lim, aTag = NULL, gTag = NULL) {
  if(is.null(aTag)) {
    aTag <- "ass"
    message("fwstab: Default tag \"ass\" is used to search model for assimilation.")}
  if(is.null(gTag)) {
    gTag <- "growth"
    message("fwstab: Default tag \"growth\" is used to search model for secondary production.")}
  # TODO: how to deal with words containing more than comp name and tag.
  AE <- rep(NA, length = lim$NComponents)
  names(AE) <- toupper(lim$Components$name)
  GE <- rep(NA, length = lim$NComponents)
  names(GE) <- toupper(lim$Components$name)
  AP <- getTag(vars = vars, tag = aTag)
  ### TEMPORARY CODE - start ###
  #AP["BAC"] <- sum(FM[,"BAC"], na.rm = TRUE)
  ### TEMPORARY CODE - stop ###
  GP <- getTag(vars = vars, tag = gTag)
  inAE <- names(AE) %in% names(AP)
  inAP <- names(AP) %in% names(AE)
  AE[inAE] <- AP[inAP] / colSums(FM[,names(AE)[inAE]], na.rm = TRUE)
  inGE <- names(GE) %in% names(GP)
  inGP <- names(GP) %in% names(GE)
  GE[inGE] <- GP[inGP] / AP[names(GE)[inGE]]
  ### TEMPORARY CODE - start ###
  #temp <- AE
  #temp <- GE[names(AE)]
  #GE <- temp
  ### TEMPORARY CODE - stop ###
  CE <- list(AE = AE, GE = GE)
  return(CE)
}

#' Get mortality rates from a LIM.
#'
#' This function calculates the mortality rates (per unit time) of all organisms from a LIM.
#' @references LIM package reference
#' @param BM (required) A named vector containing biomasses of all LIM components.
#' @param web (required) A named vector with the flow values.
#' @param mTag (optional) Tag assigned to the variables/flows containing the natural
#' flow mortality. Default is set to "mort". Not case sensitive.
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
#' @return Returns a named vector with mortality rates (per unit time).
#' @export
getMR <- function(BM, web, mTag = NULL) {
  if(is.null(mTag)) {
    mTag <- "mort"
    message("fwstab: Default tag \"mort\" is used to search model for mortality.")}
  names(BM) <- toupper(names(BM))
  MR <- rep(NA, length = length(BM))
  names(MR) <- names(BM)
  MP <- getTag(web, mTag)
  MR[names(MP)] <- MP / BM[names(MP)]
  return(MR)
}

getDeadInfo <- function(dead, readLIM, web, FM = NULL, defTag = NULL) {
  if(is.null(FM)) {
    FM <- getFlowMatrix(readLIM, web)
  }
  if(is.null(defTag)) {
    defTag <- "def"
    message("fwstab: Default tag \"def\" used to search model for defecation.")}

  dead <- adjustDeadInput(dead)

  dead$def <- rep("noDef", length(dead$names))
  allSinks <- readLIM$flows[,"to"]
  names(allSinks) <- readLIM$flows[,"name"]
  defComps <- getTag(allSinks, defTag)
  dead$def[dead$names %in% readLIM$compnames[defComps]] <- "Def"

  DM <- matrix(1, nrow = length(readLIM$compnames), ncol = length(readLIM$compnames))
  rownames(DM) <- readLIM$compnames
  colnames(DM) <- readLIM$compnames
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
        sum(getTag(web[same], tag = defTag), na.rm = T) /
        FM[flows[i,"from"],flows[i,"to"]]
    }
  }
  dead$frac <- DM

  return(dead)
}

# model$type
# model$LIM
# model$web)
# model$aTag
# model$gTag
extractLIMdata <- function(model) {
  FM <- getFlowMatrix(readLIM = model$LIM, web = model$web, lim = model$setup)

  BM <- model$LIM$comp[,"val"]
  names(BM) <- toupper(model$LIM$comp[,"name"])
  ### TEMPORARY CODE - start ###
  #if(model$site == "dist") {
  #  BM["DOCPHYTO_S"] <- 0.28
  #  BM["DOCOTHER_S"] <- 6.91
  #  BM["CARC"] <- sum(FM[,"CARC"])
  #} else if(model$site == "undist") {
  #  BM["DOCPHYTO_S"] <- 0.28
  #  BM["DOCPHYTO_S"] <- 4.46
  #  BM["CARC"] <- sum(FM[,"CARC"])
  #} else if(model$site == "ref") {
  #  BM["DOCPHYTO_S"] <- 0.65
  #  BM["DOCPHYTO_S"] <- 9.47
  #  BM["CARC"] <- sum(FM[,"CARC"])
  #} else {stop("Error! Unknown site.")}
  ### TEMPORARY CODE - stop ###

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
      message("fwstab: Default tag \"dead\" is used to search model for dead compartments.")}
    model$dead <- list(names = c(fwnames[grepl(toupper(deadTag), fwnames)]))
    ### TEMPORARY CODE - start ###
    #deadnames = c("PHYTO_S", "PHYTO_W", "SLAB_S", "SLAB_W", "REFRAC",
    #          "DOCPHYTO_S", "DOCOTHER_S", "CARC")
    #model$dead <- list(names = deadnames)
    ### TEMPORARY CODE - stop ###
  }
  dead <- getDeadInfo(
    dead = model$dead, readLIM = model$LIM,
    web = model$web, FM = FM, defTag = model$defTag)

  MR <- getMR(BM = BM, web = model$web, mTag = model$mTag)

  remove <- which(BM == 0)
  if(length(remove) > 0) {
    message("fwstab: Internal components with biomass of zero are removed.")
    FM <- FM[-remove, -remove]
    BM <- BM[-remove]
    CE$AE <- CE$AE[-remove]
    CE$GE <- CE$GE[-remove]
    MR <- MR[-remove]
    dead$frac <- dead$frac[-remove, -remove]
    l <- !(dead$names %in% names(BM)[remove])
    dead$names <- dead$names[l]
    dead$def <- dead$def[l]
  }

  return(list(
    FM = FM, BM = BM, AE = CE$AE, GE = CE$GE,
    externals = externals, dead = dead, MR = MR)
  )
}
