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
  # TODO: refactor, DRY
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
