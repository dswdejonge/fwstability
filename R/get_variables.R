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
#' @return Returns an energy-flux model in the format needed for the getJacobian function.
#' @export
getsVariables <- function(readLIM, web = NULL) {
  lim <- Setup(readLIM)
  vars <- numeric(lim$NVariables) # vals
  pars <- readLIM$pars$val # parvec
  vareq <- readLIM$vars
  if(is.null(web)) {
    if(!is.null(lim$Cost) || !is.null(lim$Profit)) {
      web <- Linp(lim)$X
    } else {
      web <- Lsei(lim, parsimonious = TRUE)$X
    }
  }
  # Possible to remove loop?
  for (i in 1:LIM$NVariables) {
    # refresh vector with variable parsimonious values
    varvec <- vals # vector with all variables values in right order: needs to be refreshed every time

    # get subset with the same equation nr
    subset <- vareq[vareq$nr == i,]

    # takes parameter, variable or flow nr from subset to use as index
    # to find the corresponding parsimonious values and
    # multiplies them by 'val' which is the coefficient (like 1 or -1) and sums them
    sum <-
      sum(parvec[subset$par1]*subset$val, na.rm = TRUE) +
      sum(parvec[subset$par2]*subset$val, na.rm = TRUE) +
      sum(parvec[subset$par3]*subset$val, na.rm = TRUE) +
      sum(parvec[subset$par4]*subset$val, na.rm = TRUE) +
      sum(varvec[subset$var]*subset$val, na.rm = TRUE) +
      sum(flowvec[subset$flow]*subset$val, na.rm = TRUE)

    # replaces 0 in data frame to actual sum value
    vals[i] = sum
  }
  vals <- data.frame(LIM$Variables, vals)
  return(vals)
}
