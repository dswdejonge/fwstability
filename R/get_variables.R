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
