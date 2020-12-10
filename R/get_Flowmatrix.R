# The function in this file is an adapted version of the 'Flowmatrix' function from the LIM package
# written by Soetaert & van Oevelen (2015).
# Copyright: Karline Soetaert, Dick van Oevelen
# License: GPL (≥ 2)]

#' Get Flowmatrix from LIM
#'
#' This function extracts flow data from the LIM and produces a Flowmatrix i.e. flows from sources
#' in rows to the sinks in columns.
#' @references \code{\link{LIM}} package, Soetaert & van Oevelen 2015.
#' @param readLIM (required) A linear inverse model read into the R environment with the R-package LIM.
#' Can be achieved with \code{readLIM <- Read("path_to_input_file")}
#' @param web (optional) The solved (food) web problem, i.e. the values of the unknowns;
#' if not specified, the model is solved first, using Ldei.
#' @param lim (optional) The set-up linear problem.
#' Can be achieved with \code{lim <- Setup(readLIM)}.
#' @param verbose (optional) Should the function provide messages?
#' @details This function is an adapted version of the \code{Flowmatrix} function from the LIM package.
#' The major difference is that this function also provides the right answer if multiple parallel flows occur
#' between the same two compartments.
#' @return Returns a named numeric matrix.
#' @export
getFlowMatrix <- function(readLIM, web = NULL, lim = NULL, verbose = T) {
  if(is.null(lim)) {
    lim <- LIM::Setup(readLIM)
  }
  flows <- readLIM$flows[,1:2]
  flowmatrix <- lim$Flowmatrix
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
  X <- as.vector(web)
  X2 <- X
  dup <- which(duplicated(flows))
  for(i in dup) {
    a <- which(flows[,1] == flows[i, 1])
    b <- which(flows[,2] == flows[i, 2])
    same <- c(a, b)[duplicated(c(a, b))]
    X[i] <- sum(X2[same])
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
