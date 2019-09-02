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


