#' Get a Flowmatrix from a LIM
#'
#' This function extracts flow data from the LIM and produces a Flowmatrix i.e. flows from sources
#' in rows to the sinks in columns.
#' @references LIM package reference
#' @param lim (required) A lim model.
#' @details The LIM must be set up in a specific way.
#' @return Returns an energy-flux model in the format needed for the getJacobian function.
#' @export
getFlowMatrix <- function (originalFM = NULL, flows = NULL, answers = NULL) {
  # Check function input
  if(is.null(originalFM)) {stop("Supply flowmatrix with flownumbers from LIM$Flowmatrix")}
  if(is.null(flows)) {stop("Supply the flows dataframe from readLIM$flows")}
  if(is.null(answers)) {stop("Supply a vector with all flow values")}

  # Create a vector which will store all total flow values
  totalflows <- numeric(length(flows[,1]))
  # Create a vector with all flownumbers
  flownrs <- c(1:length(flows[,1]))

  for(flownr in flownrs) {
    same <- findParallelFlows(flows = flows, flownr = flownr)
    totalflow <- sum(answers[same], na.rm = TRUE) # sum the flow values of all flows with the same source and sink.
    totalflows[flownr] <- totalflow # store the value at the index of the original flow
  }

  # Find the locations in the matrix with flownumbers
  indices <- match(flownrs, originalFM)
  indices <- indices[!is.na(indices)] # remove NA's (no matches) from the list.
  # Replace the flownumbers in the original matrix with the corresponding
  # total flow values (so if the flownr is 3 in the matrix, but it should be flow
  # 2 and 3, it will insert the totalflow of flow 2 and 3).
  FM <- replace(originalFM, indices, totalflows[originalFM[indices]])
  # Return a flowmatrix with flow values.
  return(FM)
}

function (lim, web = NULL)
{
  flowmatrix <- lim$Flowmatrix
  if (is.null(web))
    if (!is.null(lim$Cost) || !is.null(lim$Profit))
      web <- Linp(lim)$X
    else web <- Lsei.lim(lim, parsimonious = TRUE)$X
    X <- as.vector(web)
    Xpos <- pmax(0, X)
    ii <- which(flowmatrix > 0, arr.ind = TRUE)
    flowmatrix[ii] <- Xpos[lim$Flowmatrix[ii]]
    Xneg <- -1 * pmin(0, X)
    if (sum(Xneg) > 0)
      flowmatrix[ii[, c(2, 1)]] <- flowmatrix[ii[, c(2, 1)]] +
      Xneg[lim$Flowmatrix[ii]]
    return(flowmatrix)
}
<environment: namespace:LIM>
