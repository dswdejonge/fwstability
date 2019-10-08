#' Assess stability effect of each compartment.
#'
#' This function assesses the affect on stability of removing individual compartments from
#' the Jacobian matrix.
#' @inheritParams getStability
#' @references \itemize{
#' \item{
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced carbon cycles and
#' the absence of a relation between ecosystem complexity and stability. Ecol. Lett. 17,
#' 651–661. https://doi.org/10.1111/ele.12266
#' }
#' }
#' @return Returns a dataframe with absolute and relative changes in stability for removal
#' of each food web compartment from the Jacobian matrix.
#' @details If the change in stability (delta) is negative, the system becomes more stable if
#' the respective food web compartment is removed from the Jacobian matrix. The opposite
#' is true if delta is positive. \cr
#' The 'scalar' method might be somewhat slower than the default 'eigenvalue' method.
#' @export
assessComps <- function(JM, method = "eigenvalue",
                        mortalities = NULL, dead = NULL) {
  checkMformat(JM)
  ini_s <- getStability(JM, method, mortalities, dead)
  df <- data.frame(
    removed_comp = rownames(JM),
    ini_s = ini_s,
    new_s = NA)

  for(i in 1:dim(JM)[1]){
    m <- JM[-i,-i]
    df$new_s[i] <- getStability(m, method, mortalities, dead)
  }
  df$delta <- df$new_s - df$ini_s
  df$delta_frac <- df$delta / df$ini_s
  return(df)
}


#' Assess stabiliy effect of a fixed link alteration
#'
#' This function assesses the effect on stability of altering each interaction strength
#' in a Jacobian matrix according to a fixed function.
#' @inheritParams getStability
#' @param func (required) Function describing how to alter each interaction strength to
#' assess the effect on stability of the respective link. Default is doubling each
#' interaction strength.
#' @references \itemize{
#' \item{
#' de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of Interaction
#' Strengths, and Stability in Real Ecosystems. Science (80-. ). 269, 1257–1260.
#' https://doi.org/10.1126/science.269.5228.1257
#' }
#' }
#' @return Returns a dataframe with absolute and relative changes in stability for alteration
#' of each interaction strength in the Jacobian matrix.
#' @details If the change in stability (delta) is negative, the system becomes more stable if
#' the respective food web compartment is removed from the Jacobian matrix. The opposite
#' is true if delta is positive. \cr
#' Using the method 'scalar' might be somewhat slower than the method 'eigenvalue'.
#' @export
assessLinksFixed <- function(JM, method = "eigenvalue",
                             mortalities = NULL, dead = NULL,
                             func = function(x) {return(x * 2)}) {
  checkMformat(JM)
  eg <- expand.grid(rownames(JM), colnames(JM))
  ini_s <- getStability(JM, method, mortalities, dead)
  df <- data.frame(
    eff.of = eg$Var1,
    eff.on = eg$Var2,
    ini_s = ini_s,
    new_s = NA)
  for(i in 1:length(JM)) {
    m <- JM
    m[i] <- func(JM[i])
    df$new_s[i] <- getStability(m, method, mortalities, dead)
  }
  df$delta <- df$new_s - df$ini_s
  df$delta_frac <- df$delta / df$ini_s

  return(df)
}

#' Assess probability of destabilizing system per interaction
#'
#' This function determines the probability that a system becomes unstable if the
#' strengths in an interaction are varied.
#' @inheritParams getStability
#' @param perms (required) Default is 100. Number of times a pair of interaction strengths
#' is varied.
#' @param threshold (required) Default is 0.01. The Jacobian matrix is set at this fraction below
#' the stability threshold before starting permutations.
#' @param seed (required) Default is 1. Set seed to allow reproducability of results.
#' @references \itemize{
#' \item{
#' de Ruiter, P.C., Neutel, A.M., Moore, J.C., 1995. Energetics, Patterns of Interaction
#' Strengths, and Stability in Real Ecosystems. Science (80-. ). 269, 1257–1260.
#' https://doi.org/10.1126/science.269.5228.1257
#' }
#' }
#' @return Returns a dataframe with the probability (0 to 1) of destabilizing the system for each
#' pair of trophic interactions.
#' @details First, the diagonal of the Jacobian matrix is altered so that the matrix is
#' just stable. This is done by subtracting the maximum real part of the eigenvalues from the
#' diagonal values, and subsequently reducing the diagonal values by a fraction
#' (in argument \code{threshold}). The default threshold of 0.01 is used in De Ruiter et al. (1995).
#' Every pair of interaction strengths is then randomly varied between 0 and 2aij a certain
#' number of times. The probability that the food web becomes unstable is the total count of
#' unstable matrices divided by the total number of permutations. \cr
#' Using the method 'scalar' to asses stablity might be somewhat slower than the method 'eigenvalue'.
#' @seealso \code{getStability}
#' @export
assessLinksPerm <- function(JM, method = "eigenvalue",
                            mortalities = NULL, dead = NULL,
                            perms = 100, threshold = 0.01, seed = 1) {
  set.seed(seed)
  checkMformat(JM)
  pairs <- which(lower.tri(JM), arr.ind = TRUE)

  cJM <- JM
  diag(cJM) <- diag(cJM) - getStability(JM)
  diag(cJM) <- diag(cJM) - abs(diag(cJM)) * threshold

  if(getStability(cJM, method, mortalities, dead) > 0) {
    stop("Initial stability needs to be negative in order for this test to work.")
  }

  probs <- vector(length = length(pairs[,1]))
  for(i in 1:length(pairs[,1])){
    x <- pairs[i,1]; y <- pairs[i,2]
    counter <- 0
    for(j in 1:perms){
      m <- cJM
      m[x,y] <- runif(1, min=0, max=2*abs(m[x,y]))*sign(m[x,y])
      m[y,x] <- runif(1, min=0, max=2*abs(m[x,y]))*sign(m[y,x])
      if(getStability(m, method, mortalities, dead) > 0) {
        counter <- counter+1
      }
    }
    probs[i] <- counter / perms
  }
  df <- data.frame(
    comp1 = rownames(JM)[pairs[,1]],
    comp2 = colnames(JM)[pairs[,2]],
    Pinstability = probs)
  return(df)
}

#' Get flux size diversity.

#' Calculates the diversity in flux weights.
#' @param FM (required) A square flow matrix with flows from source in rows to sink in columns.
#' @references \itemize{
#' \item{
#' van Altena, C., Hemerik, L., de Ruiter, P.C., 2016. Food web stability and weighted
#' connectance: the complexity-stability debate revisited. Theor. Ecol. 9, 49–58.
#' https://doi.org/10.1007/s12080-015-0291-7
#' }
#' \item{
#' Boit, A.; Gaedke, U., 2014. Benchmarking successional progress in a quantitative food web.
#' Plos One 9(2):e90404. doi:10.1371/journal.pone.0090404
#' }
#' \item{Ulanowicz, R.E., 1997. Limitations on the connectivity of ecosystem flow networks. In:
#' Rinaldo, A.; Marani, A. (eds) Biological models. Instituto Veneto de Scienze, Lettre ed Arti,
#' Venica, pp 125-143.}
#' }
#' @details Cannibalistic flows (from i to i) are included in calculations.
#' @return Returns a double.
#' @export
fluxSizeDiversity <- function(FM){
  FMsum <- sum(FM)
  m <- FM / FMsum * log(FM / FMsum)
  H <- -sum(m, na.rm = TRUE)
  return(H)
}


#' Get averague mutual information.
#'
#' Calculates average mutual information.
#' @param FM (required) A square flow matrix with flows from source in rows to sink in columns.
#' @references \itemize{
#' \item{
#' van Altena, C., Hemerik, L., de Ruiter, P.C., 2016. Food web stability and weighted
#' connectance: the complexity-stability debate revisited. Theor. Ecol. 9, 49–58.
#' https://doi.org/10.1007/s12080-015-0291-7
#' }
#' \item{
#' Boit, A.; Gaedke, U., 2014. Benchmarking successional progress in a quantitative food web.
#' Plos One 9(2):e90404. doi:10.1371/journal.pone.0090404
#' }
#' \item{Ulanowicz, R.E., 1997. Limitations on the connectivity of ecosystem flow networks. In:
#' Rinaldo, A.; Marani, A. (eds) Biological models. Instituto Veneto de Scienze, Lettre ed Arti,
#' Venica, pp 125-143.}
#' }
#' @details Cannibalistic flows (from i to i) are not included in calculations.
#' @return Returns a double.
#' @export
averageMutualInfo <- function(FM){
  FMsum <- sum(FM)
  outgoingSum <- rowSums(FM)
  incomingSum <- colSums(FM)
  product <- expand.grid(outgoingSum, incomingSum)
  product$product <- product$Var1 * product$Var2
  productM <- matrix(product$product, nrow = dim(FM)[1], ncol = dim(FM)[2])
  m <- FM / FMsum * log(FM * FMsum / productM)
  diag(m) <- 0
  A <- sum(m, na.rm = TRUE)
  return(A)
}

#' Get effective or topological connectance per node.
#'
#' This function calculates effective or topological connectance per node, which can
#' plotted in the window of vitality.
#' @param FM (required) A flow matrix with flows from source in rows to sink in columns.
#' @param type \itemize{
#' \item{"eff" calculates the effective connectance per node (m).}
#' \item{"topo" calculates the topological connectance per node (m*), which is the special
#' case where all flows have equal weights.}
#' }
#' @references \itemize{
#' \item{
#' van Altena, C., Hemerik, L., de Ruiter, P.C., 2016. Food web stability and weighted
#' connectance: the complexity-stability debate revisited. Theor. Ecol. 9, 49–58.
#' https://doi.org/10.1007/s12080-015-0291-7
#' }
#' \item{Ulanowicz, R.E., 1997. Limitations on the connectivity of ecosystem flo networks. In:
#' Rinaldo, A.; Marani, A. (eds) Biological models. Instituto Veneto de Scienze, Lettre ed Arti,
#' Venica, pp 125-143.}
#' }
#' @return Returns a double.
#' @export
getConnPerNode <- function(FM, type = c("eff", "topo")) {
  if(length(type) == 1 && type == "topo"){
    FM[which(FM > 0)] <- 1
  }
  m <- exp((fluxSizeDiversity(FM) - averageMutualInfo(FM)) / 2)
  return(m)
}

#' Get weighted connectance.
#'
#' This function calculates weighted connectance, capturing system complexity and
#' the skewness of flows.
#' @param FM (required) A flow matrix with flows from source in rows to sink in columns.
#' @references \itemize{
#' \item{
#' van Altena, C., Hemerik, L., de Ruiter, P.C., 2016. Food web stability and weighted
#' connectance: the complexity-stability debate revisited. Theor. Ecol. 9, 49–58.
#' https://doi.org/10.1007/s12080-015-0291-7
#' }
#' \item{
#' Boit, A.; Gaedke, U., 2014. Benchmarking successional progress in a quantitative food web.
#' Plos One 9(2):e90404. doi:10.1371/journal.pone.0090404
#' }
#' }
#' @details Cw can vary between 1/number of trophic groups and 1. Larger values of Cw indicate a more even distribution
#' of fluxes. Cw is used as a measure of system complexity.
#' @return Returns a double.
#' @export
getCw <- function(FM) {
  m <- getConnPerNode(FM)
  Cw <- m / length(FM[,1])
  if(Cw < 1/dim(FM)[1] | Cw > 1) {
    stop("Something is wrong. Cw should lie between 1/S to 1, where S is number of network nodes.")
  }
  return(Cw)
}

#' Get loop weight.
#'
#' This function calculates the loop weight.
#' @param IS (required) Numeric vector of length \code{k} with interaction strengths.
#' @param d (optional) Numeric vector of length \code{k} with natural death rates in
#' the same unit as IS.
#' @references \itemize{
#' \item{
#' Neutel, A.M., Heesterbeek, J.A.P., Van De Koppel, J., Hoenderboom, G.,
#' Vos, A., Kaldeway, C., Berendse, F., De Ruiter, P.C., 2007. Reconciling
#' complexity with stability in naturally assembling food webs. Nature 449,
#' 599–602. https://doi.org/10.1038/nature06154
#' }
#' }
#' @details Loop weight is defined as the geometric mean of all absolute
#' interaction strengths in a loop of length k. It combines information on
#' the feedback of the loop, with natural death rates and the length of
#' the loop. The feedback is scaled by the natural death rates in the reference
#' because the stability measure is also scaled by death rates. Depending on
#' your stability measure you can decide to omit inclusion of death rates.
#' @return Returns a double.
#' @export
getLoopWeight <- function(IS, d = NULL) {
  k <- length(IS)
  if(!is.null(d) & !all(d > 0)) {
    d <- NULL
    warning("Deat rates should all be positive. Did not include death rate
            in calculation.")
  }
  if(is.null(d)) {
    LW <- abs(prod(IS)) ^ (1 / k)
  } else {
    LW <- abs(prod(IS) / prod(d)) ^ (1 / k)
  }
  return(LW)
}

#' Get loop feedback
#'
#' This function calculates the feedback of a loop.
#' @param IS (required) Numeric vector with interaction strengths in a loop.
#' @references \itemize{
#' \item{
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced
#' carbon cycles and the absence of a relation between ecosystem
#' complexity and stability. Ecol. Lett. 17, 651–661.
#' https://doi.org/10.1111/ele.12266
#' }
#' }
#' @details Feedback of a loop is the product of all interaction strengths
#' in a loop (so feedback is not additative but multiplicative).
#' @return Returns a double.
#' @export
getFeedback <- function(IS) {
  fdb <- prod(IS)
  return(fdb)
}

#' Maximum number of loops
#'
#' This function calculates the maximum number of loops, assuming full connectance.
#' @param n (required) Integer. Number of network compartments.
#' @param k (optional) Integer. Length of loop. Default NULL searches of loops of all length
#' (k=2 to k = n).
#' @details A fully connected network is assumed. The total number of loops of length k
#' can be found as n! / (n - k)!. If k is not given the total number of possible loops (i.e.
#' of length k = 2 to k = n) is found.
#' @return Returns the maximum number of possible loops.
#' @export
maxNrLoops <- function(n, k = NULL) {
  if(is.null(k)) {
    x <- factorial(n) * sum(1/factorial(0:(n-2)))
  } else {
    x <- factorial(n) / factorial(n - k)
  }
  return(x)
}

#' Assess feedback characteristics of network
#'
#' This function can find all loops or loops of a specific length, and consequently
#' calculate the feedback and loop weight of each loop.
#' @param JM (required) Jacobian matrix with interaction strengths.
#' @param MR (optional) Natural mortality/death rates for scaling the
#' maximum loop weight, same unit as Jacobian matrix.
#' @param compnames (optional) Vector with compartment names in same order
#' as the Jacobian matrix. If it is not included the names of JM
#' are used as compartment names. If the JM is not named, the output will
#' simply include the index of compartments to indicate loops.
#' @param findLoops (optional) Default is FALSE. You need a text file with the indices
#' of each loop per line. If you do not have this yet, this function redirects
#' to a recursive depth-first-search function to find all loops, or all loops
#' of length k (see function \code{dfs()}.
#' @param k (optional) Integer. Can be used if findLoops is TRUE,
#' and indicates that you only want to search for loops of length k.
#' @param file (optional) Default is "allLoops.txt". Is the text file
#' with indices of all loops. Should be present in the working directory.
#' @param verbose (optional) Default is TRUE.
#' Set to FALSE if you don't want messages printed.
#' @references \itemize{
#' \item{
#' Neutel, A.M., Heesterbeek, J.A.P., Van De Koppel, J., Hoenderboom, G.,
#' Vos, A., Kaldeway, C., Berendse, F., De Ruiter, P.C., 2007. Reconciling
#' complexity with stability in naturally assembling food webs. Nature 449,
#' 599–602. https://doi.org/10.1038/nature06154
#' }
#' \item{
#' Neutel, A.M., Thorne, M.A.S., 2014. Interaction strengths in balanced
#' carbon cycles and the absence of a relation between ecosystem
#' complexity and stability. Ecol. Lett. 17, 651–661.
#' https://doi.org/10.1111/ele.12266
#' }
#' }
#' @details Natural death/mortality rates (MR) can be found as the absolute
#' values on the diagonal.
#' Diagonal values for detritus are also regarded 'death' rates, as they
#' represent self-dampening effects. \cr
#'
#' The depth-first-search algorithm can find all loops in your network.
#' Please be aware that increasing the size of your network, exponentially
#' increases the computation time to find all loops. If you have a large
#' network you might want to limit your search to loops of length k = 2
#' or k = 3, as those are found to be the most important in determining
#' overall feedback and stability in your system (Neutel & Thorne 2014). \cr
#'
#' Feedback of a loop is the product of all interaction strengths
#' in a loop (so feedback is not additative but multiplicative).
#'
#' Loop weight is defined as the geometric mean of all absolute
#' interaction strengths in a loop of length k. It combines information on
#' the feedback of the loop, with natural death rates and the length of
#' the loop. The feedback is scaled by the natural death rates in the reference
#' because the stability measure is also scaled by death rates. Depending on
#' your stability measure you can decide to omit inclusion of death rates.
#' @return Returns a dataframe with all loops noted as
#' compName1->compName2->compNamek (column "loop"),
#' with feedbacks (column "fdb") and loop weights (column "lw") of
#' those loops.
#' @export
assessFeedback <- function(JM, MR = NULL, compnames = NULL,
                           findLoops = F, k = NULL,
                           file = NULL, verbose = T) {
  # path = getwd(), file = "allLoops.txt"
  if(findLoops) {
    maxL <- maxNrLoops(dim(JM)[1], k)
    if(verbose) {message(paste0(
      "Assuming full connectance, there are at most ",maxL," loops.\n
      Finding all loops in specified network..."))}
    AM <- abs(JM)
    AM[which(AM > 0)] <- 1
    # paste0(path,"/",file)
    file <- dfs(AM, k, output = file, verbose)
  }
  if(verbose) {message("Read loop data...")}
  allLoops <- readLines(paste0(file))
  allLoops <- strsplit(allLoops, " ")
  allLoops <- lapply(allLoops, function(x) as.numeric(x))

  if(verbose) {message("Find feedback and loop weight for all loops...")}
  if(is.null(compnames)) {
    loopnames <- unlist(lapply(allLoops, function(x) paste(x, collapse = "->")))
  } else {
    loopnames <- lapply(allLoops, function(x) compnames[x])
    loopnames <- unlist(lapply(loopnames, function(x) paste(x, collapse = "->")))
  }
  result <- data.frame(loop = loopnames, fdb = NA, lw = NA)
  for(i in 1:length(loopnames)){
    pathway <- allLoops[[i]]
    coords <- matrix(
      data = c(pathway[1:length(pathway)-1], pathway[2:length(pathway)]),
      ncol = 2, nrow = length(pathway)-1)
    IS <- JM[coords]
    result$fdb[i] <- getFeedback(IS)
    result$lw[i] <- getLoopWeight(IS, d = MR)
  }
  if(verbose) {message("Done.")}
  return(result)
}
