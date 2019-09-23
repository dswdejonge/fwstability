# Check data format of JM
checkJMformat <- function(JM) {
  if(is.null(colnames(JM)) | is.null(rownames(JM))) {
    stop("Jacobian matrix must have named rows and columns.")
  } else if(!all(colnames(JM) == rownames(JM))) {
    stop("Jacobian matrix must have same names in rows and columns.")
  }
}

#' Assess stabiliy effect of each compartment.
#'
#' This function assesses the affect on stability of removing one compartment from
#' the Jacobian matrix.
#' @param JM (required) A Jacobian matrix.
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
#' is true if delta is positive.
#' @export
assessComps <- function(JM) {
  checkJMformat(JM)
  ini_s <- getStability(JM)
  df <- data.frame(
    removed_comp = rownames(JM),
    ini_s = ini_s,
    new_s = NA)

  for(i in 1:dim(JM)[1]){
    m <- JM[-i,-i]
    df$new_s[i] <- getStability(m)
  }
  df$delta <- df$new_s - df$ini_s
  df$delta_frac <- df$delta / df$ini_s
  return(df)
}

# if delta is negative, the system becomes more stable.
# if delta is positive, the system becomes less stable.
# show in vignette that doubling small interaction strength can have larger
# effect on stability than doubling large interaction strengths, thus a
# function + 10 would be more difficult to interpret in stead of relative change.
# per link (not per link pair) also intraspecific interaction.

#' Assess stabiliy effect of a fixed link alteration
#'
#' This function assesses the effect on stability of altering each interaction strength
#' in a Jacobian matrix according to a fixed function.
#' @param JM (required) A Jacobian matrix.
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
#' is true if delta is positive.
#' @export
assessLinksFixed <- function(JM,
                        func = function(x) {return(x * 2)}) {
  checkJMformat(JM)
  eg <- expand.grid(rownames(JM), colnames(JM))
  ini_s <- getStability(JM)
  df <- data.frame(
    eff.of = eg$Var1,
    eff.on = eg$Var2,
    ini_s = ini_s,
    new_s = NA)
  for(i in 1:length(JM)) {
    m <- JM
    m[i] <- func(JM[i])
    df$new_s[i] <- getStability(m)
  }
  df$delta <- df$new_s - df$ini_s
  df$delta_frac <- df$delta / df$ini_s

  return(df)
}

# Vary the values of one pair between 0 and Xaij a number of times,
# and calculate the probability that the food web becomes unstable.
# Method in deruiter1995: 0 - 2ij, 1% below threshold
# Set the diagonal so that the mean is 1% below the
# critical values (so that it is just stable).
# critical diagonal found by substracting I*maxEV from the original matrix
# per pair. no intraspecific interaction.
# allow choice of min and max (now 0 to 2ij)
assessLinksPerm <- function(JM,
                            perms = 100, threshold = 0.01, seed = 1) {
  set.seed(seed)
  checkJMformat(JM)
  pairs <- which(lower.tri(JM), arr.ind = TRUE)

  cJM <- JM
  diag(cJM) <- diag(cJM) - getStability(JM)
  diag(cJM) <- diag(cJM) - abs(diag(cJM)) * threshold

  if(getStability(cJM) > 0) {
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
      if(getStability(m) > 0) {
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

# wrapper function assesslinks
# method = "default"
# method = "permutations"
assessLinks <- function(JM, method = "default") {
  if(method == "default") {
    df <- assessLinksFixed(JM)
  } else if(method == "permutations") {
    df <- assessLinksPerm(JM)
  } else {
    stop("Unknown method.")
  }
  return(df)
}

# i = j works
fluxSizeDiversity <- function(FM){
  FMsum <- sum(FM)
  m <- FM / FMsum * log(FM / FMsum)
  H <- -sum(m, na.rm = TRUE)
  return(H)
}

# i cannot equal j
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
#' \item{Ulanowicz 1997.}
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
#' This function calculates weighted connectance, capturing the skewness of flows.
#' @param FM (required) A flow matrix with flows from source in rows to sink in columns.
#' @references \itemize{
#' \item{
#' van Altena, C., Hemerik, L., de Ruiter, P.C., 2016. Food web stability and weighted
#' connectance: the complexity-stability debate revisited. Theor. Ecol. 9, 49–58.
#' https://doi.org/10.1007/s12080-015-0291-7
#' }
#' }
#' @details Cw can vary between 0 and 1. Larger values of Cw indicate a more even distribution
#' of fluxes.
#' @return Returns a double.
#' @export
getCw <- function(FM) {
  m <- getConnPerNode(FM)
  Cw <- m / length(FM[,1])
  #Cw <- exp((fluxSizeDiversity(FM) - averageMutualInfo(FM)) / 2) / length(FM[,1])
  if(Cw < 0 | Cw > 1) {
    stop("Something is wrong. Cw should vary within 0 to 1.")
  }
  return(Cw)
}

#' Get loop weight.
#'
#' This function calculates the loop weight.
#' @param IS (required) Numeric vector of length \code{k} with interaction strengths.
#' @param d (optional) Natural death rates in same unit as IS. Required to
#' replicate method cited in reference.
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
#' @param IS (required) Numeric vector of length \code{k} with interaction strengths.
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

# MR should be in t-1 (same unit as JM)
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
    if(verbose) {message(paste0("Finding all ",maxL," loops in network..."))}
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
    result$mlw[i] <- getLoopWeight(IS, d = MR)
  }
  if(verbose) {message("Done.")}
  return(result)
}
