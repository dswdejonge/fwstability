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

getLoopWeight <- function(IS, d, k) {
  LW <- abs(prod(IS) / prod(d)) ^ (1/k)
  return(LW)
}

getFeedback <- function(IS) {
  fdb <- prod(IS)
  return(fdb)
}

assessFeedback <- function(JM, findLoops = T,
                           loopFile = "allLoops.txt",
                           assess = c("fdb", "mlw")) {

  # find cycles if findLoops = T

  # read cycle data
  allLoops <- readLines(loopFile)
  allLoops <- strsplit(allLoops, " ")
  allLoops <- lapply(allLoops, function(x) as.numeric(x))

  # create storage matrix with compname -> compname -> compname fdw mlw
  result <- matrix(nrow = dim(allLoops)[1], ncol = length(assess))
  colnames(result) <- assess

  # find fdb and mlw
  for(i in 1:dim(allLoops)[1]){
    pathway <- allLoops[i,]
    coords <- matrix(
      data = c(pathway[1:length(pathway)-1], pathway[2:length(pathway)]),
      ncol = 2, nrow = length(pathway)-1)
    IS <- JM[coords]
    result[i, "fdb"] <- getFeedback(IS)
    result[i, "mlw"] <- getLoopWeight(IS, d, k)
  }
  allLoops <- cbind(allLoops, result)
  return(allLoops)
}
