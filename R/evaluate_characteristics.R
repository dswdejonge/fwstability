# Check data format of JM
checkJMformat <- function(JM) {
  if(is.null(colnames(JM)) | is.null(rownames(JM))) {
    stop("Jacobian matrix must have named rows and columns.")
  } else if(!all(colnames(JM) == rownames(JM))) {
    stop("Jacobian matrix must have same names in rows and columns.")
  }
}

# Assess absolute and relative change in stability if a compartment disappears.
# if delta is negative, the system becomes more stable.
# if delta is positive, the system becomes less stable.
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

# Assess the effect of altering each link on the overall stability.
# Default alteration is doubling the interaction strength.
# if delta is negative, the system becomes more stable.
# if delta is positive, the system becomes less stable.
# show in vignette that doubling small interaction strength can have larger
# effect on stability than doubling large interaction strengths, thus a
# function + 10 would be more difficult to interpret in stead of relative change.
# per link (not per link pair) also intraspecific interaction.
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

fluxSizeDiversity <- function(FM, cannibalism = FALSE){
  FMsum <- sum(FM)
  m <- FM / FMsum * log(FM / FMsum)
  if(cannibalism == FALSE){
    diag(m) <- 0
  }
  H <- -sum(m, na.rm = TRUE)
  return(H)
}

averageMutualInfo <- function(FM, cannibalism = FALSE){
  FMsum <- sum(FM)
  outgoingSum <- rowSums(FM)
  incomingSum <- colSums(FM)
  product <- expand.grid(outgoingSum, incomingSum)
  product$product <- product$Var1 * product$Var2
  productM <- matrix(product$product, nrow = dim(FM)[1], ncol = dim(FM)[2])
  m <- FM / FMsum * log(FM * FMsum / productM)
  if(cannibalism == FALSE){
    diag(m) <- 0
  }
  A <- sum(m, na.rm = TRUE)
  return(A)
}

#' Get weighted connectance.
#'
#' This function calculates weighted connectance, capturing the skewness of flows.
#' @param FM (required) A flow matrix with flows from source in rows to sink in columns.
#' @param cannibalism (optional) If cannibalism occurs in the food web.
#' @references \itemize{
#' \item{van Altena et al. 2016
#' }
#' }
#' @return Returns a double.
#' @export
getCw <- function(FM, cannibalism = FALSE) {
  c = cannibalism
  Cw <- exp((fluxSizeDiversity(FM, c) - averageMutualInfo(FM, c)) / 2) / length(FM[,1])
  return(Cw)
}
