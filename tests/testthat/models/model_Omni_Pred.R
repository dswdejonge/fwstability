# Detritus, DF, P, O
# Meiofauna have no feedback to labile detritus.
# Defecation is split over lab and ref, mortality always to ref.
# Use netto FM
consumptionRate <- function(MR, BM, P, AE, GE) {
  F <- (MR * BM + P) / (AE * GE)
  return(F)
}

topDownBalancing <- function(PM, MR, BM, AE, GE){
  FM <- PM
  finished <- numeric()
  N <- which(rowSums(PM) == 0)
  while(length(N) > 0) {
    for(n in N) {
      Q <- consumptionRate(MR[n], BM[n], sum(FM[n,], na.rm = T), AE[n], GE[n])
      FM[,n] <- PM[,n]*BM / sum(PM[,n]*BM) * Q
    }
    finished <- c(finished, N)
    left <- PM[,-finished]
    if(is.null(dim(left))) {break}
    N <- which(rowSums(left) == 0)
    N <- N[!N %in% finished]
  }
  return(FM)
}

fwnames <- c("LAB", "REF", "meiDF", "macDF", "meiP", "macP", "meiO", "macO")
AM <- matrix(c(
  0, 0, 1, 1, 0, 0, 1, 1,
  0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 1, 0, 1,
  0, 0, 0, 0, 0, 1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 0, 1,
  0, 0, 0, 0, 0, 1, 0, 0
), byrow = T, ncol = length(fwnames), nrow= length(fwnames))
rownames(AM) <- fwnames ; colnames(AM) <- fwnames

MR <- c(NA, NA, 1, 0.5, 1, 0.5, 1, 0.5)
BM <- c(50, 100, 30, 25, 20, 20, 35, 10)
AE <- c(NA, NA, rep(0.9, 6))
GE <- c(NA, NA, rep(0.9, 6))
names(BM) <- fwnames ; names(AE) <- fwnames ; names(GE) <- fwnames ; names(MR) <- fwnames
# topDownBalancing needs apex predator
FM <- topDownBalancing(PM = AM, MR = MR, BM = BM, AE = AE, GE = GE)
# Detritus feedback
defecation  <- colSums(FM[,3:length(fwnames)])*(1-AE[3:length(fwnames)])
mortalities <- colSums(FM[,3:length(fwnames)])*
  AE[3:length(fwnames)]*GE[3:length(fwnames)] - rowSums(FM[3:length(fwnames),], na.rm = T)
ldet <- c(0, 0, 0.3, 0.4, 0.3, 0.5)
rdet <- 1-ldet
FM[3:8, 1] <- defecation * ldet
FM[3:8, 2] <- defecation * rdet + mortalities
#FM[1:2, 1:2] <- 0
# Fraction defecation matrix
frac <- matrix(
  c(0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, defecation*rdet/FM[3:8,2],
    rep(0,8*6)),
  ncol = length(fwnames), nrow = length(fwnames))
rownames(frac) <- fwnames; colnames(frac) <- fwnames

model <- list(
  type = "EF",
  FM = FM,
  BM = BM,
  AE = AE,
  GE = GE,
  dead = list(names = fwnames[1:2], frac = frac),
  netto = T
)

## Expected answer
FMn <- FM - t(FM)
FMn[which(FMn < 0)] <- 0
FMn[1:2,] <- FM[1:2,] ; FMn[,1:2] <- FM[,1:2]
JM <- matrix(c(
  # LAB
  0,
  (sum(FM[2,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * ldet, na.rm = T) ) / BM[2],
  (FM[3,1] - FM[1,3] + sum(FM[3,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * ldet, na.rm = T) ) / BM[3],
  (FM[4,1] - FM[1,4] + sum(FM[4,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * ldet, na.rm = T) ) / BM[4],
  (FM[5,1] - FM[1,5] + sum(FM[5,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * ldet, na.rm = T) ) / BM[5],
  (FM[6,1] - FM[1,6] + sum(FM[6,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * ldet, na.rm = T) ) / BM[6],
  (FM[7,1] - FM[1,7] + sum(FM[7,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * ldet, na.rm = T) ) / BM[7],
  (FM[8,1] - FM[1,8] + sum(FM[8,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * ldet, na.rm = T) ) / BM[8],

  # REF
  (sum(FM[1,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rdet, na.rm = T) ) / BM[1],
  0,
  (FM[3,2] - FM[2,3] + sum(FM[3,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rdet, na.rm = T) ) / BM[3],
  (FM[4,2] - FM[2,4] + sum(FM[4,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rdet, na.rm = T) ) / BM[4],
  (FM[5,2] - FM[2,5] + sum(FM[5,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rdet, na.rm = T) ) / BM[5],
  (FM[6,2] - FM[2,6] + sum(FM[6,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rdet, na.rm = T) ) / BM[6],
  (FM[7,2] - FM[2,7] + sum(FM[7,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rdet, na.rm = T) ) / BM[7],
  (FM[8,2] - FM[2,8] + sum(FM[8,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rdet, na.rm = T) ) / BM[8],

  # meiDF
  AE[3] * GE[3] * FM[1,3] / BM[1],
  AE[3] * GE[3] * FM[2,3] / BM[2],
  0,
  0,
  -FMn[3,5] / BM[5],
  -FMn[3,6] / BM[6],
  -FMn[3,7] / BM[7],
  -FMn[3,8] / BM[8],

  # macDF
  AE[4] * GE[4] * FM[1,4] / BM[1],
  AE[4] * GE[4] * FM[2,4] / BM[2],
  0,
  0,
  -FMn[4,5] / BM[5],
  -FMn[4,6] / BM[6],
  -FMn[4,7] / BM[7],
  -FMn[4,8] / BM[8],

  # meiP
  0,
  0,
  AE[5] * GE[5] * FMn[3,5] / BM[3],
  0,
  0,
  -FMn[5,6] / BM[6],
  AE[5] * GE[5] * FMn[7,5] / BM[7],
  -FMn[5,8] / BM[8],

  # macP
  0,
  0,
  AE[6] * GE[6] * FMn[3,6] / BM[3],
  AE[6] * GE[6] * FMn[4,6] / BM[4],
  AE[6] * GE[6] * FMn[5,6] / BM[5],
  0,
  AE[6] * GE[6] * FMn[7,6] / BM[7],
  AE[6] * GE[6] * FMn[8,6] / BM[8],

  # meiO
  AE[7] * GE[7] * FM[1,7] / BM[1],
  AE[7] * GE[7] * FM[2,7] / BM[2],
  AE[7] * GE[7] * FMn[3,7] / BM[3],
  0,
  -FMn[7,5] / BM[5],
  -FMn[7,6] / BM[6],
  0,
  -FMn[7,8] / BM[8],

  # macO
  AE[8] * GE[8] * FM[1,8] / BM[1],
  AE[8] * GE[8] * FM[2,8] / BM[2],
  AE[8] * GE[8] * FMn[3,8] / BM[3],
  AE[8] * GE[8] * FMn[4,8] / BM[4],
  AE[8] * GE[8] * FMn[5,8] / BM[5],
  -FMn[8,6] / BM[6],
  AE[8] * GE[8] * FMn[7,8] / BM[7],
  0
), nrow = dim(AM)[1], ncol = dim(AM)[2])
rownames(JM) <- fwnames ; colnames(JM) <- fwnames
