# Detritus, DF, P, O
# Use netto FM

set.seed(1994)
fwnames <- c("LAB", "REF", "meiDF", "macDF", "meiP", "macP", "meiO", "macO")
AM <- matrix(c(
  0, 0, 1, 1, 0, 0, 1, 1,
  0, 0, 1, 0, 0, 0, 0, 0,
  0, 1, 0, 0, 1, 1, 1, 1,
  0, 1, 0, 0, 0, 1, 0, 1,
  1, 1, 0, 0, 1, 1, 1, 1,
  1, 1, 0, 0, 0, 1, 0, 1,
  1, 1, 0, 0, 1, 1, 0, 1,
  1, 1, 0, 0, 0, 1, 0, 1
), byrow = T, ncol = length(fwnames), nrow = length(fwnames))
rownames(AM) <- fwnames ; colnames(AM) <- fwnames

FM <- AM
FM[which(FM == 1)] <- runif(length(which(FM == 1)), min = 0.01, max = 10)
BM <- runif(dim(FM)[1], min = 1, max = 100)
AE <- runif(dim(FM)[1], min = 0.1, max = 1)
GE <- runif(dim(FM)[1], min = 0.1, max = 1)
names(BM) <- fwnames ; names(AE) <- fwnames ; names(GE) <- fwnames

model <- list(
  type = "EF",
  FM = FM,
  BM = BM,
  AE = AE,
  GE = GE,
  dead = list(names = fwnames[1:2], def = c("Def", "Def")),
  netto = T
)

dlab <- FM[,1] / rowSums(FM[,1:2]) ; dlab <- dlab[3:length(fwnames)]
rlab <- FM[,2] / rowSums(FM[,1:2]) ; rlab <- rlab[3:length(fwnames)]

## Expected answer
FMn <- FM - t(FM)
FMn[which(FMn < 0)] <- 0
FMn[1:2,] <- FM[1:2,] ; FMn[,1:2] <- FM[,1:2]
JM <- matrix(c(
  # LAB
  0,
  (FM[2,1] - FM[1,2] + sum(FM[2,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * dlab) ) / BM[2],
  (FM[3,1] - FM[1,3] + sum(FM[3,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * dlab) ) / BM[3],
  (FM[4,1] - FM[1,4] + sum(FM[4,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * dlab) ) / BM[4],
  (FM[5,1] - FM[1,5] + sum(FM[5,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * dlab) ) / BM[5],
  (FM[6,1] - FM[1,6] + sum(FM[6,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * dlab) ) / BM[6],
  (FM[7,1] - FM[1,7] + sum(FM[7,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * dlab) ) / BM[7],
  (FM[8,1] - FM[1,8] + sum(FM[8,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * dlab) ) / BM[8],

  # REF
  (FM[1,2] - FM[2,1] + sum(FM[1,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rlab) ) / BM[1],
  0,
  (FM[3,2] - FM[2,3] + sum(FM[3,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rlab) ) / BM[3],
  (FM[4,2] - FM[2,4] + sum(FM[4,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rlab) ) / BM[4],
  (FM[5,2] - FM[2,5] + sum(FM[5,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rlab) ) / BM[5],
  (FM[6,2] - FM[2,6] + sum(FM[6,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rlab) ) / BM[6],
  (FM[7,2] - FM[2,7] + sum(FM[7,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rlab) ) / BM[7],
  (FM[8,2] - FM[2,8] + sum(FM[8,3:length(fwnames)] * (1-AE[3:length(fwnames)]) * rlab) ) / BM[8],

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
  -FMn[5,7] / BM[7],
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
  AE[7] * GE[7] * FMn[5,7] / BM[5],
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
