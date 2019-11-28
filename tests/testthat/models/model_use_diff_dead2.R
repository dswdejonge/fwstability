# Species deposit different material in different dead compartments
#
#   <----\<----\  mort i & def j
#  d --> i --> j
#  r <---/----/   def i & mort j

fwnames <- c("D", "R", "I", "J")
BM <- c(30, 30, 20, 10) ; names(BM) <- fwnames
AE <- c(NA, NA, 0.2, 0.3) ; names(AE) <- fwnames
GE <- c(NA, NA, 0.2, 0.3) ; names(GE) <- fwnames
MR <- c(NA, NA, 0.8, 0.7) ; names(MR) <- fwnames
# Top down balancing
Imort <- MR[3]*BM[3]
Jmort <- MR[4]*BM[4]
Fij <- (Jmort)/(AE[4]*GE[4])
Fdi <- (Imort+Fij)/(AE[3]*GE[3])
Idef <- Fdi*(1-AE[3])
Jdef <- Fij*(1-AE[4])

FM <- matrix(c(
  0, 0, Fdi, 0,
  0, 0, 0, 0,
  Imort, Idef, 0, Fij,
  Jdef, Jmort, 0, 0),
  nrow = 4, ncol = 4, byrow = T)
rownames(FM) <- fwnames
colnames(FM) <- fwnames

# model
frac <- matrix(c(
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 1, 0, 0,
  1, 0, 0, 0),
  nrow = 4, ncol = 4, byrow = T)
rownames(frac) <- fwnames
colnames(frac) <- fwnames

dead <- list(names = c("D","R"), frac = frac)
model <- list(
  type = "EF", FM = FM, BM = BM, AE = AE, GE = GE,
  dead = dead
)
# Expected answer
JM <- matrix(c(0,
               0,
               (FM[3,1] - FM[1,3] + FM[3,4]*(1-AE[4])) / BM[3],
               (FM[4,1] - FM[1,4]) / BM[4],

               (FM[1,3]*(1-AE[3])) / BM[1],
               0,
               (FM[3,2] - FM[2,3]) / BM[3],
               (FM[4,2] - FM[2,4]) / BM[4],

               AE[3] * GE[3] * FM[1,3] / BM[1],
               0,
               0,
               -FM[3,4] / BM[4],

               AE[4] * GE[4] * FM[1,4] / BM[1],
               AE[4] * GE[4] * FM[2,4] / BM[2],
               AE[4] * GE[4] * FM[3,4] / BM[3],
               0
), nrow = 4, ncol = 4)
rownames(JM) <- fwnames
colnames(JM) <- fwnames
