# Species deposit material in different dead compartments
#
#   <----\        def + mort
#  d --> i --> j
#  r <-------/    def+ mort

fwnames <- c("D", "R", "I", "J")
FM <- matrix(c(0, 0, 0, 0,  0, 0, 0, 0,  20, 0, 0, 0,  0, 0, 10, 0), nrow = 4, ncol = 4)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 30, 20, 10) ; names(BM) <- fwnames
AE <- c(NA, NA, 0.8, 0.7) ; names(AE) <- fwnames
GE <- c(NA, NA, 0.8, 0.7) ; names(GE) <- fwnames
# Detritus feedback
defecation <- colSums(FM[,3:4])*(1-AE[3:4])
mortalities <- colSums(FM[,3:4])*AE[3:4]*GE[3:4] - rowSums(FM[3:4,])
Dldet <- c(0.5, 0.5) ; Mldet <- c(0.2, 0.3)
Drdet <- 1-Dldet ; Mrdet <- 1-Mldet
FM[3:4,1] <- defecation*Dldet + mortalities*Mldet
FM[3:4,2] <- defecation*Drdet + mortalities*Mrdet
frac <- matrix(
  c(0, 0, defecation*Dldet/FM[3:4,1],
    0, 0, defecation*Drdet/FM[3:4,2],
    rep(0,8)),
  ncol = 4, nrow = 4)
rownames(frac) <- fwnames; colnames(frac) <- fwnames
dead <- list(names = c("D","R"), frac = frac)
# Combine to model
model <- list(
  type = "EF", FM = FM, BM = BM, AE = AE, GE = GE,
  dead = dead
  )
# Expected answer
JM <- matrix(c(0,
               0,
               (FM[3,1] - FM[1,3] + FM[3,4]*(1-AE[4])*Dldet[2]) / BM[3],
               (FM[4,1] - FM[1,4]) / BM[4],

               FM[1,3]*(1-AE[3])*Dldet[1] / BM[1],
               0,
               (FM[3,2] - FM[2,3] + FM[3,4]*(1-AE[4])*Drdet[2]) / BM[3],
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
