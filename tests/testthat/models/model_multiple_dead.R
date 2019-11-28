### Dead compartments, no externals
# Meiofauna feeds on labile detritus.
# Meiofauna and macrofauna deposit carcasses into labile detritus (mortality).
# Meiofauna and macrofauna defecate into refractory detritus.
#
#
# LABILE -7-> MEIO -4->  MACRO
#  ^______3______/ <____2___/  # mortality
# REFRAC <--4---/ <----2---/   # defecation
#
## Proper data format
fwnames <- c("LABILE", "REFRAC", "MEIO", "MACRO")
FM <- matrix(c(
  0, 0, 7, 0,
  0, 0, 0, 0,
  0, 0, 0, 4,
  0, 0, 0, 0),
  nrow = 4, ncol = 4, byrow = T)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10, 5) ; names(BM) <- fwnames
AE <- c(NA, NA, 0.8, 0.3) ; names(AE) <- fwnames
GE <- c(NA, NA, 0.9, 0.3) ; names(GE) <- fwnames
# Detritus feedback
defecation <- colSums(FM[,3:4])*(1-AE[3:4])
mortalities <- colSums(FM[,3:4])*AE[3:4]*GE[3:4] - rowSums(FM[3:4,])
meiLab <- 0.2
meiRef <- 1-meiLab
macLab <- 0.3
macRef <- 1-macLab
FM[3:4,1] <- defecation*c(meiLab, macLab) + mortalities*c(0.9, 0.7)
FM[3:4,2] <- defecation*c(meiRef, macRef) + mortalities*c(0.1, 0.3)
# Fraction defecation matrix
frac <- matrix(
  c(0, 0, defecation*c(meiLab, macLab)/FM[3:4,1],
    0, 0, defecation*c(meiRef, macRef)/FM[3:4,2],
    rep(0,8)),
  ncol = 4, nrow = 4)
rownames(frac) <- fwnames; colnames(frac) <- fwnames
dead <- list(names = c("LABILE", "REFRAC"), frac = frac)
model <- list(
  type = "EF", FM = FM, BM = BM, AE = AE, GE = GE, dead = dead
)
# Expected answer
JM <- matrix(c(0,
               0,
               ((FM[3,1] - FM[1,3]) + FM[3,4]*(1-AE[4])*macLab) / BM[3],
               (FM[4,1] - FM[1,4]) / BM[4],

               ((FM[1,2] - FM[2,1]) + FM[1,3]*(1-AE[3])*meiRef) / BM[1],
               0,
               ((FM[3,2] - FM[2,3]) + FM[3,4]*(1-AE[4])*macRef) / BM[3],
               (FM[4,2] - FM[2,4]) / BM[4],

               AE[3] * GE[3] * FM[1,3] / BM[1],
               AE[3] * GE[3] * FM[2,3] / BM[2],
               0,
               -FM[3,4] / BM[4],

               AE[4] * GE[4] * FM[1,4] / BM[1],
               AE[4] * GE[4] * FM[2,4] / BM[2],
               AE[4] * GE[4] * FM[3,4] / BM[3],
               0
), nrow = 4, ncol = 4)
rownames(JM) <- fwnames ; colnames(JM) <- fwnames
