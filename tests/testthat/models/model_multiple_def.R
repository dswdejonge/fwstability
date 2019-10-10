### Dead compartments, no externals
# Meiofauna feeds on labile detritus.
# Meiofauna and macrofauna defecate into both labile and refractory detritus.
#
#
# LABILE --> MEIO -->  MACRO
#  ^____________/ <_______/  # defecation 50% + mortality
# REFRAC <-----/ <-------/   # defecation 50%
#
## Proper data format
fwnames <- c("LABILE", "REFRAC", "MEIO", "MACRO")
BM <- c(30, 20, 10, 5) ; names(BM) <- fwnames
AE <- c(NA, NA, 0.3, 0.3) ; names(AE) <- fwnames
GE <- c(NA, NA, 0.3, 0.3) ; names(GE) <- fwnames
MR <- c(1, 1, 0.8, 0.7) ; names(MR) <- fwnames
# Top down balancing
Fd <- 0.5
meiMort <- MR[3]*BM[3]
macMort <- MR[4]*BM[4]
Fmei_mac <- meiMort/(AE[4]*GE[4])
Flab_mei <- (meiMort+Fmei_mac)/(AE[3]*GE[3])
meiDef <- (1-AE[3])*Flab_mei
macDef <- (1-AE[4])*Fmei_mac
Fmei_lab <- meiMort+meiDef*Fd
Fmei_ref <- meiDef*(1-Fd)
Fmac_lab <- macMort+macDef*Fd
Fmac_ref <- macDef*(1-Fd)
FM <- matrix(c(
  0, 0, Flab_mei, 0,
  0, 0, 0, 0,
  Fmei_lab, Fmei_ref, 0, Fmei_mac,
  Fmac_lab, Fmac_ref, 0, 0),
  nrow = 4, ncol = 4, byrow = T)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
# Model
frac <- matrix(c(
  0, 0, 0, 0,
  0, 0, 0, 0,
  meiDef*Fd/Fmei_lab, meiDef*(1-Fd)/Fmei_ref, 0, 0,
  macDef*Fd/Fmac_lab, macDef*(1-Fd)/Fmac_ref, 0, 0),
  nrow = 4, ncol = 4, byrow = T)
rownames(frac) <- fwnames
colnames(frac) <- fwnames
dead <- list(names = c("LABILE", "REFRAC"), def = c("Def", "Def"), frac = frac)
model <- list(
  type = "EF", FM = FM, BM = BM, AE = AE, GE = GE, dead = dead
)
# Expected answer
JM <- matrix(c(0,
               0,
               (FM[3,1] - FM[1,3] + FM[3,4]*(1-AE[4])*Fd) / BM[3],
               (FM[4,1] - FM[1,4]) / BM[4],

               (FM[1,2] - FM[2,1] + FM[1,3]*(1-AE[3])*(1-Fd)) / BM[1],
               0,
               (FM[3,2] - FM[2,3] + FM[3,4]*(1-AE[4])*(1-Fd)) / BM[3],
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



