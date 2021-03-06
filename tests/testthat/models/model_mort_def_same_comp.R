### Dead compartments, no externals
# Meiofauna feeds on labile detritus.
# Meiofauna and macrofauna defecate into both labile and refractory detritus.
# Meiofauna and macrofauna both also deposit carbon due to mortality into the labile compartment.
#
#
# LABILE -7-> MEIO -4->  MACRO
#  ^______3______/ <____2___/  # defecation + mortality: 2/3 for meio, 2/4 for macro.
# REFRAC <--4---/ <----2---/   # defecation
#
## Proper data format
fwnames <- c("LABILE", "REFRAC", "MEIO", "MACRO")
FM <- matrix(c(
  0, 0, 7, 0,
  0, 0, 0, 0,
  3, 4, 0, 4,
  2, 2, 0, 0),
  nrow = 4, ncol = 4, byrow = T)
rownames(FM) <- fwnames
colnames(FM) <- fwnames
BM <- c(30, 20, 10, 5) ; names(BM) <- fwnames
AE <- c(NA, NA, 0.3, 0.3) ; names(AE) <- fwnames
GE <- c(NA, NA, 0.3, 0.3) ; names(GE) <- fwnames
DM <- matrix(0, nrow = 4, ncol = 4)
rownames(DM) <- fwnames
colnames(DM) <- fwnames
DM["MEIO", "LABILE"] <- 2/3
DM["MACRO", "LABILE"] <- 1/2
DM["MEIO", "REFRAC"] <- 1
DM["MACRO", "REFRAC"] <- 1
FDM <- FM * DM
dead <- list(names = c("LABILE", "REFRAC"), frac = DM)
model <- list(
  type = "EF", FM = FM, BM = BM, AE = AE, GE = GE, dead = dead
)
# Expected answer
JM <- matrix(c(0,
               0,
               (FM[3,1] - FM[1,3] + FM[3,4]*(1-AE[4])*(FDM[4,1]/(FDM[4,1]+FDM[4,2]))) / BM[3],
               (FM[4,1] - FM[1,4]) / BM[4],

               (FM[1,2] - FM[2,1] + FM[1,3]*(1-AE[3])*(FDM[3,2]/(FDM[3,1]+FDM[3,2]))) / BM[1],
               0,
               (FM[3,2] - FM[2,3] + FM[3,4]*(1-AE[4])*(FDM[4,2]/(FDM[4,1]+FDM[4,2]))) / BM[3],
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
